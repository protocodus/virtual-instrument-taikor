#include "ModalRuntime.h"

#include <algorithm>
#include <chrono>
#include <cmath>

namespace taikor::candidates::modal
{
namespace
{
// The one-sided Hertz potential and its derivatives.
struct Hertz
{
    double k = 0.0, alpha = 1.5;
    [[nodiscard]] double phi (double x) const { return x > 0.0 ? k / (alpha + 1.0) * std::pow (x, alpha + 1.0) : 0.0; }
    [[nodiscard]] double force (double x) const { return x > 0.0 ? k * std::pow (x, alpha) : 0.0; }
    [[nodiscard]] double stiffness (double x) const { return x > 0.0 ? k * alpha * std::pow (x, alpha - 1.0) : 0.0; }
};

// Discrete gradient (Phi(eta + r) - Phi(eta)) / r; within a relative 1e-7 the
// midpoint derivative (second-order accurate, no cancellation).
double discreteGradient (const Hertz& h, double eta, double r)
{
    const double scale = std::max (std::abs (eta), std::abs (eta + r));
    if (std::abs (r) <= 1.0e-7 * scale || r == 0.0)
        return h.force (eta + 0.5 * r);
    return (h.phi (eta + r) - h.phi (eta)) / r;
}

double discreteGradientSlope (const Hertz& h, double eta, double r)
{
    const double scale = std::max (std::abs (eta), std::abs (eta + r));
    if (std::abs (r) <= 1.0e-7 * scale || r == 0.0)
        return 0.5 * h.stiffness (eta + 0.5 * r);
    return (h.force (eta + r) * r - (h.phi (eta + r) - h.phi (eta))) / (r * r);
}

// Solves r = U - W F(r), F(r) = discrete gradient, for the gap increment r.
// f(r) = r - U + W F(r) is increasing (Phi convex), so the root is unique; it
// is bracketed by [U - W F(U), U] and found by Newton safeguarded by bisection.
double solveContact (const Hertz& h, double eta, double U, double W, double& r)
{
    if (eta <= 0.0 && eta + U <= 0.0)
    {
        r = U;
        return 0.0;
    }
    const double fU = W * discreteGradient (h, eta, U);
    if (fU <= 0.0)
    {
        r = U;
        return 0.0;
    }
    double hi = U, lo = U - fU;
    double x = lo;
    const auto f = [&] (double t) { return t - U + W * discreteGradient (h, eta, t); };
    double flo = f (lo);
    if (flo >= 0.0)
    {
        r = lo;
        return discreteGradient (h, eta, lo);
    }
    const double tol = 1.0e-15 * (std::abs (eta) + std::abs (U) + 1.0e-30);
    for (int iteration = 0; iteration < 200; ++iteration)
    {
        const double fx = f (x);
        if (fx == 0.0)
            break;
        if (fx < 0.0)
            lo = x;
        else
            hi = x;
        const double slope = 1.0 + W * discreteGradientSlope (h, eta, x);
        double next = x - fx / slope;
        if (! (next > lo && next < hi))
            next = 0.5 * (lo + hi);
        if (std::abs (next - x) <= tol || hi - lo <= tol)
        {
            x = next;
            break;
        }
        x = next;
    }
    r = x;
    return discreteGradient (h, eta, x);
}
} // namespace

Stepper::Stepper (Bank bank, double sampleRate, const RuntimeOptions& o) : data (std::move (bank)), options (o)
{
    n = static_cast<int> (data.modes.size());
    const auto un = static_cast<std::size_t> (n);
    for (auto* v : { &y, &q, &c, &xt, &u, &gB, &gR, &h1L, &h2L, &h1R, &h2R, &v1, &v2, &w, &ggB, &ggR, &sv })
        v->assign (un, 0.0);
    for (int k = 0; k < n; ++k)
        c[static_cast<std::size_t> (k)] = data.modes[static_cast<std::size_t> (k)].contact;
    clustered.assign (un, 0);
    for (const auto& cl : data.clusters)
        for (int k : cl.members)
            clustered[static_cast<std::size_t> (k)] = 1;
    setRate (base, 1.0 / sampleRate);
    setRate (over, 1.0 / (sampleRate * std::max (1, options.oversample)));
    for (int k = 0; k < n; ++k)
    {
        const auto uk = static_cast<std::size_t> (k);
        const auto& m = data.modes[uk];
        const double mr = base.muR[uk], mi = base.muI[uk];
        // p = Re{H A}, A = y + j (y mu_r - q)/mu_i.
        h1L[uk] = m.mic[0].real() - m.mic[0].imag() * mr / mi;
        h2L[uk] = m.mic[0].imag() / mi;
        h1R[uk] = m.mic[1].real() - m.mic[1].imag() * mr / mi;
        h2R[uk] = m.mic[1].imag() / mi;
        // dy/dt = Re{s A} = -alpha y - w (y mu_r - q)/mu_i, weighted by c.
        v1[uk] = c[uk] * (-m.alpha - m.omega * mr / mi);
        v2[uk] = c[uk] * m.omega / mi;
    }
    std::size_t largest = 0;
    for (const auto& head : data.berger)
        for (const auto& g : head.groups)
            largest = std::max (largest, g.rows.size());
    headQ.assign (largest, 0.0);
    headDs.assign (largest, 0.0);
    current = &base;
    // The feedthrough kernel is binned at the substep this stepper uses.
    const std::size_t span = data.feedthrough[0].size();
    if (span > 0 && std::abs (data.feedthroughBin - over.dt) <= 1.0e-9 * over.dt)
    {
        useFeedthrough = true;
        forceHistory.assign (span, 0.0);
        forceAge = span;
    }
}

void Stepper::pushForce (double force)
{
    if (! useFeedthrough)
        return;
    historyHead = historyHead + 1 == forceHistory.size() ? 0 : historyHead + 1;
    forceHistory[historyHead] = force;
    forceAge = force != 0.0 ? 0 : std::min (forceAge + 1, forceHistory.size());
}

// p -= sum_l F_l W_l, F_l the force of the substep l substeps before the latest
// (ModalFeedthrough.cpp).
void Stepper::subtractFeedthrough (double& left, double& right) const
{
    if (! useFeedthrough || forceAge >= forceHistory.size())
        return;
    const double* w0 = data.feedthrough[0].data();
    const double* w1 = data.feedthrough[1].data();
    const std::size_t span = forceHistory.size();
    double sl = 0.0, sr = 0.0;
    std::size_t lag = 0;
    for (std::size_t i = historyHead + 1; i-- > 0; ++lag)
    {
        sl += forceHistory[i] * w0[lag];
        sr += forceHistory[i] * w1[lag];
    }
    for (std::size_t i = span; i-- > historyHead + 1; ++lag)
    {
        sl += forceHistory[i] * w0[lag];
        sr += forceHistory[i] * w1[lag];
    }
    left -= sl;
    right -= sr;
}

void Stepper::setRate (Rate& r, double dt)
{
    const auto un = static_cast<std::size_t> (n);
    r.dt = dt;
    r.xi = 0.5 * dt * dt;
    r.xiStick = data.tipMass > 0.0 ? 0.5 * dt * dt / data.tipMass : 0.0;
    for (auto* v : { &r.a, &r.d, &r.bOverXi, &r.muR, &r.muI })
        v->assign (un, 0.0);
    r.cec = 0.0;
    for (int k = 0; k < n; ++k)
    {
        const auto uk = static_cast<std::size_t> (k);
        const auto& m = data.modes[uk];
        const double R = std::exp (-m.alpha * dt), O = std::cos (m.omega * dt);
        const double den = 1.0 + 2.0 * R * O + R * R;
        const double a = (1.0 - 2.0 * R * O + R * R) / den;
        const double b = 2.0 * (1.0 - R * R) / den;
        r.a[uk] = a;
        r.d[uk] = 1.0 / (1.0 + a + b);
        r.bOverXi[uk] = b / r.xi;
        const std::complex<double> zeta = std::polar (R, m.omega * dt);
        const std::complex<double> mu = (zeta - 1.0) / (zeta + 1.0);
        r.muR[uk] = mu.real();
        r.muI[uk] = mu.imag();
    }
    // Clusters: B = S (dt A) S with S = diag(sqrt(b*_k / (alpha_k dt))), so the
    // diagonal is each mode's own pre-warped b*_k and B stays positive
    // semi-definite (a congruence of A); G = (I + diag(a) + B)^{-1}.
    r.clusterG.clear();
    r.clusterBOverXi.clear();
    for (const auto& cl : data.clusters)
    {
        const std::size_t p = cl.members.size();
        std::vector<double> scale (p), m (p * p), b (p * p);
        for (std::size_t i = 0; i < p; ++i)
        {
            const auto k = static_cast<std::size_t> (cl.members[i]);
            const double alpha = data.modes[k].alpha;
            const double bk = (1.0 / r.d[k]) - 1.0 - r.a[k];
            scale[i] = alpha > 0.0 ? std::sqrt (bk / (alpha * dt)) : 1.0;
        }
        for (std::size_t i = 0; i < p; ++i)
            for (std::size_t j = 0; j < p; ++j)
            {
                b[i * p + j] = scale[i] * scale[j] * dt * cl.decay[i * p + j];
                m[i * p + j] = b[i * p + j] + (i == j ? 1.0 + r.a[static_cast<std::size_t> (cl.members[i])] : 0.0);
            }
        // Invert the small symmetric positive definite matrix by Gauss-Jordan.
        std::vector<double> inv (p * p, 0.0);
        for (std::size_t i = 0; i < p; ++i)
            inv[i * p + i] = 1.0;
        for (std::size_t col = 0; col < p; ++col)
        {
            std::size_t pivot = col;
            for (std::size_t row = col + 1; row < p; ++row)
                if (std::abs (m[row * p + col]) > std::abs (m[pivot * p + col]))
                    pivot = row;
            for (std::size_t j = 0; j < p; ++j)
            {
                std::swap (m[col * p + j], m[pivot * p + j]);
                std::swap (inv[col * p + j], inv[pivot * p + j]);
            }
            const double diag = m[col * p + col];
            for (std::size_t j = 0; j < p; ++j)
            {
                m[col * p + j] /= diag;
                inv[col * p + j] /= diag;
            }
            for (std::size_t row = 0; row < p; ++row)
                if (row != col)
                {
                    const double factor = m[row * p + col];
                    for (std::size_t j = 0; j < p; ++j)
                    {
                        m[row * p + j] -= factor * m[col * p + j];
                        inv[row * p + j] -= factor * inv[col * p + j];
                    }
                }
        }
        for (auto& value : b)
            value /= r.xi;
        r.clusterG.push_back (std::move (inv));
        r.clusterBOverXi.push_back (std::move (b));
    }
    r.gc.assign (un, 0.0);
    applyG (r, c.data(), r.gc.data());
    r.cec = 0.0;
    for (int k = 0; k < n; ++k)
        r.cec += r.xi * c[static_cast<std::size_t> (k)] * r.gc[static_cast<std::size_t> (k)];
}

void Stepper::applyG (const Rate& r, const double* in, double* out) const
{
    for (int k = 0; k < n; ++k)
        out[k] = r.d[static_cast<std::size_t> (k)] * in[k];
    for (std::size_t ci = 0; ci < data.clusters.size(); ++ci)
    {
        const auto& members = data.clusters[ci].members;
        const auto& g = r.clusterG[ci];
        const std::size_t p = members.size();
        for (std::size_t i = 0; i < p; ++i)
        {
            double sum = 0.0;
            for (std::size_t j = 0; j < p; ++j)
                sum += g[i * p + j] * in[members[j]];
            out[members[i]] = sum;
        }
    }
}

void Stepper::berger (const std::vector<double>& x, double& sB, double& sR)
{
    double* g[2] = { gB.data(), gR.data() };
    double* s[2] = { &sB, &sR };
    for (int h = 0; h < 2; ++h)
    {
        std::fill (g[h], g[h] + n, 0.0);
        *s[h] = 0.0;
        const auto& head = data.berger[static_cast<std::size_t> (h)];
        if (head.beta <= 0.0)
            continue;
        double total = 0.0;
        for (const auto& group : head.groups)
        {
            const std::size_t rows = group.rows.size();
            for (std::size_t i = 0; i < rows; ++i)
            {
                double v = 0.0;
                for (const auto& [k, coef] : group.rows[i])
                    v += coef * x[static_cast<std::size_t> (k)];
                headQ[i] = v;
                headDs[i] = 2.0 * group.diagonal[i] * v;
                total += group.diagonal[i] * v * v;
            }
            for (const auto& l : group.lowRank)
            {
                const double* lp = l.data();
                double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
                std::size_t i = 0;
                for (; i + 4 <= rows; i += 4)
                {
                    t0 += lp[i] * headQ[i];
                    t1 += lp[i + 1] * headQ[i + 1];
                    t2 += lp[i + 2] * headQ[i + 2];
                    t3 += lp[i + 3] * headQ[i + 3];
                }
                for (; i < rows; ++i)
                    t0 += lp[i] * headQ[i];
                const double t = (t0 + t1) + (t2 + t3);
                total -= t * t;
                const double twoT = 2.0 * t;
                for (std::size_t j = 0; j < rows; ++j)
                    headDs[j] -= twoT * lp[j];
            }
            for (std::size_t i = 0; i < rows; ++i)
                for (const auto& [k, coef] : group.rows[i])
                    g[h][k] += coef * headDs[i];
        }
        const double scale = std::sqrt (2.0 * head.beta);
        for (int k = 0; k < n; ++k)
            g[h][k] *= scale;
        *s[h] = total;
    }
}

double Stepper::step (const Rate& r)
{
    const double xi = r.xi;
    const bool useBerger = options.berger;
    for (int k = 0; k < n; ++k)
        xt[static_cast<std::size_t> (k)] = y[static_cast<std::size_t> (k)] + q[static_cast<std::size_t> (k)];
    stepS[0] = stepS[1] = 0.0;
    if (useBerger)
        berger (xt, stepS[0], stepS[1]);
    // u = 2 G (q - a y); G g_B, G g_R; G c is precomputed per rate.
    for (int k = 0; k < n; ++k)
        w[static_cast<std::size_t> (k)] = q[static_cast<std::size_t> (k)] - r.a[static_cast<std::size_t> (k)] * y[static_cast<std::size_t> (k)];
    applyG (r, w.data(), u.data());
    if (useBerger)
    {
        applyG (r, gB.data(), ggB.data());
        applyG (r, gR.data(), ggR.data());
    }
    double gu0 = 0.0, gu1 = 0.0, e00 = 0.0, e01 = 0.0, e11 = 0.0, gc0 = 0.0, gc1 = 0.0, cu = 0.0;
    for (int k = 0; k < n; ++k)
    {
        const auto uk = static_cast<std::size_t> (k);
        const double value = 2.0 * u[uk];
        u[uk] = value;
        cu += c[uk] * value;
        if (useBerger)
        {
            const double b = gB[uk], rr = gR[uk];
            gu0 += b * value;
            gu1 += rr * value;
            e00 += xi * b * ggB[uk];
            e01 += xi * b * ggR[uk];
            e11 += xi * rr * ggR[uk];
            gc0 += xi * b * r.gc[uk];
            gc1 += xi * rr * r.gc[uk];
        }
    }
    // (I + E/2) psibar = psi + Gu/2 - F Gc/2  ->  psibar = P0 + F P1.
    const double a00 = 1.0 + 0.5 * e00, a01 = 0.5 * e01, a11 = 1.0 + 0.5 * e11;
    const double det = a00 * a11 - a01 * a01;
    const double r0 = psi[0] + 0.5 * gu0, r1 = psi[1] + 0.5 * gu1;
    const double p00 = (a11 * r0 - a01 * r1) / det, p01 = (a00 * r1 - a01 * r0) / det;
    const double t0 = -0.5 * gc0, t1 = -0.5 * gc1;
    const double p10 = (a11 * t0 - a01 * t1) / det, p11 = (a00 * t1 - a01 * t0) / det;
    const double U = cu - (p00 * gc0 + p01 * gc1) - 2.0 * qs;
    const double W = r.cec + (p10 * gc0 + p11 * gc1) + r.xiStick;
    double force = 0.0, increment = U;
    if (struck && data.tipMass > 0.0)
        force = solveContact (Hertz { data.hertzStiffness, data.hertzExponent }, eta, U, W, increment);
    const double psibar0 = p00 + force * p10, psibar1 = p01 + force * p11;
    double wbar = 0.0, dissipation = 0.0;
    for (int k = 0; k < n; ++k)
    {
        const auto uk = static_cast<std::size_t> (k);
        double s = u[uk] - xi * force * r.gc[uk];
        if (useBerger)
            s -= xi * (psibar0 * ggB[uk] + psibar1 * ggR[uk]);
        sv[uk] = s;
        if (! clustered[uk])
            dissipation += r.bOverXi[uk] * s * s;
        y[uk] += s;
        q[uk] = s - q[uk];
        wbar += c[uk] * y[uk];
    }
    for (std::size_t ci = 0; ci < data.clusters.size(); ++ci)
    {
        const auto& members = data.clusters[ci].members;
        const auto& b = r.clusterBOverXi[ci];
        const std::size_t p = members.size();
        for (std::size_t i = 0; i < p; ++i)
            for (std::size_t j = 0; j < p; ++j)
                dissipation += sv[static_cast<std::size_t> (members[i])] * b[i * p + j] * sv[static_cast<std::size_t> (members[j])];
    }
    if (useBerger)
    {
        psi[0] = 2.0 * psibar0 - psi[0];
        psi[1] = 2.0 * psibar1 - psi[1];
    }
    if (data.tipMass > 0.0 && struck)
    {
        const double sz = 2.0 * qs + r.xiStick * force;
        z += sz;
        qs = sz - qs;
    }
    eta = wbar - z;
    stepDissipation = dissipation;
    return force;
}

void Stepper::convert (const Rate& from, const Rate& to)
{
    for (int k = 0; k < n; ++k)
    {
        const auto uk = static_cast<std::size_t> (k);
        const double ai = (y[uk] * from.muR[uk] - q[uk]) / from.muI[uk];
        q[uk] = y[uk] * to.muR[uk] - ai * to.muI[uk];
    }
    qs *= to.dt / from.dt;
    current = &to;
}

double Stepper::totalEnergy() const
{
    const Rate& r = *current;
    double e = 0.0;
    for (int k = 0; k < n; ++k)
    {
        const auto uk = static_cast<std::size_t> (k);
        e += (q[uk] * q[uk] + r.a[uk] * y[uk] * y[uk]) / r.xi;
    }
    e += 0.5 * (psi[0] * psi[0] + psi[1] * psi[1]);
    if (data.tipMass > 0.0 && r.xiStick > 0.0)
        e += qs * qs / r.xiStick;
    e += Hertz { data.hertzStiffness, data.hertzExponent }.phi (eta);
    return e;
}

double Stepper::stickVelocity() const { return 2.0 * qs / current->dt; }

std::array<double, 2> Stepper::freeDecayIntegral() const
{
    const Rate& r = *current;
    std::array<double, 2> out { 0.0, 0.0 };
    for (int k = 0; k < n; ++k)
    {
        const auto uk = static_cast<std::size_t> (k);
        const auto& m = data.modes[uk];
        const std::complex<double> amplitude (y[uk], (y[uk] * r.muR[uk] - q[uk]) / r.muI[uk]);
        const std::complex<double> s (-m.alpha, m.omega);
        for (std::size_t j = 0; j < 2; ++j)
            out[j] -= (m.mic[j] * amplitude / s).real();
    }
    return out;
}

void Stepper::setModeState (int k, double displacement, double velocity)
{
    const auto uk = static_cast<std::size_t> (k);
    const auto& m = data.modes[uk];
    // Analytic amplitude A = y + j A_i with dy/dt = -alpha y - w A_i; q = y mu_r - A_i mu_i.
    const double ai = -(velocity + m.alpha * displacement) / m.omega;
    y[uk] = displacement;
    q[uk] = displacement * current->muR[uk] - ai * current->muI[uk];
}

void Stepper::strike()
{
    double wbar = 0.0;
    for (int k = 0; k < n; ++k)
        wbar += c[static_cast<std::size_t> (k)] * y[static_cast<std::size_t> (k)];
    // On the head surface, moving into the drum (-z) at the impact speed.
    z = wbar;
    eta = 0.0;
    qs = -0.5 * current->dt * data.impactSpeed;
    struck = true;
    time = 0.0;
    contactReport = ContactReport {};
    peakPsi[0] = peakPsi[1] = 0.0;
    inEnergyWindow = true;
    energyStarted = false;
    inContact = false;
    lastSampleContact = false;
    etaRate = data.impactSpeed;
}

void Stepper::process (float* left, float* right, int samples)
{
    thread_local std::vector<double> l, r;
    l.resize (static_cast<std::size_t> (samples));
    r.resize (static_cast<std::size_t> (samples));
    processDouble (l.data(), r.data(), samples);
    for (int i = 0; i < samples; ++i)
    {
        left[i] = static_cast<float> (l[static_cast<std::size_t> (i)]);
        right[i] = static_cast<float> (r[static_cast<std::size_t> (i)]);
    }
}

void Stepper::processDouble (double* left, double* right, int samples)
{
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < samples; ++i)
    {
        bool oversampled = false;
        if (struck && data.tipMass > 0.0)
        {
            if (options.forceRate == 2)
                oversampled = true;
            else if (options.forceRate == 0)
                oversampled = lastSampleContact || eta > 0.0 || eta + 3.0 * base.dt * std::max (etaRate, 0.0) > 0.0;
        }
        bool contact = false;
        const auto track = [&] (double force, double dt, double subTime)
        {
            if (force > 0.0)
            {
                contact = true;
                if (! inContact)
                    ++contactReport.episodes;
                inContact = true;
                if (contactReport.firstContact < 0.0)
                    contactReport.firstContact = subTime;
                contactReport.lastContact = subTime;
                contactReport.peakForce = std::max (contactReport.peakForce, force);
                contactReport.impulse += force * dt;
            }
            else
                inContact = false;
            contactReport.maxPenetration = std::max (contactReport.maxPenetration, eta);
            for (int h = 0; h < 2; ++h)
                contactReport.peakTensionRise[h] = std::max (contactReport.peakTensionRise[h],
                                                             data.berger[static_cast<std::size_t> (h)].tensionPerS * stepS[h]);
        };
        if (oversampled)
        {
            if (current != &over)
                convert (base, over);
            if (inEnergyWindow && ! energyStarted)
            {
                contactReport.energyStart = totalEnergy();
                energyStarted = true;
            }
            for (int sub = 0; sub < std::max (1, options.oversample); ++sub)
            {
                const double force = step (over);
                if (inEnergyWindow)
                    contactReport.dissipated += stepDissipation;
                track (force, over.dt, time + (sub + 1) * over.dt);
                pushForce (force);
            }
            if (inEnergyWindow)
                contactReport.energyEnd = totalEnergy();
            convert (over, base);
            ++contactReport.oversampledSamples;
        }
        else
        {
            if (struck && inEnergyWindow && energyStarted)
            {
                inEnergyWindow = false;
                contactReport.energyWindowClosed = true;
            }
            const double force = step (base);
            if (force > 0.0)
                ++contactReport.baseRateContactSteps;
            track (force, base.dt, time + base.dt);
            for (int sub = 0; sub < std::max (1, options.oversample); ++sub)
                pushForce (force);
        }
        lastSampleContact = contact;
        double pl = 0.0, pr = 0.0, rate = 0.0;
        for (int k = 0; k < n; ++k)
        {
            const auto uk = static_cast<std::size_t> (k);
            pl += h1L[uk] * y[uk] + h2L[uk] * q[uk];
            pr += h1R[uk] * y[uk] + h2R[uk] * q[uk];
            rate += v1[uk] * y[uk] + v2[uk] * q[uk];
        }
        etaRate = rate - 2.0 * qs / base.dt;
        subtractFeedthrough (pl, pr);
        if (left != nullptr)
            left[i] = pl;
        if (right != nullptr)
            right[i] = pr;
        if (struck)
        {
            time += base.dt;
            // Occasional check of the SAV variable against sqrt(2 beta) S(y).
            if (options.berger && (++sampleCount % 64) == 0)
            {
                double sB = 0.0, sR = 0.0;
                berger (y, sB, sR);
                for (int h = 0; h < 2; ++h)
                {
                    const double exact = std::sqrt (2.0 * data.berger[static_cast<std::size_t> (h)].beta) * (h == 0 ? sB : sR);
                    peakPsi[h] = std::max ({ peakPsi[h], std::abs (psi[h]), std::abs (exact) });
                    if (peakPsi[h] > 0.0)
                        contactReport.maxPsiDrift = std::max (contactReport.maxPsiDrift, std::abs (psi[h] - exact) / peakPsi[h]);
                }
            }
        }
    }
    contactReport.exitVelocity = stickVelocity();
    renderSeconds += std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
}

} // namespace taikor::candidates::modal
