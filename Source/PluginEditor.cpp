#include "PluginEditor.h"

#include <BinaryData.h>

#include <cmath>
#include <initializer_list>
#include <limits>

namespace
{
// Design size and the range the editor may be resized through. The plug-in
// tests read the same numbers, so a layout change that breaks the contract
// fails the suite rather than the host.
constexpr int designWidth = 1280;
constexpr int designHeight = 880;
constexpr int minimumWidth = 1024;
constexpr int minimumHeight = 704;
constexpr int maximumWidth = 1472;
constexpr int maximumHeight = 1012;

constexpr float pi = 3.14159265358979f;

// Indigo lacquer frames opaque washi panels; artwork never competes with a
// parameter label. Vermilion is reserved for selection, strikes and focus.
const juce::Colour backgroundTop { 0xff21383d };
const juce::Colour backgroundBottom { 0xff152b31 };
const juce::Colour panelColour { 0xfff1e9d8 };
const juce::Colour panelRaised { 0xfffbf5e9 };
const juce::Colour panelEdge { 0xffb9ad95 };
const juce::Colour hideColour { 0xffe6d4ad };
const juce::Colour washiColour { 0xfffaf1df };
const juce::Colour accentColour { 0xffa53e2e };
const juce::Colour accentDim { 0xff813326 };
const juce::Colour brassColour { 0xff806334 };
const juce::Colour textColour { 0xff25383b };
const juce::Colour mutedText { 0xff5b625b };

juce::Font controlFont (float height, bool bold = false)
{
    return juce::Font (juce::FontOptions (height).withStyle (bold ? "Bold" : "Regular"));
}

float editorScale (const juce::Component& component)
{
    if (const auto* editor = component.findParentComponentOfClass<juce::AudioProcessorEditor>())
        return static_cast<float> (editor->getWidth()) / designWidth;
    return 1.0f;
}

juce::Font fitFont (juce::Font font, const juce::String& text, float width)
{
    const auto textWidth = juce::GlyphArrangement::getStringWidth (font, text);
    if (textWidth > width && width > 0.0f)
        font.setHeight (font.getHeight() * width / textWidth);
    return font;
}

juce::Colour roleColour (TaikorKnob::VisualRole role) noexcept
{
    switch (role)
    {
        case TaikorKnob::VisualRole::Drum:       return brassColour;
        case TaikorKnob::VisualRole::Stroke:     return accentColour;
        case TaikorKnob::VisualRole::Microphone: return juce::Colour { 0xff375f60 };
        case TaikorKnob::VisualRole::Master:     return textColour;
    }
    return accentColour;
}

juce::Font displayFont (float height, int style = juce::Font::plain)
{
    return juce::Font (juce::FontOptions (
        juce::Font::getDefaultSerifFontName(), height, style));
}

// Scientific pitch notation puts middle C (MIDI 60) at C4, which makes this
// instrument's reference note - MIDI 48 - C3. Both the octave strip and the
// pads derive their labels from that one constant rather than each carrying
// their own idea of it, because they disagreed by an octave when they did:
// JUCE's note-name helper takes the octave number to give middle C, and it was
// being handed the reference note's octave instead.
constexpr int referenceOctaveNumber = taikor::referenceNote / 12 - 1;
constexpr int octaveNumberForMiddleC =
    referenceOctaveNumber + (60 - taikor::referenceNote) / 12;

static_assert (referenceOctaveNumber == 3,
               "the documented mapping calls MIDI 48 C3");
static_assert (octaveNumberForMiddleC == 4,
               "middle C must be C4 if MIDI 48 is C3");

juce::String octaveName (int octaveOffset)
{
    return "C" + juce::String (referenceOctaveNumber + octaveOffset);
}

// The drum an octave plays, read straight off the engine's own table rather
// than restated here. Each octave is a different instrument of the family, so
// the octave row names instruments: there is nothing for the editor to have its
// own opinion about.
juce::String octaveDescription (int octaveOffset)
{
    const auto& name = taikor::getDrumDescription (octaveOffset).displayName;
    return juce::String (name.data(), name.size());
}

juce::String octaveSummary (int octaveOffset)
{
    const auto& summary = taikor::getDrumDescription (octaveOffset).summary;
    return juce::String (summary.data(), summary.size());
}
} // namespace

// ---------------------------------------------------------------------------
// Look and feel
// ---------------------------------------------------------------------------

TaikorLookAndFeel::TaikorLookAndFeel()
{
    setColour (juce::Slider::rotarySliderFillColourId, accentColour);
    setColour (juce::Slider::rotarySliderOutlineColourId, panelEdge);
    setColour (juce::Slider::textBoxTextColourId, textColour);
    setColour (juce::Slider::textBoxOutlineColourId, juce::Colours::transparentBlack);
    setColour (juce::Slider::textBoxBackgroundColourId, juce::Colours::transparentBlack);
    setColour (juce::Label::textColourId, textColour);
    setColour (juce::Label::outlineColourId, juce::Colours::transparentBlack);
    setColour (juce::Label::backgroundColourId, juce::Colours::transparentBlack);
    setColour (juce::TooltipWindow::backgroundColourId, panelColour);
    setColour (juce::TooltipWindow::textColourId, textColour);
    setColour (juce::TooltipWindow::outlineColourId, panelEdge);
}

void TaikorLookAndFeel::drawRotarySlider (juce::Graphics& g, int x, int y, int width,
                                          int height, float sliderPos,
                                          float rotaryStartAngle, float rotaryEndAngle,
                                          juce::Slider& slider)
{
    const auto scale = editorScale (slider);
    const auto bounds = juce::Rectangle<int> (x, y, width, height).toFloat().reduced (3.0f * scale);
    const auto radius = juce::jmin (29.0f, juce::jmin (bounds.getWidth(), bounds.getHeight()) * 0.5f);
    const auto centre = bounds.getCentre();
    const auto angle = rotaryStartAngle + sliderPos * (rotaryEndAngle - rotaryStartAngle);
    const auto fill = slider.findColour (juce::Slider::rotarySliderFillColourId);
    const auto opacity = slider.isEnabled() ? 1.0f : 0.42f;

    // Engraved divisions remain visible independently of the value arc.
    g.setColour (mutedText.withAlpha (0.55f * opacity));
    for (int tick = 0; tick <= 10; ++tick)
    {
        const auto a = rotaryStartAngle + static_cast<float> (tick) / 10.0f
                                            * (rotaryEndAngle - rotaryStartAngle);
        const auto outer = radius + 2.0f * scale;
        const auto inner = outer - (tick % 5 == 0 ? 3.0f : 1.5f) * scale;
        g.drawLine (centre.x + std::sin (a) * inner, centre.y - std::cos (a) * inner,
                    centre.x + std::sin (a) * outer, centre.y - std::cos (a) * outer, 0.8f);
    }
    juce::Path track;
    track.addCentredArc (centre.x, centre.y, radius, radius, 0.0f,
                         rotaryStartAngle, rotaryEndAngle, true);
    g.setColour (panelEdge.withAlpha (opacity));
    g.strokePath (track, juce::PathStrokeType (2.0f));
    if (angle > rotaryStartAngle)
    {
        juce::Path value;
        value.addCentredArc (centre.x, centre.y, radius, radius, 0.0f,
                             rotaryStartAngle, angle, true);
        g.setColour (fill.withAlpha (opacity));
        g.strokePath (value, juce::PathStrokeType (2.2f));
    }

    const auto disc = juce::Rectangle<float> (radius * 1.66f, radius * 1.66f).withCentre (centre);
    g.setColour (textColour.withAlpha (0.15f * opacity));
    g.fillEllipse (disc.translated (0.0f, 2.5f).expanded (1.0f));
    g.setGradientFill (juce::ColourGradient (textColour.brighter (0.15f), disc.getTopLeft(),
                                            backgroundBottom, disc.getBottomRight(), false));
    g.fillEllipse (disc);
    g.setColour (brassColour.withAlpha (opacity));
    g.drawEllipse (disc, 1.0f);
    g.setColour (washiColour.withAlpha (0.15f * opacity));
    g.drawEllipse (disc.reduced (2.0f), 0.7f);

    juce::Path pointer;
    pointer.addRoundedRectangle (-1.3f, -radius * 0.65f, 2.6f, radius * 0.47f, 1.0f);
    g.setColour (washiColour.withAlpha (opacity));
    g.fillPath (pointer, juce::AffineTransform::rotation (angle).translated (centre.x, centre.y));
    if (slider.hasKeyboardFocus (true) || slider.isMouseOverOrDragging())
    {
        g.setColour (accentColour.withAlpha (slider.hasKeyboardFocus (true) ? 1.0f : 0.45f));
        g.drawEllipse (disc.expanded (2.0f), 1.5f);
    }
}

void TaikorLookAndFeel::drawButtonBackground (juce::Graphics& g, juce::Button& button,
                                              const juce::Colour& backgroundColour,
                                              bool isHighlighted, bool isDown)
{
    const auto bounds = button.getLocalBounds().toFloat().reduced (1.0f);
    auto fill = backgroundColour;
    if (isDown)
        fill = fill.darker (0.12f);
    else if (isHighlighted)
        fill = fill.interpolatedWith (button.getToggleState() ? washiColour : panelEdge, 0.18f);
    g.setColour (fill);
    g.fillRoundedRectangle (bounds, 3.0f);
    g.setColour (button.getToggleState() ? accentDim : panelEdge);
    g.drawRoundedRectangle (bounds, 3.0f, 0.8f);
    if (button.hasKeyboardFocus (true))
    {
        g.setColour (button.getToggleState() ? washiColour : accentColour);
        g.drawRoundedRectangle (bounds.reduced (2.5f), 2.0f, 1.5f);
    }
}

juce::Font TaikorLookAndFeel::getTextButtonFont (juce::TextButton& button, int buttonHeight)
{
    auto font = controlFont (juce::jmin (24.0f * editorScale (button),
                                        static_cast<float> (buttonHeight) * 0.82f), true);
    auto widest = button.getButtonText();
    if (const auto* choices = dynamic_cast<TaikorChoiceSwitch*> (button.getParentComponent()))
        for (const auto* child : choices->getChildren())
            if (const auto* choice = dynamic_cast<const juce::TextButton*> (child))
                if (juce::GlyphArrangement::getStringWidth (font, choice->getButtonText())
                    > juce::GlyphArrangement::getStringWidth (font, widest))
                    widest = choice->getButtonText();
    return fitFont (font, widest, static_cast<float> (button.getWidth()) - 8.0f);
}

void TaikorLookAndFeel::drawButtonText (juce::Graphics& g, juce::TextButton& button,
                                        bool, bool)
{
    g.setFont (getTextButtonFont (button, button.getHeight()));
    g.setColour (button.findColour (button.getToggleState()
                                        ? juce::TextButton::textColourOnId
                                        : juce::TextButton::textColourOffId)
                     .withMultipliedAlpha (button.isEnabled() ? 1.0f : 0.5f));
    g.drawText (button.getButtonText(), button.getLocalBounds().reduced (4, 2),
                juce::Justification::centred, false);
}

juce::Label* TaikorLookAndFeel::createSliderTextBox (juce::Slider& slider)
{
    auto* label = juce::LookAndFeel_V4::createSliderTextBox (slider);
    label->setFont (controlFont (24.0f));
    label->setBorderSize (juce::BorderSize<int> (0));
    label->setMinimumHorizontalScale (1.0f);
    label->setColour (juce::Label::textWhenEditingColourId, textColour);
    label->setColour (juce::Label::backgroundWhenEditingColourId, panelRaised);
    label->setColour (juce::Label::outlineWhenEditingColourId, accentColour);
    label->setColour (juce::TextEditor::highlightColourId, textColour);
    label->setColour (juce::TextEditor::highlightedTextColourId, washiColour);
    label->onEditorShow = [label]
    {
        if (auto* editor = label->getCurrentTextEditor())
            editor->applyFontToAllText (fitFont (
                label->getLookAndFeel().getLabelFont (*label), label->getText(),
                static_cast<float> (editor->getWidth()) - 10.0f));
    };
    return label;
}

juce::Font TaikorLookAndFeel::getLabelFont (juce::Label& label)
{
    if (const auto* slider = dynamic_cast<juce::Slider*> (label.getParentComponent()))
        return fitFont (controlFont (24.0f * editorScale (*slider)), label.getText(),
                        static_cast<float> (label.getWidth()) - 2.0f);
    return juce::LookAndFeel_V4::getLabelFont (label);
}

// ---------------------------------------------------------------------------
// Stroke pad
// ---------------------------------------------------------------------------

TaikorPad::TaikorPad (taikor::Articulation articulationToUse, int octaveOffsetToUse)
    : juce::Button (juce::String (
          taikor::getArticulationDisplayName (articulationToUse).data(),
          taikor::getArticulationDisplayName (articulationToUse).size())),
      articulation (articulationToUse),
      octaveOffset (octaveOffsetToUse)
{
    refreshNoteText();
    setWantsKeyboardFocus (true);
}

void TaikorPad::refreshNoteText()
{
    const auto note = taikor::midiNoteFor (articulation, octaveOffset);
    keyText = juce::MidiMessage::getMidiNoteName (note, true, true,
                                                  octaveNumberForMiddleC);
    noteText = keyText + " (" + juce::String (note) + ")";
    drumNameText = octaveName (octaveOffset) + " " + taikor::getDrumDescription (octaveOffset).displayName.data();
}

void TaikorPad::setSelected (bool shouldBeSelected)
{
    if (selected == shouldBeSelected)
        return;
    selected = shouldBeSelected;
    repaint();
}

void TaikorPad::setOctaveOffset (int newOctaveOffset)
{
    if (octaveOffset == newOctaveOffset)
        return;
    octaveOffset = newOctaveOffset;
    refreshNoteText();
    repaint();

    // The pad now plays a different note, so anything holding the old wording
    // has to be told. A reader that had already fetched the help text would
    // otherwise go on announcing the octave that was selected when it asked.
    if (auto* handler = getAccessibilityHandler())
        handler->notifyAccessibilityEvent (juce::AccessibilityEvent::titleChanged);
}

juce::String TaikorPad::accessibleHelpText() const
{
    const auto& metadata = taikor::getArticulationMetadata (articulation);
    return drumNameText + " " + juce::String (metadata.description.data(), metadata.description.size())
         + ". Plays " + noteText;
}

void TaikorPad::triggerFlash()
{
    flashLevel = 1.0f;
    repaint();
}

void TaikorPad::advanceFlash()
{
    if (flashLevel <= 0.0f)
        return;
    flashLevel *= 0.80f;
    if (flashLevel < 0.01f)
        flashLevel = 0.0f;
    repaint();
}

void TaikorPad::paintButton (juce::Graphics& g, bool isMouseOver, bool isButtonDown)
{
    const auto bounds = getLocalBounds().toFloat().reduced (1.0f);
    const auto active = isButtonDown || flashLevel > 0.05f;
    const auto fill = active ? accentColour : (selected ? textColour : panelRaised);
    const auto ink = active || selected ? washiColour : textColour;
    g.setColour (fill);
    g.fillRoundedRectangle (bounds, 3.0f);
    g.setColour (selected ? textColour : panelEdge);
    g.drawRoundedRectangle (bounds, 3.0f, 0.8f);
    if (isMouseOver && ! active)
    {
        g.setColour (accentColour.withAlpha (selected ? 0.22f : 0.09f));
        g.fillRoundedRectangle (bounds, 3.0f);
    }

    const auto mapSize = juce::jlimit (19.0f, 26.0f, bounds.getHeight() * 0.49f);
    const auto map = juce::Rectangle<float> (mapSize, mapSize).withCentre (
        { bounds.getX() + bounds.getWidth() * 0.25f, bounds.getCentreY() });
    g.setColour (ink.withAlpha (0.65f));
    g.drawEllipse (map, 1.0f);
    auto strike = map.getCentre();
    if (articulation == taikor::Articulation::Ka)
        strike.x = map.getRight() - mapSize * 0.13f;
    else if (articulation == taikor::Articulation::DonRim)
        strike.x = map.getRight();
    g.setColour (selected || active ? juce::Colour (0xffdfba80) : accentColour);
    g.fillEllipse (strike.x - 2.1f, strike.y - 2.1f, 4.2f, 4.2f);
    if (articulation == taikor::Articulation::Tsu)
    {
        g.setColour (ink);
        g.drawLine (map.getX() + 3.0f, map.getBottom() - 3.0f,
                    map.getRight() - 3.0f, map.getY() + 3.0f, 1.3f);
    }
    else if (articulation == taikor::Articulation::DonRim)
    {
        juce::Path rim;
        rim.addCentredArc (map.getCentreX(), map.getCentreY(), mapSize * 0.61f,
                           mapSize * 0.61f, 0.0f, -0.5f, 1.55f, true);
        g.strokePath (rim, juce::PathStrokeType (1.5f));
    }
    g.setColour (ink);
    g.setFont (controlFont (32.0f * editorScale (*this), true));
    g.drawText (keyText, bounds.withTrimmedLeft (bounds.getWidth() * 0.47f),
                juce::Justification::centredLeft, false);
    if (hasKeyboardFocus (true))
    {
        g.setColour (selected || active ? washiColour : accentColour);
        g.drawRoundedRectangle (bounds.reduced (3.0f), 2.0f, 1.8f);
    }
}

std::unique_ptr<juce::AccessibilityHandler> TaikorPad::createAccessibilityHandler()
{
    class PadHandler final : public juce::AccessibilityHandler
    {
    public:
        explicit PadHandler (TaikorPad& padToUse)
            : juce::AccessibilityHandler (
                  padToUse, juce::AccessibilityRole::button,
                  juce::AccessibilityActions().addAction (
                      juce::AccessibilityActionType::press,
                      [&padToUse] { padToUse.triggerClick(); })),
              pad (padToUse)
        {
        }

        // Asked of the pad each time rather than captured when the handler is
        // built. The octave strip moves every pad to a different MIDI note and
        // the handler outlives that change, so a captured string would announce
        // whichever octave happened to be selected when the editor opened.
        juce::String getHelp() const override { return pad.accessibleHelpText(); }

    private:
        TaikorPad& pad;
    };

    return std::make_unique<PadHandler> (*this);
}

// ---------------------------------------------------------------------------
// Drum row selector
// ---------------------------------------------------------------------------

TaikorDrumButton::TaikorDrumButton (int octaveOffsetToUse, juce::Image atlas)
    : juce::Button (octaveName (octaveOffsetToUse) + " "
                    + octaveDescription (octaveOffsetToUse)),
      noteName (octaveName (octaveOffsetToUse)),
      drumName (octaveDescription (octaveOffsetToUse))
{
    // Individual source windows keep the generated studies and their stands
    // intact despite their different sizes and spacing in the atlas.
    const std::array<juce::Rectangle<float>, 4> crops {{
        { 20.0f / 2172.0f, 32.0f / 724.0f, 600.0f / 2172.0f, 668.0f / 724.0f },
        { 664.0f / 2172.0f, 140.0f / 724.0f, 480.0f / 2172.0f, 550.0f / 724.0f },
        { 1228.0f / 2172.0f, 112.0f / 724.0f, 425.0f / 2172.0f, 574.0f / 724.0f },
        { 1704.0f / 2172.0f, 340.0f / 724.0f, 451.0f / 2172.0f, 340.0f / 724.0f },
    }};

    const auto cropIndex = static_cast<std::size_t> (
        juce::jlimit (0, 3, octaveOffsetToUse - taikor::lowestOctaveOffset));
    const auto crop = crops[cropIndex];
    if (atlas.isValid())
    {
        const auto source = juce::Rectangle<int> (
            juce::roundToInt (crop.getX() * static_cast<float> (atlas.getWidth())),
            juce::roundToInt (crop.getY() * static_cast<float> (atlas.getHeight())),
            juce::roundToInt (crop.getWidth() * static_cast<float> (atlas.getWidth())),
            juce::roundToInt (crop.getHeight() * static_cast<float> (atlas.getHeight())))
                                .getIntersection (atlas.getBounds());
        drumPainting = atlas.getClippedImage (source);
    }

    setTitle (noteName + " " + drumName);
    setDescription (octaveSummary (octaveOffsetToUse));
    setTooltip (drumName + " - " + octaveSummary (octaveOffsetToUse)
                + ". Selects this drum for the live head readout.");
}

void TaikorDrumButton::paintButton (juce::Graphics& g, bool isMouseOver,
                                    bool isButtonDown)
{
    const auto bounds = getLocalBounds().toFloat().reduced (1.0f);
    g.setColour (getToggleState() ? juce::Colour (0xffe5dcc8) : panelColour);
    g.fillRoundedRectangle (bounds, 3.0f);
    if (isMouseOver || isButtonDown)
    {
        g.setColour (accentColour.withAlpha (isButtonDown ? 0.14f : 0.06f));
        g.fillRoundedRectangle (bounds, 3.0f);
    }
    if (getToggleState())
    {
        g.setColour (accentColour);
        g.fillRect (bounds.withWidth (3.0f));
    }
    const auto art = bounds.withWidth (bounds.getWidth() * 0.23f).reduced (5.0f, 2.0f);
    if (drumPainting.isValid())
        g.drawImage (drumPainting, art, juce::RectanglePlacement::centred);

    const auto scale = editorScale (*this);
    auto copy = bounds.withTrimmedLeft (bounds.getWidth() * 0.24f).reduced (3.0f, 0.0f);
    copy = copy.withHeight (44.0f * scale).withCentre (copy.getCentre());
    g.setColour (accentColour);
    g.setFont (controlFont (18.0f * scale, true));
    g.drawText (noteName, copy.removeFromTop (18.0f * scale),
                juce::Justification::centredLeft, false);
    g.setColour (textColour);
    g.setFont (fitFont (displayFont (26.0f * scale, juce::Font::bold), drumName, copy.getWidth()));
    g.drawText (drumName, copy, juce::Justification::centredLeft, false);
    if (hasKeyboardFocus (true))
    {
        g.setColour (accentColour);
        g.drawRoundedRectangle (bounds.reduced (1.5f), 2.0f, 1.5f);
    }
}

// ---------------------------------------------------------------------------
// Knob
// ---------------------------------------------------------------------------

TaikorKnob::TaikorKnob (juce::String name, ValueStyle style, VisualRole roleToUse)
    : role (roleToUse)
{
    setName (name);
    slider.setSliderStyle (juce::Slider::RotaryHorizontalVerticalDrag);
    slider.setTextBoxStyle (juce::Slider::TextBoxBelow, false, 84, 20);
    slider.setColour (juce::Slider::rotarySliderFillColourId, roleColour (roleToUse));
    slider.setColour (juce::Slider::textBoxTextColourId, textColour);
    slider.setColour (juce::Slider::textBoxOutlineColourId,
                      juce::Colours::transparentBlack);
    slider.setColour (juce::Slider::textBoxBackgroundColourId,
                      juce::Colours::transparentBlack);

    switch (style)
    {
        case ValueStyle::Plain:       break;
        case ValueStyle::Percent:     slider.setTextValueSuffix (" %"); break;
        case ValueStyle::Semitones:   slider.setTextValueSuffix (" st"); break;
        case ValueStyle::Degrees:     slider.setTextValueSuffix (" deg"); break;
        case ValueStyle::Centimetres: slider.setTextValueSuffix (" cm"); break;
        case ValueStyle::Decibels:    slider.setTextValueSuffix (" dB"); break;
    }

    addAndMakeVisible (slider);

    label.setText (name, juce::dontSendNotification);
    label.setJustificationType (juce::Justification::centred);
    label.setBorderSize (juce::BorderSize<int> (1, 1, 1, 1));
    label.setColour (juce::Label::textColourId, mutedText);
    label.setFont (juce::Font (juce::FontOptions (10.5f).withStyle ("Bold")));
    label.setInterceptsMouseClicks (false, false);
    addAndMakeVisible (label);
}

void TaikorKnob::setLabelText (const juce::String& text, const juce::String& description)
{
    label.setText (text, juce::dontSendNotification);
    setTooltip (description);
    slider.setTooltip (description);
    slider.setTitle (text);
    slider.setDescription (description);
}

void TaikorKnob::resized()
{
    auto bounds = getLocalBounds().reduced (1, 0);
    label.setFont (controlFont (10.5f, true));
    label.setBounds (bounds.removeFromTop (16));
    const auto valueHeight = juce::roundToInt (30.0f * editorScale (*this));
    const bool wide = static_cast<float> (bounds.getWidth())
                    > static_cast<float> (bounds.getHeight()) * 1.6f;
    slider.setTextBoxStyle (wide ? juce::Slider::TextBoxRight : juce::Slider::TextBoxBelow,
                            false, wide ? juce::roundToInt (static_cast<float> (bounds.getWidth()) * 0.58f)
                                        : bounds.getWidth(), valueHeight);
    slider.setBounds (bounds);
}

// ---------------------------------------------------------------------------
// Indexed switches
// ---------------------------------------------------------------------------

namespace
{
std::vector<TaikorChoiceSwitch::Choice> indexedChoices (
    const juce::RangedAudioParameter& parameter)
{
    std::vector<TaikorChoiceSwitch::Choice> result;
    const auto names = parameter.getAllValueStrings();
    for (int index = 0; index < names.size(); ++index)
        result.push_back ({ names[index], static_cast<float> (index) });
    return result;
}
} // namespace

TaikorChoiceSwitch::TaikorChoiceSwitch (juce::String name,
                                        juce::RangedAudioParameter& parameter,
                                        const juce::String& description)
    : TaikorChoiceSwitch (name, parameter, description, indexedChoices (parameter))
{
}

TaikorChoiceSwitch::TaikorChoiceSwitch (juce::String name,
                                        juce::RangedAudioParameter& parameter,
                                        const juce::String& description,
                                        std::vector<Choice> choicesToUse)
    : choices (std::move (choicesToUse)),
      attachment (parameter, [this] (float value) { selectNearest (value); })
{
    setName (name);
    setTitle (name);
    setDescription (description);
    label.setText (name, juce::dontSendNotification);
    label.setFont (controlFont (11.0f, true));
    label.setColour (juce::Label::textColourId, mutedText);
    label.setJustificationType (juce::Justification::centredLeft);
    label.setInterceptsMouseClicks (false, false);
    addAndMakeVisible (label);

    for (std::size_t index = 0; index < choices.size(); ++index)
    {
        const auto& choice = choices[index];
        auto* button = buttons.add (new juce::TextButton (choice.label));
        button->setName (choice.label);
        button->setTitle (name + ": " + choice.label);
        button->setDescription (description);
        button->setTooltip (choice.label + ". " + description);
        button->setClickingTogglesState (true);
        button->setRadioGroupId (2, juce::dontSendNotification);
        button->setWantsKeyboardFocus (true);
        button->setColour (juce::TextButton::buttonColourId, panelRaised);
        button->setColour (juce::TextButton::buttonOnColourId, accentColour);
        button->setColour (juce::TextButton::textColourOffId, textColour);
        button->setColour (juce::TextButton::textColourOnId, washiColour);
        const float value = choice.value;
        button->onClick = [this, value]
        {
            attachment.setValueAsCompleteGesture (value);
        };
        addAndMakeVisible (button);
    }
    attachment.sendInitialUpdate();
}

void TaikorChoiceSwitch::selectNearest (float value)
{
    int selected = -1;
    float nearest = std::numeric_limits<float>::infinity();
    for (std::size_t index = 0; index < choices.size(); ++index)
    {
        const float distance = std::abs (choices[index].value - value);
        if (distance < nearest)
        {
            nearest = distance;
            selected = static_cast<int> (index);
        }
    }
    for (int index = 0; index < buttons.size(); ++index)
        buttons[index]->setToggleState (index == selected, juce::dontSendNotification);
}

void TaikorChoiceSwitch::resized()
{
    const auto scale = editorScale (*this);
    auto bounds = getLocalBounds();
    label.setBounds (bounds.removeFromTop (16));
    bounds.removeFromTop (juce::roundToInt (4.0f * scale));
    auto buttonHeight = juce::jmin (juce::roundToInt (36.0f * scale), bounds.getHeight());
    auto centreY = bounds.getCentreY();
    if (getHeight() > juce::roundToInt (70.0f * scale))
    {
        // Beside knobs, center the choices on the dial, leaving the same
        // caption and value clearance as the neighboring continuous controls.
        const auto dial = getLocalBounds().withTrimmedTop (16)
                              .withTrimmedBottom (juce::roundToInt (30.0f * scale));
        buttonHeight = juce::jmin (buttonHeight, dial.getHeight());
        centreY = dial.getCentreY();
    }
    bounds = bounds.withHeight (buttonHeight).withCentre ({ bounds.getCentreX(), centreY });
    const auto layout = taikor::ui::rowLayout (bounds.getWidth(), buttons.size(), 3, buttons.size());
    for (int index = 0; index < buttons.size(); ++index)
        buttons[index]->setBounds (bounds.getX() + taikor::ui::cellOffset (layout, 3, index),
                                    bounds.getY(), layout.cellSize, bounds.getHeight());
}

// ---------------------------------------------------------------------------
// Head display
// ---------------------------------------------------------------------------

TaikorHeadDisplay::TaikorHeadDisplay()
{
    setInterceptsMouseClicks (false, false);
    setTitle ("Drum head");
}

void TaikorHeadDisplay::setStrike (float normalisedRadius, float angleRadians,
                                   float level, taikor::Articulation articulation)
{
    const auto radius = taikor::ui::clamp (normalisedRadius, 0.0f, 1.0f);
    const auto bounded = taikor::ui::clamp (level, 0.0f, 1.0f);

    // Every drawn property has to be in this decision, not just the level. A
    // roll played at one velocity moves the marker around the head and changes
    // the articulation without changing the level at all, and testing the level
    // alone left the display showing the first stroke of the roll for as long
    // as it went on.
    const bool moved = std::abs (radius - strikeRadius) >= 0.002f
                    || std::abs (angleRadians - strikeAngle) >= 0.002f
                    || articulation != lastArticulation
                    || std::abs (bounded - strikeLevel) >= 0.002f;

    strikeRadius = radius;
    strikeAngle = angleRadians;
    lastArticulation = articulation;
    strikeLevel = bounded;

    if (moved)
        repaint();
}

void TaikorHeadDisplay::setMicrophones (float spread, float normalisedDistance)
{
    if (std::abs (spread - micSpread) < 0.002f
        && std::abs (normalisedDistance - micDistance) < 0.002f)
        return;
    micSpread = taikor::ui::clamp (spread, 0.0f, 1.0f);
    micDistance = taikor::ui::clamp (normalisedDistance, 0.0f, 1.0f);
    repaint();
}

void TaikorHeadDisplay::setMeasurements (float fundamentalHz, float breathingHz,
                                         float diameterCentimetres, float tailSeconds)
{
    if (std::abs (fundamentalHz - fundamental) < 0.05f
        && std::abs (breathingHz - breathing) < 0.05f
        && std::abs (diameterCentimetres - diameter) < 0.05f
        && std::abs (tailSeconds - tail) < 0.005f)
        return;
    fundamental = fundamentalHz;
    breathing = breathingHz;
    diameter = diameterCentimetres;
    tail = tailSeconds;
    repaint();
}

void TaikorHeadDisplay::paint (juce::Graphics& g)
{
    const auto bounds = getLocalBounds().toFloat();
    const auto centre = bounds.getCentre();
    const auto headRadius = juce::jmin (bounds.getWidth(), bounds.getHeight()) * 0.37f;

    // An incomplete ink circle gives the technical head map a place in the
    // same visual world as the row paintings, without obscuring live data.
    juce::Path enso;
    enso.addCentredArc (centre.x, centre.y, headRadius * 1.23f,
                        headRadius * 1.23f, -0.08f, -0.15f, pi * 1.72f, true);
    g.setColour (textColour.withAlpha (0.12f));
    g.strokePath (enso, juce::PathStrokeType (headRadius * 0.10f,
                                              juce::PathStrokeType::curved,
                                              juce::PathStrokeType::rounded));

    // The shell, seen edge on behind the head.
    g.setColour (textColour.withAlpha (0.94f));
    g.fillEllipse (centre.x - headRadius * 1.10f, centre.y - headRadius * 1.10f,
                   headRadius * 2.20f, headRadius * 2.20f);
    g.setColour (brassColour.withAlpha (0.48f));
    g.drawEllipse (centre.x - headRadius * 1.10f, centre.y - headRadius * 1.10f,
                   headRadius * 2.20f, headRadius * 2.20f, 1.4f);

    // The hide.
    juce::ColourGradient hide { hideColour.brighter (0.10f), centre.x,
                                centre.y - headRadius, hideColour.darker (0.30f),
                                centre.x, centre.y + headRadius, false };
    g.setGradientFill (hide);
    g.fillEllipse (centre.x - headRadius, centre.y - headRadius,
                   headRadius * 2.0f, headRadius * 2.0f);

    // The tack ring that holds the head on.
    g.setColour (juce::Colour { 0xff554638 });
    for (int tack = 0; tack < 32; ++tack)
    {
        const auto angle = static_cast<float> (tack) * 2.0f * pi / 32.0f;
        const auto point = taikor::ui::headPointFor (0.945f, angle);
        g.fillEllipse (centre.x + point.x * headRadius - 1.6f,
                       centre.y + point.y * headRadius - 1.6f, 3.2f, 3.2f);
    }

    // The ripple where the last stroke landed. Its size follows the stroke's
    // own radius, so a Ka really does land out by the tacks.
    if (strikeLevel > 0.004f)
    {
        const auto point = taikor::ui::headPointFor (strikeRadius, strikeAngle);
        const auto x = centre.x + point.x * headRadius;
        const auto y = centre.y + point.y * headRadius;

        for (int ring = 0; ring < 3; ++ring)
        {
            const auto spread = headRadius * (0.10f + 0.11f * static_cast<float> (ring))
                              * (0.6f + 0.9f * strikeLevel);
            const auto alpha = strikeLevel * (0.55f - 0.15f * static_cast<float> (ring));
            g.setColour (accentColour.withAlpha (juce::jmax (0.0f, alpha)));
            g.drawEllipse (x - spread, y - spread, spread * 2.0f, spread * 2.0f, 1.6f);
        }

        g.setColour (accentColour.withAlpha (juce::jmin (1.0f, strikeLevel)));
        const auto dot = headRadius * 0.045f;
        g.fillEllipse (x - dot, y - dot, dot * 2.0f, dot * 2.0f);
    }

    g.setColour (hideColour.darker (0.55f));
    g.drawEllipse (centre.x - headRadius, centre.y - headRadius,
                   headRadius * 2.0f, headRadius * 2.0f, 1.2f);

    // The close pair, where the microphone controls put it. The dots sit at the
    // radius and arc the model actually reads the head at - the same constants
    // resolveDrumFor uses - and they fade as the pair backs off. The arc kept
    // the engine's old 2.2-radian spread after the model narrowed to a true
    // close pair, so at full spread the display showed capsules a hundred and
    // twenty-six degrees apart while the drum being heard used fifty.
    const auto micRingRadius = headRadius * (0.10f + 0.68f * micSpread);
    constexpr float micReference = 0.60f;
    const auto separation = 0.9f * micSpread;
    const auto proximity = 1.0f - micDistance;

    for (int side = 0; side < 2; ++side)
    {
        const auto angle = micReference
                         + (side == 0 ? 0.5f : -0.5f) * separation;
        const auto point = taikor::ui::headPointFor (1.0f, angle);
        const auto x = centre.x + point.x * micRingRadius;
        const auto y = centre.y + point.y * micRingRadius;
        const auto size = headRadius * (0.045f + 0.030f * proximity);

        g.setColour (roleColour (TaikorKnob::VisualRole::Microphone)
                         .withAlpha (0.35f + 0.5f * proximity));
        g.fillEllipse (x - size, y - size, size * 2.0f, size * 2.0f);
        g.setColour (roleColour (TaikorKnob::VisualRole::Microphone));
        g.drawEllipse (x - size, y - size, size * 2.0f, size * 2.0f, 1.2f);

        g.setFont (juce::Font (juce::FontOptions (9.0f).withStyle ("Bold")));
        g.drawText (side == 0 ? "L" : "R",
                    juce::Rectangle<float> (x - size, y - size, size * 2.0f, size * 2.0f),
                    juce::Justification::centred, false);
    }


}

std::unique_ptr<juce::AccessibilityHandler>
TaikorHeadDisplay::createAccessibilityHandler()
{
    class HeadHandler final : public juce::AccessibilityHandler
    {
    public:
        explicit HeadHandler (TaikorHeadDisplay& display)
            : juce::AccessibilityHandler (display, juce::AccessibilityRole::staticText),
              owner (display)
        {
        }

        juce::String getHelp() const override
        {
            const auto stroke =
                taikor::getArticulationDisplayName (owner.lastArticulation);
            return "Drum head. " + juce::String (owner.diameter, 1)
                 + " centimetre head sounding "
                 + (owner.fundamental > 0.0f
                        ? juce::String (owner.fundamental, 1) + " hertz"
                        : juce::String ("no membrane tone at this sample rate"))
                 + " with a " + juce::String (owner.breathing, 1)
                 + " hertz breathing mode, and a " + juce::String (owner.tail, 2)
                 + " second tail. Last stroke: "
                 + juce::String (stroke.data(), stroke.size()) + " at "
                 + juce::String (juce::roundToInt (owner.strikeRadius * 100.0f))
                 + " per cent of the radius.";
        }

    private:
        TaikorHeadDisplay& owner;
    };

    return std::make_unique<HeadHandler> (*this);
}

// ---------------------------------------------------------------------------
// Status and metering
// ---------------------------------------------------------------------------

TaikorStatusDisplay::TaikorStatusDisplay()
{
    setInterceptsMouseClicks (false, false);
    setTitle ("Engine status");
}

void TaikorStatusDisplay::setStatus (int activeVoices, bool ready, double sampleRate)
{
    if (voices == activeVoices && isReady == ready
        && std::abs (rate - sampleRate) < 0.5)
        return;
    voices = activeVoices;
    isReady = ready;
    rate = sampleRate;
    repaint();
}

void TaikorStatusDisplay::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();
    bounds = bounds.reduced (10.0f, 4.0f);
    const auto statusDot = juce::jmin (6.0f, bounds.getHeight() * 0.32f);
    g.setColour (isReady ? brassColour : mutedText.withAlpha (0.42f));
    g.fillEllipse (bounds.getX(), bounds.getCentreY() - statusDot * 0.5f,
                   statusDot, statusDot);
    bounds.removeFromLeft (statusDot + 7.0f);
    g.setColour (isReady ? textColour : mutedText);

    const juce::String rateText = rate > 0.0
        ? juce::String (rate / 1000.0, 1) + " kHz"
        : juce::String ("- kHz");
    const auto status = juce::String (juce::jmax (0, voices)) + " voices  "
                      + juce::String::fromUTF8 ("\xc2\xb7") + "  " + rateText;
    g.setFont (fitFont (controlFont (20.0f * editorScale (*this)), status, bounds.getWidth()));
    g.drawText (status, bounds, juce::Justification::centredLeft, false);
}

std::unique_ptr<juce::AccessibilityHandler>
TaikorStatusDisplay::createAccessibilityHandler()
{
    class StatusHandler final : public juce::AccessibilityHandler
    {
    public:
        explicit StatusHandler (TaikorStatusDisplay& display)
            : juce::AccessibilityHandler (display, juce::AccessibilityRole::staticText),
              owner (display)
        {
        }

        juce::String getHelp() const override
        {
            return juce::String (juce::jmax (0, owner.voices)) + " voices sounding at "
                 + (owner.rate > 0.0 ? juce::String (owner.rate / 1000.0, 1)
                                     : juce::String ("unknown"))
                 + " kilohertz.";
        }

    private:
        TaikorStatusDisplay& owner;
    };

    return std::make_unique<StatusHandler> (*this);
}

TaikorMeter::TaikorMeter()
{
    setInterceptsMouseClicks (false, false);
    setTitle ("Output meter");
    leftBallistics.reset();
    rightBallistics.reset();
}

void TaikorMeter::setLevels (float leftLinear, float rightLinear)
{
    constexpr float updateRate = 30.0f;
    const auto attack = taikor::ui::onePoleCoefficient (0.012f, updateRate);
    const auto release = taikor::ui::onePoleCoefficient (0.240f, updateRate);
    const auto peakFall = taikor::ui::decayMultiplier (-12.0f, 1.0f, updateRate);
    constexpr float hold = 18.0f;

    leftBallistics.update (taikor::ui::meterPositionForLinear (leftLinear, floorDecibels),
                           attack, release, peakFall, hold);
    rightBallistics.update (
        taikor::ui::meterPositionForLinear (rightLinear, floorDecibels),
        attack, release, peakFall, hold);

    const bool clipping = leftLinear >= 0.999f || rightLinear >= 0.999f;
    const int clipState = clipping ? 1 : 0;
    if (clipState != announcedClipState)
    {
        announcedClipState = clipState;
        if (clipping)
            if (auto* handler = getAccessibilityHandler())
                handler->notifyAccessibilityEvent (juce::AccessibilityEvent::valueChanged);
    }

    if (std::abs (leftBallistics.level - lastPaintedLeft) < 0.003f
        && std::abs (rightBallistics.level - lastPaintedRight) < 0.003f
        && std::abs (leftBallistics.peak - lastPaintedLeftPeak) < 0.003f
        && std::abs (rightBallistics.peak - lastPaintedRightPeak) < 0.003f)
        return;

    lastPaintedLeft = leftBallistics.level;
    lastPaintedRight = rightBallistics.level;
    lastPaintedLeftPeak = leftBallistics.peak;
    lastPaintedRightPeak = rightBallistics.peak;
    repaint();
}

void TaikorMeter::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();
    bounds = bounds.reduced (1.0f, 3.0f);
    auto channelLabels = bounds.removeFromLeft (22.0f * editorScale (*this));
    const auto barHeight = (bounds.getHeight() - 3.0f) * 0.5f;

    const auto drawBar = [&] (juce::Rectangle<float> area,
                              const taikor::ui::MeterBallistics& ballistics)
    {
        g.setColour (textColour.withAlpha (0.15f));
        g.fillRect (area);

        const auto filled = area.withWidth (area.getWidth() * ballistics.level);
        const auto hot = ballistics.level
                       > taikor::ui::meterPositionForLinear (0.708f, floorDecibels);
        g.setColour (hot ? accentColour : roleColour (TaikorKnob::VisualRole::Microphone));
        g.fillRect (filled);

        const auto peakX = area.getX() + area.getWidth() * ballistics.peak;
        g.setColour (textColour.withAlpha (0.85f));
        g.fillRect (juce::Rectangle<float> (peakX - 1.0f, area.getY(), 2.0f,
                                            area.getHeight()));
    };

    g.setColour (mutedText);
    g.setFont (controlFont (16.0f * editorScale (*this), true));
    g.drawText ("L", channelLabels.removeFromTop (barHeight), juce::Justification::centredLeft, false);
    channelLabels.removeFromTop (3.0f);
    g.drawText ("R", channelLabels, juce::Justification::centredLeft, false);
    drawBar (bounds.removeFromTop (barHeight), leftBallistics);
    bounds.removeFromTop (3.0f);
    drawBar (bounds.removeFromTop (barHeight), rightBallistics);
}

std::unique_ptr<juce::AccessibilityHandler> TaikorMeter::createAccessibilityHandler()
{
    class MeterHandler final : public juce::AccessibilityHandler
    {
    public:
        explicit MeterHandler (TaikorMeter& meterToUse)
            : juce::AccessibilityHandler (meterToUse,
                                          juce::AccessibilityRole::staticText),
              owner (meterToUse)
        {
        }

        juce::String getHelp() const override
        {
            const auto decibels = [] (float position)
            {
                return juce::String (
                    juce::Decibels::gainToDecibels (
                        taikor::ui::linearForMeterPosition (position, floorDecibels)),
                    1);
            };
            return "Output " + decibels (owner.leftBallistics.level) + " dB left, "
                 + decibels (owner.rightBallistics.level) + " dB right.";
        }

    private:
        TaikorMeter& owner;
    };

    return std::make_unique<MeterHandler> (*this);
}

// ---------------------------------------------------------------------------
// Editor
// ---------------------------------------------------------------------------

TaikorAudioProcessorEditor::TaikorAudioProcessorEditor (TaikorAudioProcessor& processorToUse)
    : juce::AudioProcessorEditor (&processorToUse),
      audioProcessor (processorToUse),
      tooltipWindow (this, 650)
{
    setLookAndFeel (&lookAndFeel);

    drumAtlas = juce::ImageFileFormat::loadFrom (
        BinaryData::taikordrumstudies_png,
        static_cast<std::size_t> (BinaryData::taikordrumstudies_pngSize));
    backgroundPainting = juce::ImageFileFormat::loadFrom (
        BinaryData::taikorcoastalprint_png,
        static_cast<std::size_t> (BinaryData::taikorcoastalprint_pngSize));

    logoLabel.setText ("TAIKOR", juce::dontSendNotification);
    logoLabel.setFont (displayFont (46.0f, juce::Font::bold));
    logoLabel.setColour (juce::Label::textColourId, textColour);
    logoLabel.setJustificationType (juce::Justification::centredLeft);
    addAndMakeVisible (logoLabel);

    editionLabel.setText ("PHYSICALLY\nMODELED TAIKO", juce::dontSendNotification);
    editionLabel.setFont (juce::Font (juce::FontOptions (10.5f)));
    editionLabel.setColour (juce::Label::textColourId, mutedText);
    editionLabel.setJustificationType (juce::Justification::centredLeft);
    addAndMakeVisible (editionLabel);

    addAndMakeVisible (statusDisplay);
    addAndMakeVisible (meter);
    limiterLabel.setText ("OUTPUT  /  LIMIT -1 dB", juce::dontSendNotification);
    limiterLabel.setFont (juce::Font (juce::FontOptions (9.5f).withStyle ("Bold")));
    limiterLabel.setColour (juce::Label::textColourId, mutedText);
    limiterLabel.setJustificationType (juce::Justification::centredLeft);
    addAndMakeVisible (limiterLabel);

    panicButton.setColour (juce::TextButton::buttonColourId, accentDim);
    panicButton.setColour (juce::TextButton::textColourOffId, washiColour);
    panicButton.setTooltip ("Silence every sounding stroke immediately");
    panicButton.onClick = [this] { audioProcessor.requestPanic(); };
    addAndMakeVisible (panicButton);

    sealLabel.setText (juce::String::fromUTF8 ("\xe9\xbc\x93"),
                       juce::dontSendNotification);
    sealLabel.setFont (displayFont (19.0f, juce::Font::bold));
    sealLabel.setColour (juce::Label::backgroundColourId, accentColour);
    sealLabel.setColour (juce::Label::textColourId, washiColour);
    sealLabel.setJustificationType (juce::Justification::centred);
    sealLabel.setInterceptsMouseClicks (false, false);
    addAndMakeVisible (sealLabel);

    addAndMakeVisible (headDisplay);

    headCaption.setText ("LIVE HEAD RESPONSE", juce::dontSendNotification);
    headCaption.setFont (controlFont (11.0f, true));
    headCaption.setColour (juce::Label::textColourId, mutedText);
    headCaption.setJustificationType (juce::Justification::centred);
    addAndMakeVisible (headCaption);

    gridCaption.setText ("THE DRUM ENSEMBLE", juce::dontSendNotification);
    gridCaption.setFont (displayFont (19.0f, juce::Font::bold));
    gridCaption.setColour (juce::Label::textColourId, textColour);
    gridCaption.setJustificationType (juce::Justification::centredLeft);
    addAndMakeVisible (gridCaption);

    for (std::size_t artIdx = 0; artIdx < taikor::articulationCount; ++artIdx)
    {
        const auto articulation = static_cast<taikor::Articulation> (artIdx);
        const auto& metadata = taikor::getArticulationMetadata (articulation);
        auto& label = strokeLabels[artIdx];
        label.setText (juce::String (metadata.displayName.data(),
                                     metadata.displayName.size()).toUpperCase(),
                       juce::dontSendNotification);
        label.setFont (controlFont (11.0f, true));
        label.setColour (juce::Label::textColourId, mutedText);
        label.setJustificationType (juce::Justification::centred);
        label.setTooltip (juce::String (metadata.description.data(),
                                        metadata.description.size()));
        addAndMakeVisible (label);
    }

    std::size_t padIdx = 0;
    for (int octave = taikor::lowestOctaveOffset; octave <= taikor::highestOctaveOffset; ++octave)
    {
        for (std::size_t artIdx = 0; artIdx < taikor::articulationCount; ++artIdx)
        {
            const auto articulation = static_cast<taikor::Articulation> (artIdx);
            auto pad = std::make_unique<TaikorPad> (articulation, octave);
            pad->onClick = [this, articulation, octave]
            {
                selectOctave (octave);
                audioProcessor.triggerFromUi (articulation, octave);
            };
            addAndMakeVisible (*pad);
            pads[padIdx++] = std::move (pad);
        }
    }

    for (int index = 0; index < octaveCount; ++index)
    {
        const auto octave = taikor::lowestOctaveOffset + index;
        auto button = std::make_unique<TaikorDrumButton> (octave, drumAtlas);
        button->setClickingTogglesState (true);
        button->setRadioGroupId (1, juce::dontSendNotification);
        button->onClick = [this, octave] { selectOctave (octave); };
        addAndMakeVisible (*button);
        octaveButtons[static_cast<std::size_t> (index)] = std::move (button);
    }

    const auto deckLabel = [this] (juce::Label& label, const juce::String& text)
    {
        label.setText (text, juce::dontSendNotification);
        label.setFont (displayFont (17.0f, juce::Font::bold));
        label.setColour (juce::Label::textColourId, textColour);
        label.setJustificationType (juce::Justification::centredLeft);
        addAndMakeVisible (label);
    };

    deckLabel (drumDeckLabel, "BODY & TUNING");
    deckLabel (strokeDeckLabel, "THE STROKE");
    deckLabel (microphoneDeckLabel, "MICROPHONES & OUTPUT");

    namespace ids = taikor::parameters;
    addKnob (sizeKnob, ids::headDiameter,
             "Head diameter. Pitch follows one over the radius, so this moves the "
             "whole drum without changing the ratios between its modes.");
    addKnob (depthKnob, ids::bodyDepth,
             "Body depth. A shallow body has a stiffer air spring, which pushes the "
             "breathing mode further above the fundamental.");
    addKnob (tensionKnob, ids::tension,
             "Head tension. Wave speed is the square root of tension over the head's "
             "areal density, so this and the head material together set the pitch.");
    addKnob (headMaterialKnob, ids::headMaterial,
             "Head material, from a thin synthetic film to a thick cowhide. Sets both "
             "the head's weight and how much it loses per cycle.");
    addKnob (shellMaterialKnob, ids::shellMaterial,
             "Shell material, from light laminated ply to dense carved zelkova. Moves "
             "the body's ring modes, their Q, and how much the rim absorbs.");
    addKnob (headDampingKnob, ids::headDamping,
             "Extra loss in the head on top of the material's own.");
    addKnob (pitchKnob, ids::pitch,
             "Musical transposition, applied as head tension because that is what "
             "tuning a drum is.");

    addKnob (hardnessKnob, ids::bachiHardness,
             "Bachi hardness, from a felt-wrapped beater to seasoned oak. Sets the "
             "Hertz contact stiffness and therefore how long the stick stays down.");
    addKnob (strikePositionKnob, ids::strikePosition,
             "Moves every stroke towards the centre or towards the rim, on top of "
             "the position its own articulation already asks for.");
    addKnob (strikeAzimuthKnob, ids::strikeAzimuth,
             "Turns the strike around the head. CC16 overrides this angle for "
             "sample-accurate left and right hand placement.");
    performerSwitch = std::make_unique<TaikorChoiceSwitch> (
        "PERFORMER", *audioProcessor.parameters.getParameter (ids::performer),
        "Four repeatable players, each with a stable touch and strike character.");
    addAndMakeVisible (*performerSwitch);
    addKnob (velocityDepthKnob, ids::velocityDepth,
             "How far MIDI velocity moves the impact speed. The timbre follows on "
             "its own: contact time goes as impact speed to the minus one fifth.");
    velocityCurveSwitch = std::make_unique<TaikorChoiceSwitch> (
        "CURVE", *audioProcessor.parameters.getParameter (ids::velocityCurve),
        "Shapes MIDI velocity before impact speed: Soft opens up quiet playing, "
        "Linear leaves it unchanged, and Hard asks for a firmer hit.",
        std::vector<TaikorChoiceSwitch::Choice> {
            { "Soft", -1.0f }, { "Linear", 0.0f }, { "Hard", 1.0f } });
    addAndMakeVisible (*velocityCurveSwitch);
    addKnob (tensionModKnob, ids::tensionModulation,
             "Attack pitch glide. A hard stroke stretches the head, raising its "
             "tension until the stroke decays.");
    addKnob (strikeNoiseKnob, ids::strikeNoise,
             "Level of the broadband contact noise the stick makes on the hide.");
    addKnob (humaniseKnob, ids::humanise,
             "Per-stroke variation in position, angle, impact speed and contact time. "
             "At 0, subtle speed and contact differences still keep repeated hits alive.");
    {
        std::vector<TaikorChoiceSwitch::Choice> players;
        for (int count = 1; count <= taikor::maximumEnsembleSize; ++count)
            players.push_back ({ juce::String (count), static_cast<float> (count) });
        ensembleSizeSwitch = std::make_unique<TaikorChoiceSwitch> (
            "PLAYERS", *audioProcessor.parameters.getParameter (ids::ensembleSize),
            "Number of players on separate copies of each drum, from 1 to 8. "
            "One keeps the solo sound. Two sit left/right; three add centre; "
            "four sit at 100% left, 50% left, 50% right and 100% right. "
            "Moves use 15 ms smoothing; removed tails keep their positions. "
            "Width at 50% preserves the stage; 0% sums it to mono.",
            std::move (players));
        addAndMakeVisible (*ensembleSizeSwitch);
    }
    addKnob (ensembleVariationKnob, ids::ensembleVariation,
             "Differences between the players: timing, hit placement and, for each "
             "companion, its own hide, tension, shell and damping. At 0 every "
             "player is the same drum struck at the same instant; turn up for a "
             "looser, more varied ensemble. Inactive with one player.");
    ensembleVariationKnob.slider.setTitle ("Ensemble Variation");
    drumLayoutSwitch = std::make_unique<TaikorChoiceSwitch> (
        "DRUM LAYOUT", *audioProcessor.parameters.getParameter (ids::octaveBody),
        "1 Drum retunes one design across four rows. 4 Drums uses four taiko families. "
        "Each row rings independently; there is no shared room or inter-drum resonance.");
    addAndMakeVisible (*drumLayoutSwitch);

    addKnob (micDistanceKnob, ids::micDistance,
             "How far the close pair stands off the head. Near in it reads the shape "
             "of the membrane and the image opens; further back only what the drum "
             "radiates survives, and the pair narrows.");
    addKnob (micSpreadKnob, ids::micSpread,
             "How far apart the pair is spread across the head. This is where the "
             "stereo comes from: two points see different signs of every mode with a "
             "circumferential order.");
    addKnob (widthKnob, ids::stereoWidth, "Width trim on the finished pair.");
    addKnob (driveKnob, ids::drive, "Gentle output-stage saturation.");
    addKnob (outputHighPassKnob, ids::outputHighPass,
             "Gentle 6 dB/octave output high-pass filter, before the limiter. "
             "Off bypasses the filter; turn up to remove low frequencies, up to 500 Hz.");
    addKnob (outputKnob, ids::output, "Output level before the always-on stereo limiter (-1 dBFS peak ceiling).");

    selectOctave (0);

    setResizable (true, true);
    setResizeLimits (minimumWidth, minimumHeight, maximumWidth, maximumHeight);
    getConstrainer()->setFixedAspectRatio (static_cast<double> (designWidth)
                                          / static_cast<double> (designHeight));
    setSize (designWidth, designHeight);

    timerCallback();
    startTimerHz (30);
}

TaikorAudioProcessorEditor::~TaikorAudioProcessorEditor()
{
    stopTimer();
    setLookAndFeel (nullptr);
}

void TaikorAudioProcessorEditor::addKnob (TaikorKnob& knob,
                                          const juce::String& parameterId,
                                          const juce::String& description)
{
    addAndMakeVisible (knob);
    // The knob's own caption is already set from its constructor; this attaches
    // the description used by tooltips and by assistive technology.
    knob.setLabelText (knob.getName().isNotEmpty() ? knob.getName()
                                                   : parameterId.toUpperCase(),
                       description);
    knob.slider.setName (knob.getName());
    attachments.push_back (std::make_unique<SliderAttachment> (
        audioProcessor.parameters, parameterId, knob.slider));

    // Double-click resets a knob to its own parameter's declared default -
    // the value a freshly-inserted instance opens with - rather than to
    // whichever value JUCE's slider would otherwise fall back to. Reading it
    // from the parameter keeps every knob in lockstep with its default should
    // that ever change, instead of a second table of defaults kept here.
    if (const auto* parameter = audioProcessor.parameters.getParameter (parameterId))
    {
        knob.slider.setDoubleClickReturnValue (
            true, parameter->convertFrom0to1 (parameter->getDefaultValue()));
        knob.slider.setTooltip (knob.slider.getTooltip() + " (double-click to reset)");
    }
}

void TaikorAudioProcessorEditor::selectOctave (int octaveOffset)
{
    selectedOctave = juce::jlimit (taikor::lowestOctaveOffset,
                                   taikor::highestOctaveOffset, octaveOffset);

    for (int index = 0; index < octaveCount; ++index)
    {
        const auto octave = taikor::lowestOctaveOffset + index;
        auto& button = octaveButtons[static_cast<std::size_t> (index)];
        if (button != nullptr)
            button->setToggleState (octave == selectedOctave,
                                    juce::dontSendNotification);
    }

    for (auto& pad : pads)
        if (pad != nullptr)
            pad->setSelected (pad->getOctaveOffset() == selectedOctave);
}

TaikorAudioProcessorEditor::LayoutAreas
TaikorAudioProcessorEditor::calculateLayout() const
{
    const auto scale = static_cast<float> (getWidth()) / designWidth;
    const auto rect = [scale] (int x, int y, int w, int h)
    {
        return juce::Rectangle<float> (static_cast<float> (x), static_cast<float> (y),
                                        static_cast<float> (w), static_cast<float> (h))
            .transformedBy (juce::AffineTransform::scale (scale)).toNearestInt();
    };
    LayoutAreas areas;
    areas.header = rect (0, 0, 1280, 96);
    areas.artwork = rect (24, 116, 272, 300);
    areas.head = rect (24, 432, 272, 260);
    areas.switchDeck = rect (24, 708, 272, 148);
    areas.gridArea = rect (312, 116, 944, 300);
    areas.drumDeck = rect (312, 432, 464, 260);
    areas.strokeDeck = rect (792, 432, 464, 260);
    areas.microphoneDeck = rect (312, 708, 944, 148);
    return areas;
}

void TaikorAudioProcessorEditor::paint (juce::Graphics& g)
{
    const auto areas = calculateLayout();
    const auto scale = static_cast<float> (getWidth()) / designWidth;
    g.setGradientFill (juce::ColourGradient (backgroundTop, 0.0f, 0.0f,
                                            backgroundBottom, 0.0f,
                                            static_cast<float> (getHeight()), false));
    g.fillAll();
    g.setColour (panelColour);
    g.fillRect (areas.header);

    for (const auto area : { areas.gridArea, areas.drumDeck, areas.strokeDeck,
                             areas.microphoneDeck, areas.head, areas.switchDeck })
    {
        g.setColour (juce::Colours::black.withAlpha (0.16f));
        g.fillRoundedRectangle (area.toFloat().translated (0.0f, 3.0f * scale), 4.0f * scale);
        g.setColour (panelColour);
        g.fillRoundedRectangle (area.toFloat(), 3.0f * scale);
        g.setColour (panelRaised.withAlpha (0.8f));
        g.drawRoundedRectangle (area.toFloat().reduced (0.5f), 3.0f * scale, 0.7f);
    }
    if (backgroundPainting.isValid())
    {
        juce::Graphics::ScopedSaveState saved (g);
        g.reduceClipRegion (areas.artwork);
        g.drawImage (backgroundPainting, areas.artwork.toFloat(), juce::RectanglePlacement::fillDestination);
    }
    g.setColour (mutedText);
    g.setFont (controlFont (20.0f * scale));
    const auto hint = areas.gridArea.toFloat().reduced (18.0f * scale, 11.0f * scale)
        .withHeight (40.0f * scale).withTrimmedLeft (570.0f * scale);
    g.drawText ("Click a pad to play  /  MIDI C3-D#6", hint,
                juce::Justification::centredRight, false);
}

void TaikorAudioProcessorEditor::resized()
{
    const auto areas = calculateLayout();
    const auto scale = static_cast<float> (getWidth()) / designWidth;
    const auto px = [scale] (int n) { return juce::roundToInt (static_cast<float> (n) * scale); };
    const auto rect = [&px] (int x, int y, int w, int h)
    {
        return juce::Rectangle<int> (px (x), px (y), px (w), px (h));
    };
    sealLabel.setBounds (rect (28, 16, 64, 64));
    sealLabel.setFont (displayFont (50.0f * scale, juce::Font::bold));
    logoLabel.setBounds (rect (104, 8, 326, 80));
    logoLabel.setFont (displayFont (80.0f * scale, juce::Font::bold));
    editionLabel.setBounds (rect (440, 25, 160, 46));
    editionLabel.setFont (controlFont (18.0f * scale));
    statusDisplay.setBounds (rect (608, 30, 234, 36));
    limiterLabel.setBounds (rect (850, 18, 252, 24));
    limiterLabel.setFont (controlFont (18.0f * scale, true));
    meter.setBounds (rect (850, 44, 252, 34));
    panicButton.setBounds (rect (1130, 28, 124, 40));

    auto grid = areas.gridArea.reduced (px (16), px (6));
    gridCaption.setBounds (grid.removeFromTop (px (40)).withWidth (px (520)));
    gridCaption.setFont (displayFont (34.0f * scale, juce::Font::bold));
    grid.removeFromTop (px (8));
    const int gap = px (7);
    const int rowHeaderWidth = px (280);
    auto headings = grid.removeFromTop (px (27));
    headings.removeFromLeft (rowHeaderWidth + gap);
    const auto columns = taikor::ui::rowLayout (headings.getWidth(), 4, gap, 4);
    for (int col = 0; col < 4; ++col)
    {
        auto& label = strokeLabels[static_cast<std::size_t> (col)];
        label.setBounds (headings.getX() + taikor::ui::cellOffset (columns, gap, col),
                          headings.getY(), columns.cellSize, headings.getHeight());
    }
    const int rowHeight = (grid.getHeight() - gap * 3) / 4;
    for (int row = 0; row < 4; ++row)
    {
        auto line = juce::Rectangle<int> (grid.getX(), grid.getY() + row * (rowHeight + gap),
                                         grid.getWidth(), rowHeight);
        octaveButtons[static_cast<std::size_t> (row)]->setBounds (line.removeFromLeft (rowHeaderWidth));
        line.removeFromLeft (gap);
        const auto cells = taikor::ui::rowLayout (line.getWidth(), 4, gap, 4);
        for (int col = 0; col < 4; ++col)
            pads[static_cast<std::size_t> (row * 4 + col)]->setBounds (
                line.getX() + taikor::ui::cellOffset (cells, gap, col), line.getY(),
                cells.cellSize, line.getHeight());
    }

    auto head = areas.head.reduced (px (10), px (6));
    headCaption.setBounds (head.removeFromTop (px (34)));
    headCaption.setFont (controlFont (22.0f * scale, true));
    head.removeFromTop (px (12));
    headDisplay.setBounds (head);

    struct Cell
    {
        juce::Component* component;
        int span;
    };
    const auto layoutDeck = [&px, scale] (juce::Rectangle<int> area, juce::Label& label,
                                          std::initializer_list<Cell> cells, int columnsCount)
    {
        auto working = area.reduced (px (12), px (6));
        label.setBounds (working.removeFromTop (px (42)).withTrimmedLeft (px (4)));
        label.setFont (displayFont (32.0f * scale, juce::Font::bold));
        working.removeFromTop (px (16));
        int totalColumns = 0;
        for (const auto& cell : cells)
            totalColumns += cell.span;
        const auto rows = (totalColumns + columnsCount - 1) / columnsCount;
        const auto knobRowHeight = (working.getHeight() - px (8) * (rows - 1)) / rows;
        const auto cols = taikor::ui::rowLayout (working.getWidth(), columnsCount, px (5), columnsCount);
        int index = 0;
        for (const auto& cell : cells)
        {
            cell.component->setBounds (
                working.getX() + taikor::ui::cellOffset (cols, px (5), index % columnsCount),
                working.getY() + (index / columnsCount) * (knobRowHeight + px (8)),
                cols.cellSize * cell.span + px (5) * (cell.span - 1), knobRowHeight);
            index += cell.span;
        }
    };
    layoutDeck (areas.drumDeck, drumDeckLabel,
                { { &sizeKnob, 1 }, { &depthKnob, 1 }, { &tensionKnob, 1 },
                  { &headMaterialKnob, 1 }, { &shellMaterialKnob, 1 },
                  { &headDampingKnob, 1 }, { &pitchKnob, 1 },
                  { drumLayoutSwitch.get(), 3 } }, 5);
    layoutDeck (areas.strokeDeck, strokeDeckLabel,
                { { &hardnessKnob, 1 }, { &strikePositionKnob, 1 },
                  { &strikeAzimuthKnob, 1 }, { &velocityDepthKnob, 1 },
                  { &humaniseKnob, 1 }, { &tensionModKnob, 1 },
                  { &strikeNoiseKnob, 1 }, { &ensembleVariationKnob, 1 },
                  { velocityCurveSwitch.get(), 2 } }, 5);
    layoutDeck (areas.microphoneDeck, microphoneDeckLabel,
                { { &micDistanceKnob, 1 }, { &micSpreadKnob, 1 }, { &widthKnob, 1 },
                  { &driveKnob, 1 }, { &outputHighPassKnob, 1 }, { &outputKnob, 1 } }, 6);

    auto switches = areas.switchDeck.reduced (px (14), px (12));
    performerSwitch->setBounds (switches.removeFromTop (px (57)));
    switches.removeFromTop (px (10));
    ensembleSizeSwitch->setBounds (switches);
}

void TaikorAudioProcessorEditor::timerCallback()
{
    statusDisplay.setStatus (audioProcessor.getActiveVoiceCount(),
                             audioProcessor.isEngineReady(),
                             audioProcessor.getCurrentSampleRateForDisplay());
    meter.setLevels (audioProcessor.getOutputLevel (0),
                     audioProcessor.getOutputLevel (1));

    std::array<std::uint32_t, taikor::articulationCount> currentTriggerCounters {};
    bool receivedTrigger = false;
    for (std::size_t artIdx = 0; artIdx < taikor::articulationCount; ++artIdx)
    {
        const auto articulation = static_cast<taikor::Articulation> (artIdx);
        currentTriggerCounters[artIdx] = audioProcessor.getTriggerCounter (articulation);
        receivedTrigger = receivedTrigger
                       || currentTriggerCounters[artIdx] != observedTriggerCounters[artIdx];
    }

    // Acquire the published counters before taking the relaxed visual snapshot:
    // the audio thread writes the strike and octave first, then release-publishes
    // its counter. This ordering makes the snapshot belong to that new trigger.
    taikor::DrumVisualState visual;
    audioProcessor.getVisualState (visual);

    // Host MIDI has no opportunity to click a row header first. Follow the
    // newest sounding drum so its painting, pad flash and measurements all
    // describe the same event. UI-triggered strokes already select their row
    // immediately, so this is a no-op for the on-screen playing surface.
    if (receivedTrigger)
        selectOctave (visual.lastOctaveOffset);

    headDisplay.setStrike (visual.strikeRadius, visual.strikeAngle, visual.strikeLevel,
                           visual.lastArticulation);

    const auto engineParameters = audioProcessor.snapshotEngineParameters();
    headDisplay.setMicrophones (engineParameters.micSpread, engineParameters.micDistance);

    const auto measurements = audioProcessor.measureDrum (selectedOctave);
    headDisplay.setMeasurements (measurements.soundingHz,
                                 measurements.breathingModeHz,
                                 measurements.radiusMetres * 200.0f,
                                 measurements.tailSeconds);

    for (std::size_t padIdx = 0; padIdx < totalPadCount; ++padIdx)
    {
        if (pads[padIdx] == nullptr)
            continue;

        const auto articulation = pads[padIdx]->getArticulation();
        const auto artIdx = static_cast<std::size_t> (articulation);
        const auto counter = currentTriggerCounters[artIdx];

        if (counter != observedTriggerCounters[artIdx])
        {
            if (pads[padIdx]->getOctaveOffset() == selectedOctave)
                pads[padIdx]->triggerFlash();
        }
        else
        {
            pads[padIdx]->advanceFlash();
        }
    }

    observedTriggerCounters = currentTriggerCounters;
}
