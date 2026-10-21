# Alternative drum engines for the engine-candidate listening test, and the
# tool that renders them against the shipping engine. Built only with
# -DTAIKOR_BUILD_ENGINE_CANDIDATES=ON; nothing here is part of the plug-in.
#
# Each engine owns its folder and is its own library; every .cpp in a folder
# is compiled, so an engine can add files without touching this list, and one
# engine's work in progress never breaks another's build. Each engine also has
# its own copy of the render tool and its own verification suite.

add_library(TaikorCandidateCommon STATIC ShippingSpec.cpp)
target_include_directories(TaikorCandidateCommon PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(TaikorCandidateCommon PUBLIC TaikorDSP)
target_compile_features(TaikorCandidateCommon PUBLIC cxx_std_20)

foreach(engine Fdtd VonKarman Modal)
    file(GLOB engine_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${engine}/*.cpp")
    add_library(TaikorCandidate${engine} STATIC ${engine_sources})
    target_include_directories(TaikorCandidate${engine} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/${engine}")
    target_link_libraries(TaikorCandidate${engine} PUBLIC TaikorCandidateCommon)
    if(NOT MSVC)
        target_compile_options(TaikorCandidate${engine} PRIVATE -Wall -Wextra)
    endif()
endforeach()

# Render tools: all engines (the one the listening set is rendered with), and
# one per engine for development.
function(taikor_candidate_renderer name)
    add_executable(${name} RenderEngineCandidates.cpp)
    foreach(engine ${ARGN})
        target_link_libraries(${name} PRIVATE TaikorCandidate${engine})
        string(TOUPPER ${engine} upper)
        target_compile_definitions(${name} PRIVATE TAIKOR_WITH_${upper}=1)
    endforeach()
    target_link_libraries(${name} PRIVATE TaikorCandidateCommon)
endfunction()
taikor_candidate_renderer(TaikorRenderEngineCandidates Fdtd VonKarman Modal)
taikor_candidate_renderer(TaikorRenderCandidateFdtd Fdtd)
taikor_candidate_renderer(TaikorRenderCandidateFdtdVonKarman Fdtd VonKarman)
taikor_candidate_renderer(TaikorRenderCandidateModal Modal)

# Verification suites: every .cpp in tests/<suite>/ forms one executable.
function(taikor_candidate_suite suite)
    file(GLOB suite_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/${suite}/*.cpp")
    add_executable(TaikorCandidateTests_${suite} ${suite_sources})
    target_link_libraries(TaikorCandidateTests_${suite} PRIVATE TaikorCandidateCommon)
    foreach(engine ${ARGN})
        target_link_libraries(TaikorCandidateTests_${suite} PRIVATE TaikorCandidate${engine})
    endforeach()
endfunction()
taikor_candidate_suite(fdtd Fdtd)
taikor_candidate_suite(vonkarman VonKarman)
taikor_candidate_suite(modal Modal)
taikor_candidate_suite(integration Fdtd VonKarman)
