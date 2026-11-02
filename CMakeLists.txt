cmake_minimum_required(VERSION 3.18)

# Project configuration ---
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" PROJECT_VERSION_STRING)
string(STRIP "${PROJECT_VERSION_STRING}" PROJECT_VERSION_STRING)

string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+(\\.[0-9]+)?" PROJECT_VERSION_NUMERIC "${PROJECT_VERSION_STRING}")

if (NOT PROJECT_VERSION_NUMERIC)
    message(FATAL_ERROR "Invalid VERSION file format: '${PROJECT_VERSION_STRING}'")
endif()

project(lagrange
    VERSION ${PROJECT_VERSION_NUMERIC}
    LANGUAGES C
)

if (WIN32)
    enable_language(RC)
endif()

# === C++ Standard, compile flags, configs ===
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Release keeps full optimization and also emits debug symbols, so the
# optimized binary can be debugged and profiled. Debug symbols do not change
# the generated code.
if (MSVC)
    # /Z7 stores debug info in each .obj (no shared compile-time PDB);
    # /DEBUG makes the linker write the .pdb, and /OPT:REF /OPT:ICF restore
    # the release link optimizations that /DEBUG turns off.
    set(CMAKE_C_FLAGS_RELEASE "/O2 /Z7")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /DEBUG /OPT:REF /OPT:ICF")
else()
    set(CMAKE_C_FLAGS_RELEASE "-O3 -g")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# for multi-config generators (Visual Studio, Xcode)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

add_compile_definitions(
    $<$<CONFIG:Debug>:DEBUG>
    $<$<CONFIG:Release>:RELEASE>
    _USE_MATH_DEFINES
    _CRT_SECURE_NO_WARNINGS
    _CRT_NONSTDC_NO_DEPRECATE
)

# === Generated headers / resources ===
set(LAGRANGE_GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY ${LAGRANGE_GENERATED_DIR})

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.h.in
    ${LAGRANGE_GENERATED_DIR}/version.h
)

# === Dependencies ===
include(FetchContent)

FetchContent_Declare(
    aether
    GIT_REPOSITORY https://github.com/csooriyakumaran/aether.git
    GIT_TAG v0.0.19
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(aether)

FetchContent_Declare(
    laplace
    SOURCE_DIR "D:/laplace"
    SOURCE_SUBDIR include
    # GIT_REPOSITORY https://github.com/csooriyakumaran/laplace.git
    # GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(laplace)

add_executable(lagrange
    src/main.c
)

# Generate and attach Windows Version resource
if (WIN32)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/resources.rc.in
        ${LAGRANGE_GENERATED_DIR}/resources.rc
        @ONLY
    )

    target_sources(lagrange PRIVATE
        ${LAGRANGE_GENERATED_DIR}/resources.rc
    )
endif()

target_include_directories(lagrange PRIVATE
${laplace_SOURCE_DIR}/include)

target_link_libraries(lagrange PRIVATE aether::aether)

# === Install Rules ===
install(TARGETS lagrange
    RUNTIME DESTINATION bin
    COMPONENT Runtime
)

install(FILES README.md LICENSE VERSION
    DESTINATION .
    COMPONENT Docs
)
