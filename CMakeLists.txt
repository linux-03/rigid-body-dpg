cmake_minimum_required(VERSION 3.18)
project(ASC_ODE_BLA LANGUAGES CXX)

# Ensure modern CMake policy is applied so subprojects (pybind11) don't fail due to policy checks.
cmake_policy(VERSION 3.18)

set (CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# --- CI/pyodide: skip host-Python probing when cross-compiling to Emscripten ---
if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten" OR CMAKE_CROSSCOMPILING)
  message(STATUS "Cross-compiling detected: disabling PYBIND11_FINDPYTHON and clearing Python variables to avoid host-Python vs cross-compiler checks.")

  # Force OFF so pybind11's CMake helpers won't try to find/validate host python
  set(PYBIND11_FINDPYTHON OFF CACHE BOOL "CI: do not find host python when cross-compiling" FORCE)

  # Also clear any Python executable / include / lib variables pybind11 or FindPython might consult.
  # Use CACHE + FORCE so we override any previously set values.
  set(PYBIND11_FINDPYTHON_EXECUTABLE "" CACHE FILEPATH "CI: no host python" FORCE)
  set(Python_EXECUTABLE "" CACHE FILEPATH "CI: no host python" FORCE)
  set(Python_INCLUDE_DIRS "" CACHE STRING "CI: empty" FORCE)
  set(Python_LIBRARIES "" CACHE STRING "CI: empty" FORCE)
endif()
# --- end guard ---

# Ensure FetchContent is available for providing pybind11 as a subproject when needed
include(FetchContent)

# --- pybind11 / Python module setup (works for native SKBUILD and for Pyodide cross-build) ---
if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten" OR CMAKE_CROSSCOMPILING)
  message(STATUS "Cross-compiling for Emscripten: providing pybind11 via FetchContent (no config probing).")
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.10.4
  )
  FetchContent_MakeAvailable(pybind11)

else()
  # Native path: if SKBUILD is used (scikit-build-core) prefer the SKBUILD flow.
  if(DEFINED SKBUILD)
    message(STATUS "Native SKBUILD build: using Python and pybind11 config package")
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 CONFIG REQUIRED)
  else()
    # Native, non-SKBUILD: prefer vendored pybind11 if present, else installed config, else FetchContent.
    if(EXISTS "${CMAKE_SOURCE_DIR}/third_party/pybind11/CMakeLists.txt")
      message(STATUS "Using vendored pybind11 via add_subdirectory")
      add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/pybind11)
    else()
      find_package(pybind11 CONFIG QUIET)
      if(pybind11_FOUND)
        message(STATUS "Found installed pybind11 (config package) - using it for native build")
      else()
        message(STATUS "pybind11 not found installed - fetching via FetchContent")
        FetchContent_Declare(
          pybind11
          GIT_REPOSITORY https://github.com/pybind/pybind11.git
          GIT_TAG        v2.10.4
        )
        FetchContent_MakeAvailable(pybind11)
      endif()
    endif()
  endif()
endif()

# If DEBUG build type, add debug definition
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-DDEBUG_MODE)
endif()

include_directories(bla/src)
include_directories(bla/HPC/src)
include_directories(bla/HPC/concurrentqueue)

# Include subdirectories
add_subdirectory(helper)
add_subdirectory(Elements)
add_subdirectory(nonlinfunc)
add_subdirectory(src)
include_directories(helper)

# Add the main executable
add_executable(main main.cc)

# Link libraries to the main executable
target_link_libraries(main PRIVATE helper Elements nonlinfunc src)

# Regression tests: ctest --test-dir <build> --output-on-failure
enable_testing()
add_executable(test_core tests/test_core.cc)
target_link_libraries(test_core PRIVATE helper Elements nonlinfunc src)
foreach(test_name
    mass_scales_kinematics inertia_torque_free_precession gravity_direction_and_equivalence
    potential_consistent_with_force hidden_constraint_uses_mass anchor_on_body_other_than_zero
    exact_rotations_are_regular left_invariance body_without_beams_runs
    copy_keeps_index_and_constraints fixed_anchor_away_from_origin heavy_top_second_order
    pendulum_energy_bounded jacobian_matches_finite_differences)
  add_test(NAME ${test_name} COMMAND test_core ${test_name})
endforeach()

# pybind11 module — STRICTLY guarded so python_add_library is only used when SKBUILD defines it.
if(DEFINED SKBUILD)
  # scikit-build's helper (python_add_library) is available only when building via SKBUILD.
  # This branch will be used when scikit-build/core drives the build (e.g., pip wheel creation).
  message(STATUS "Using scikit-build helper python_add_library for wheel build")
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)
  python_add_library(rigid_body_FEM MODULE bind.cc WITH_SOABI)
  target_link_libraries(rigid_body_FEM PRIVATE pybind11::headers PUBLIC Elements nonlinfunc src)
  target_compile_definitions(rigid_body_FEM PRIVATE VERSION_INFO=${PROJECT_VERSION})
else()
  # Normal build (native non-SKBUILD or cross-compile): use pybind11_add_module.
  message(STATUS "Using pybind11_add_module for build (non-SKBUILD or cross-compile)")
  pybind11_add_module(rigid_body_FEM bind.cc)
  target_link_libraries(rigid_body_FEM PUBLIC Elements nonlinfunc src)
endif()

# (Optional) any extra packaging / install logic follows...
