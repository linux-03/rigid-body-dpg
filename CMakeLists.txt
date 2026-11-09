cmake_minimum_required(VERSION 3.18)
project(ASC_ODE_BLA LANGUAGES CXX)

# Ensure modern CMake policy is applied so subprojects (pybind11) don't fail due to policy checks.
cmake_policy(VERSION 3.18)

set (CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# pybind11: wheel builds (native and Pyodide) get it from pyproject.toml.
if(DEFINED SKBUILD)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)
else()
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG        v2.12.0
    )
    FetchContent_MakeAvailable(pybind11)
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

# Native programs; a wheel build only needs the Python module.
if(NOT DEFINED SKBUILD)
  add_executable(main main.cc)
  target_link_libraries(main PRIVATE helper Elements nonlinfunc src)
endif()

if(NOT DEFINED SKBUILD)
  # Regression tests: ctest --test-dir <build> --output-on-failure
  enable_testing()
  add_executable(test_core tests/test_core.cc)
  target_link_libraries(test_core PRIVATE helper Elements nonlinfunc src)
  foreach(test_name
      mass_scales_kinematics inertia_torque_free_precession gravity_direction_and_equivalence
      potential_consistent_with_force hidden_constraint_uses_mass anchor_on_body_other_than_zero
      exact_rotations_are_regular left_invariance body_without_beams_runs
      copy_keeps_index_and_constraints fixed_anchor_away_from_origin heavy_top_second_order
      pendulum_energy_bounded jacobian_matches_finite_differences
      free_body_invariants_are_tight spring_matches_exact_oscillator spring_default_rest_length
      spring_force_is_minus_gradient spring_pair_conserves_momentum
      spring_jacobian_matches_finite_differences coupled_jacobian_matches_finite_differences
      symplectic_body_on_spring
      symplectic_check_detects_damping symplectic_heavy_top)
    add_test(NAME ${test_name} COMMAND test_core ${test_name})
  endforeach()
endif()

# Python module
if(DEFINED SKBUILD)
  python_add_library(rigid_body_FEM MODULE bind.cc WITH_SOABI)
  target_link_libraries(rigid_body_FEM PRIVATE pybind11::headers PUBLIC Elements nonlinfunc src)
  # python/lib_rigid_body is added to the wheel by scikit-build-core.
  install(TARGETS rigid_body_FEM DESTINATION lib_rigid_body)
else()
  pybind11_add_module(rigid_body_FEM bind.cc)
  target_link_libraries(rigid_body_FEM PUBLIC Elements nonlinfunc src)
endif()
