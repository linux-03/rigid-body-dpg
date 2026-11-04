add_library(src STATIC  
    rigid_body_system.cc 
    equation_rb.cc 
    equation_rbs.cc
)
target_include_directories(src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(src PRIVATE Elements nonlinfunc helper pybind11::pybind11)

if(NOT DEFINED SKBUILD)
  add_executable(main_src main.cc)
  target_link_libraries(main_src src Elements nonlinfunc helper)
endif()

