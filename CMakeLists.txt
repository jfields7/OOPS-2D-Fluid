cmake_minimum_required(VERSION 3.0)
project(Fluid)

option(DEBUG_PRIMITIVE_SOLVER "Print debug messages in the primitive solver." OFF)
option(DEBUG_RHS "Print debug messages in the righthand side." OFF)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/fluidconfig.h.in 
               ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidconfig.h)

# Add the NumTools package
#find_package(NumTools REQUIRED)

# Generate the Parameters and ParamParser files.
set(PARAM_SRC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fluidgridparser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fluidgridparameters.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fluidioparser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fluidioparameters.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fluidmethodparser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fluidmethodparameters.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fluidicparser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fluidicparameters.cpp
   )
set(PARAM_INC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidgridparameters.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidgridparser.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidioparser.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidioparameters.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidmethodparser.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidmethodparameters.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidicparser.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidicparameters.h
   )
set(SETUP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/scripts/fluid.json)

add_custom_command(
  OUTPUT ${PARAM_INC}
         ${PARAM_SRC}
  DEPENDS ${SETUP_SRC}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/genParams.py ${SETUP_SRC}
  COMMENT "Generating custom Parameters files"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM USES_TERMINAL
)

set(FLUID_INCLUDE_FILES
    include/fluid.h
    include/geom/metric.h
    include/geom/cartesianmetric.h
    include/geom/polarmetric.h
    include/geom/cylindricalmetric.h
    include/geom/sphericalmetric.h
    include/riemann/riemannsolver.h
    include/riemann/hllesolver.h
    include/primitive/primitivesolver.h
    include/primitive/idealsolver.h
    include/primitive/enthalpysolver.h
    include/recon/recon.h
    include/recon/minmod.h
    include/recon/norecon.h
    include/recon/mp5.h
   )
set(FLUID_SOURCE_FILES
    src/main.cpp
    src/fluid.cpp
    src/recon/minmod.cpp
    src/recon/norecon.cpp
    src/recon/mp5.cpp
    src/riemann/riemannsolver.cpp
    src/riemann/hllesolver.cpp
    src/primitive/primitivesolver.cpp
    src/primitive/idealsolver.cpp
    src/primitive/enthalpysolver.cpp
   )

set(SOURCE_FILES ${FLUID_INCLUDE_FILES} ${FLUID_SOURCE_FILES} ${PARAM_INC} ${PARAM_SRC})
add_executable(Fluid ${SOURCE_FILES})
target_include_directories(Fluid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(Fluid PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(Fluid oops2d NumTools ${EXTRA_LIBS})
