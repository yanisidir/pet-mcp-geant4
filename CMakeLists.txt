cmake_minimum_required(VERSION 3.16)
project(minimal_geant4_example)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and visualization drivers" ON)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

execute_process(
  COMMAND root-config --incdir
  OUTPUT_VARIABLE ROOT_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND root-config --glibs
  OUTPUT_VARIABLE ROOT_LIBRARIES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
separate_arguments(ROOT_LIBRARIES UNIX_COMMAND "${ROOT_LIBRARIES}")

include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(minimal_geant4 minimal_geant4.cc ${sources} ${headers})
target_include_directories(minimal_geant4 PRIVATE ${ROOT_INCLUDE_DIR})
target_link_libraries(minimal_geant4 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

set(MINIMAL_SCRIPTS
  macros/run.mac
  macros/vis.mac
  macros/vis_file.mac
  macros/vis_global.mac
  macros/vis_plus_z.mac
  macros/vis_minus_z.mac
  macros/vis_face_plus_z.mac
  macros/vis_export.mac
  macros/tree.mac
)

foreach(_script ${MINIMAL_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
  )
endforeach()

install(TARGETS minimal_geant4 DESTINATION bin)
