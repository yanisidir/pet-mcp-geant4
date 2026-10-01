cmake_minimum_required(VERSION 3.16)
project(pet_mcp_geant4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET COMPONENTS Core RIO Tree)
if(NOT ROOT_FOUND)
  # ROOT installe hors du chemin de recherche CMake (snap, /opt, ...).
  find_program(ROOT_CONFIG_EXECUTABLE root-config REQUIRED)
  execute_process(COMMAND "${ROOT_CONFIG_EXECUTABLE}" --incdir
    OUTPUT_VARIABLE ROOT_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND "${ROOT_CONFIG_EXECUTABLE}" --libs
    OUTPUT_VARIABLE root_config_libs OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(ROOT_LIBRARIES UNIX_COMMAND "${root_config_libs}")
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(PET_MCP main.cc ${sources} ${headers})

target_include_directories(PET_MCP PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${ROOT_INCLUDE_DIRS})

target_link_libraries(PET_MCP PRIVATE
  ${Geant4_LIBRARIES}
  ${ROOT_LIBRARIES})

# Les macros sont copiees dans le repertoire de build, comme dans les
# exemples Geant4, pour pouvoir lancer l'application depuis celui-ci.
file(GLOB macros ${PROJECT_SOURCE_DIR}/macros/*.mac)
foreach(macro ${macros})
  get_filename_component(name ${macro} NAME)
  configure_file(${macro} ${PROJECT_BINARY_DIR}/macros/${name} COPYONLY)
endforeach()
