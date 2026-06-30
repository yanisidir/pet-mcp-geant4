cmake_minimum_required(VERSION 3.16)

project(minimal_geant4_example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Geant4 REQUIRED ui_all vis_all)
find_package(ROOT QUIET COMPONENTS Core RIO Tree)

if(ROOT_FOUND)
  set(PROJECT_ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS})
  set(PROJECT_ROOT_LIBRARIES ${ROOT_LIBRARIES})
else()
  find_program(ROOT_CONFIG_EXECUTABLE root-config REQUIRED)
  execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --cflags
    OUTPUT_VARIABLE ROOT_CONFIG_CFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs
    OUTPUT_VARIABLE ROOT_CONFIG_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(PROJECT_ROOT_CFLAGS UNIX_COMMAND
    "${ROOT_CONFIG_CFLAGS}")
  separate_arguments(PROJECT_ROOT_LIBRARIES UNIX_COMMAND
    "${ROOT_CONFIG_LIBS}")
endif()

include(${Geant4_USE_FILE})

file(GLOB PROJECT_SOURCES CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/src/*.cc")

add_executable(minimal_geant4
  minimal_geant4.cc
  ${PROJECT_SOURCES})

target_include_directories(minimal_geant4 PRIVATE
  "${PROJECT_SOURCE_DIR}/include"
  ${PROJECT_ROOT_INCLUDE_DIRS})

if(PROJECT_ROOT_CFLAGS)
  target_compile_options(minimal_geant4 PRIVATE ${PROJECT_ROOT_CFLAGS})
endif()

target_link_libraries(minimal_geant4 PRIVATE
  ${Geant4_LIBRARIES}
  ${PROJECT_ROOT_LIBRARIES})

file(GLOB MACRO_FILES CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/macros/*.mac")

file(COPY ${MACRO_FILES}
  DESTINATION "${PROJECT_BINARY_DIR}/macros")
