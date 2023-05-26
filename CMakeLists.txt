cmake_minimum_required(VERSION 3.11)
project(cdt2d VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

option(BUILD_EXECUTABLE "Build cdt2d executable" ON)
option(BUILD_LIBRARY    "Build libcdt2d library" OFF)

# Find Qt
# On macOS with homebrew, need to set Qt5_DIR=/opt/homebrew/opt/qt5/lib/cmake/Qt5
find_package(Qt5 REQUIRED COMPONENTS Widgets Core Gui)

if(BUILD_EXECUTABLE)
  add_executable(cdt2d main.cpp ui.hpp ui.cpp predicates.hpp cdt2d.hpp cdt2d.cpp)
  set_target_properties(cdt2d PROPERTIES OUTPUT_NAME cdt2dviewer)
  target_link_libraries(cdt2d Qt5::Widgets Qt5::Gui Qt5::Core)
endif()

if(BUILD_LIBRARY)
  add_library(cdt2d STATIC cdt2d.cpp)
  target_compile_definitions(cdt2d PRIVATE -DTETLIBRARY)
  # Generate position independent code
  set_target_properties(cdt2d PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
