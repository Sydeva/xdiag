cmake_minimum_required(VERSION 3.19)
project(cantedAFM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Select which application source backs the 'main' executable.
set(CANTED_AFM_APP "GSwPG" CACHE STRING "Application source for the main executable (GSwPG or strucfac)")
set_property(CACHE CANTED_AFM_APP PROPERTY STRINGS GSwPG strucfac)

# Select which GSwPG entry file backs the executable when CANTED_AFM_APP=GSwPG.
set(CANTED_AFM_GSWPG_ENTRY "main" CACHE STRING "GSwPG entry source (main or mainani)")
set_property(CACHE CANTED_AFM_GSWPG_ENTRY PROPERTY STRINGS main mainani)

if(NOT CANTED_AFM_APP STREQUAL "GSwPG" AND NOT CANTED_AFM_APP STREQUAL "strucfac")
  message(FATAL_ERROR "CANTED_AFM_APP must be one of: GSwPG, strucfac")
endif()

if(NOT CANTED_AFM_GSWPG_ENTRY STREQUAL "main" AND NOT CANTED_AFM_GSWPG_ENTRY STREQUAL "mainani")
  message(FATAL_ERROR "CANTED_AFM_GSWPG_ENTRY must be one of: main, mainani")
endif()

if(CANTED_AFM_APP STREQUAL "GSwPG")
  set(CANTED_AFM_MAIN_SOURCE "${CANTED_AFM_APP}/${CANTED_AFM_GSWPG_ENTRY}.cpp")
else()
  if(NOT CANTED_AFM_GSWPG_ENTRY STREQUAL "main")
    message(FATAL_ERROR "CANTED_AFM_GSWPG_ENTRY=${CANTED_AFM_GSWPG_ENTRY} is only valid with CANTED_AFM_APP=GSwPG")
  endif()
  set(CANTED_AFM_MAIN_SOURCE "${CANTED_AFM_APP}/main.cpp")
endif()

# Prefer the in-repo native install prefix used by the native canted build.
find_package(xdiag REQUIRED HINTS
  "${CMAKE_CURRENT_SOURCE_DIR}/../install-native"
  "${CMAKE_CURRENT_SOURCE_DIR}/../install"
)

add_executable(main "${CANTED_AFM_MAIN_SOURCE}")
target_link_libraries(main PRIVATE xdiag::xdiag)

