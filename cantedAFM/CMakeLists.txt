cmake_minimum_required(VERSION 3.19)
project(cantedAFM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Prefer the in-repo install prefix used by xdiag's default cmake setup.
find_package(xdiag REQUIRED HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../install")

add_executable(main main.cpp)
target_link_libraries(main PRIVATE xdiag::xdiag)

