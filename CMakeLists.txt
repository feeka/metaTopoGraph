cmake_minimum_required(VERSION 3.5...4.1)
project(metaTopoGraph)

if (CMAKE_VERSION VERSION_LESS "3.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else ()
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif ()

option(STATIC_BUILD "Build static execution" OFF)
option(SANITIZER    "Enable sanitizers"      OFF)

set(MEGAHIT_SRC "${CMAKE_SOURCE_DIR}/libs/megahit/src")

include_directories(${MEGAHIT_SRC})
include_directories(${CMAKE_SOURCE_DIR}/src)

file(GLOB SDBG_SOURCE "${MEGAHIT_SRC}/sdbg/*.cpp")

find_package(ZLIB REQUIRED)
find_package(OpenMP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -ftemplate-depth=3000")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")

if (OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

if (SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
endif ()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -ggdb -O1")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS}")

set(TOPO_SOURCES
    src/topo_extractor.cpp
    src/topo_json_writer.cpp
    src/main.cpp
)

add_executable(metaTopoGraph ${TOPO_SOURCES} ${SDBG_SOURCE})

target_link_libraries(metaTopoGraph ${ZLIB_LIBRARIES})
if (OpenMP_CXX_FOUND)
    target_link_libraries(metaTopoGraph OpenMP::OpenMP_CXX)
endif ()
