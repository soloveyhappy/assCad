cmake_minimum_required(VERSION 2.8.9)
project (root)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(OUTPUT_PATH_SUFFIX "-x64")
else()
    set(OUTPUT_PATH_SUFFIX "-x86")
endif()

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/${CMAKE_C_COMPILER_ID}_${CMAKE_BUILD_TYPE}${OUTPUT_PATH_SUFFIX})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

cmake_policy(SET CMP0020 OLD)
# см. https://github.com/sakra/cotire/issues/81
if(POLICY CMP0058)
    cmake_policy(SET CMP0058 NEW)
endif()
if(POLICY CMP0054)
    cmake_policy(SET CMP0054 NEW)
endif()



add_subdirectory(avoid)
add_subdirectory(main)

add_dependencies(main avoid)
