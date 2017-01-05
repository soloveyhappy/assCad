cmake_minimum_required(VERSION 2.8.9)

set(PROJECT_NAME main)
project(${PROJECT_NAME})


include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/avoid/src)

file(GLOB FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/*.h
)

set(PROJECT_FILES
    ${PROJECT_FILES}
    ${FILES}
)

ADD_EXECUTABLE(${PROJECT_NAME}
    ${PROJECT_FILES}
)

target_link_libraries(
        ${PROJECT_NAME}
        avoid
        )

#add_dependencies(main avoid)

