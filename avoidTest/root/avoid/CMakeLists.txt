cmake_minimum_required(VERSION 2.8.9)

set(PROJECT_NAME avoid)
project(${PROJECT_NAME})

#if(POLICY CMP0020)
#    cmake_policy(SET CMP0020 NEW)
#endif()

#if(POLICY CMP0043)
#  cmake_policy(SET CMP0043 OLD)
#endif()

#if(POLICY CMP0054)
#    cmake_policy(SET CMP0054 NEW)
#endif()


set(LIBAVOID_EXPORTS OFF CACHE BOOL "Build shared library")
set(LIBAVOID_NO_DLL ON CACHE BOOL "No need Dll")


include_directories(${CMAKE_SOURCE_DIR}/avoid/src)

file(GLOB FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/libavoid/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/libavoid/*.h
)

set(PROJECT_FILES
    ${PROJECT_FILES}
    ${FILES}
)

add_library(${PROJECT_NAME} STATIC ${PROJECT_FILES})

#target_compile_definitions(${PROJECT_NAME} PUBLIC _AFXDLL)


#if(LIBAVOID_EXPORTS)
#   add_library(${PROJECT_NAME} SHARED ${PROJECT_FILES})
#   target_compile_definitions(${PROJECT_NAME} PUBLIC LIBAVOID_EXPORTS)
#else()
#    add_library(${PROJECT_NAME} STATIC ${PROJECT_FILES})
   
#endif()



#target_compile_definitions(${PROJECT_NAME} PUBLIC GRAPHICSENGINE_LIBRARY)

#if (NOT MSVC)
#    target_link_libraries(
#        ${PROJECT_NAME}
#    )
#else()
#    # Для студийного компилятора запрещаем использование макросов min, max
#    target_compile_definitions(${PROJECT_NAME} PUBLIC NOMINMAX)
#    target_link_libraries(
#        ${PROJECT_NAME}
#    )
#endif()

