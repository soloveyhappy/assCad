cmake_minimum_required(VERSION 2.8.9)
set(PROJECT_NAME "editor")
project(${PROJECT_NAME})

if(POLICY CMP0020)
    cmake_policy(SET CMP0020 NEW)
endif()

# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(MINGW)
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++14 -Wall")
endif()

# Widgets finds its own dependencies (QtGui and QtCore).
find_package(Qt5Widgets REQUIRED)

# The Qt5Widgets_INCLUDES also includes the include directories for
# dependencies QtCore and QtGui
include_directories(${Qt5Widgets_INCLUDES})

# We need add -DQT_WIDGETS_LIB when using QtWidgets in Qt 5.
#add_definitions(${Qt5Widgets_DEFINITIONS})

#SET( QT_USE_QTOPENGL TRUE )
#SET( QT_USE_QTSQL    TRUE )

# Executables fail to build with Qt 5 in the default configuration
# without -fPIE. We add that here.
#set(CMAKE_CXX_FLAGS "${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

File(GLOB Sources
	${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_executable( ${PROJECT_NAME} ${Sources} )
# The Qt5Widgets_LIBRARIES variable also includes QtGui and QtCore
#target_link_libraries(${PROJECT_NAME} ${Qt5Widgets_LIBRARIES})

qt5_use_modules(${PROJECT_NAME} Core Widgets Gui Sql OpenGl)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)
