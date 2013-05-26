cmake_minimum_required (VERSION 2.6) 
project (HHParticleFilter) 

if(APPLE)
	set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
	set (CMAKE_CXX_COMPILER "clang++")
	set (CMAKE_C_COMPILER "clang")
endif()

include(cotire)
include(FindPkgConfig)

find_package(OpenCV 2.4 REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wno-sign-compare")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")

file(GLOB MYSOURCES	*.cpp *.c)
file(GLOB MYHEADERS	*.h *.hpp)

IF(APPLE)
   include_directories( /usr/local/include/ )
   INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
   set(FFMPEG_LIBS "/usr/local/lib/libavformat.dylib /usr/local/lib/libavcodec.dylib /usr/local/lib/libavutil.dylib /usr/local/lib/libswscale.dylib")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${FFMPEG_LIBS}")
ENDIF (APPLE)

add_executable(HHParticleFilter ${MYSOURCES} ${MYHEADERS} )

target_link_libraries(HHParticleFilter 
	${OpenCV_LIBS}
	)

set_target_properties(HHParticleFilter PROPERTIES COTIRE_CXX_PREFIX_HEADER_INIT "std.h")
set_target_properties(HHParticleFilter PROPERTIES COMPILE_FLAGS "-DWITHOUT_OPENCL")

cotire(HHParticleFilter)