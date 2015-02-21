CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(QP_SDP)

INCLUDE_DIRECTORIES( src )
MESSAGE("CMAKE_BINARY_DIR : ${CMAKE_BINARY_DIR}")
MESSAGE("CMAKE_SOURCE_DIR : ${CMAKE_SOURCE_DIR}")

SET(EXECUTABLE_OUTPUT_PATH	${CMAKE_BINARY_DIR})
SET(   LIBRARY_OUTPUT_PATH	${CMAKE_BINARY_DIR})
 
IF( WIN32 )
	SET( CMAKE_CONFIGURATION_TYPES "Debug;Release;" CACHE STRING "limited configs" FORCE)
	SET( CMAKE_CXX_FLAGS_DEBUG          "/D_DEBUG /MTd /Zi /Ob0 /Od /RTC1" )
	SET( CMAKE_CXX_FLAGS_RELEASE        "/MT /O2 /Ob2 /D NDEBUG")
ELSE()
#	SET ( CMAKE_BUILD_TYPE "RELEASE"     )
	SET ( CMAKE_BUILD_TYPE "DEBUG"     )
#	IF( WITH_DEBUG )
#		SET ( CMAKE_BUILD_TYPE "DEBUG"     )
#	ELSE( WITH_DEBUG )
#		SET ( CMAKE_BUILD_TYPE "RELEASE"     )
#	ENDIF( WITH_DEBUG )

	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -fmessage-length=0  -fPIC" )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-conversion -Wno-sign-compare")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-strict-aliasing -Wno-unused-parameter")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOSHORTS -DUSEGETTIME")
ENDIF( WIN32 )
 
#INCLUDE_DIRECTORIES($ENV{CPLEX_INC})
SET(CSDP_LIB	/home/mruiz/Csdp/lib/libsdp.a lapack blas gfortran m)
SET(CSDP_INC	/home/mruiz/Csdp/include)
SET(ASL_LIB	/home/mruiz/asl/sys.x86_64.Linux/amplsolver.a dl)
SET(ASL_INC	/home/mruiz/asl/sys.x86_64.Linux)

INCLUDE_DIRECTORIES( ${CSDP_INC} ${ASL_INC} )


FILE(  GLOB_RECURSE  qp_sdp_files   ./src/*.cpp )
MESSAGE( "qp_sdp_files is : ${qp_sdp_files}" )
ADD_EXECUTABLE( qp_sdp  ${qp_sdp_files} )
TARGET_LINK_LIBRARIES(qp_sdp ${CSDP_LIB} ${ASL_LIB})

