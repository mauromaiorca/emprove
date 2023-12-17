#
# scor'EM
#
#

PROJECT(emprove)
 cmake_minimum_required(VERSION 2.6)

SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread")

SET(BUILD_SHARED_LIBS OFF)
SET(ITK_DYNAMIC_LOADING OFF)

SET(SRC_DIR "./src_cpp/")
SET(BIN_DIR "./app_bin/")

get_filename_component(EXTERNAL_LIBS_ABSOLUTE_PATH "./external_cpp_libs" ABSOLUTE)
SET( EXTERNAL_LIBS  "${EXTERNAL_LIBS_ABSOLUTE_PATH}")
file(MAKE_DIRECTORY "${BIN_DIR}")
file(MAKE_DIRECTORY "${EXTERNAL_LIBS}")
set(FFTW_FFTW3_LIB_DIR ${EXTERNAL_LIBS}/fftw3)

#SET( GCC_COVERAGE_COMPILE_FLAGS  "${GCC_COVERAGE_COMPILE_FLAGS} -O2" )
SET( GCC_COVERAGE_COMPILE_FLAGS  "${GCC_COVERAGE_COMPILE_FLAGS} " )
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )



set(FFTW_DOUBLE_REQUIRED TRUE)
set(FORCE_OWN_FFTW ON)
if(AMDFFTW)
        set(FORCE_OWN_FFTW ON)
endif()

set(FFTW_FOUND TRUE)
set(OWN_FFTW_STATIC ${FFTW_FFTW3_LIB_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})
set(FFTW_INCLUDE ${FFTW_FFTW3_LIB_DIR}/include/)


if(NOT EXISTS ${OWN_FFTW_STATIC} OR NOT EXISTS ${FFTW_INCLUDE})
	message(STATUS "########################################################")
	message(STATUS "########################################################")
	message(STATUS "###        FFTW LIBRARIES NOT FOUND                  ###")
	message(STATUS "###        FFTW WILL BE DOWNLOADED AND               ###")
	message(STATUS "###        BUILT DURING COMPILE-TIME.                ###")
	message(STATUS "###                                                  ###")
	message(STATUS "###        INTERNET CONNECTION IS REQUIRED           ###")
	message(STATUS "########################################################")
	message(STATUS "########################################################")	
	
	set(FFTW_FFTW3_TAR_FILE ftp://ftp.fftw.org/pub/fftw/fftw-3.3.10.tar.gz)
	set(FFTW_EXTERNAL_PATH ${FFTW_FFTW3_LIB_DIR})
	set(ext_conf_flags_fft --enable-shared --prefix=${FFTW_EXTERNAL_PATH})
	message(STATUS "######  Download  ##################")
	message(STATUS "FFTW_FFTW3_TAR_FILE = " ${FFTW_FFTW3_TAR_FILE})
	message(STATUS "FFTW_EXTERNAL_LIBS_TAR_DIRECTORY =" ${FFTW_FFTW3_TAR_FILE})
	message(STATUS "FFTW_FFTW3_LIB_DIR = " ${FFTW_FFTW3_LIB_DIR})
	message(STATUS "ext_conf_flags_fft = " ${ext_conf_flags_fft})
	message(STATUS "FFTW_EXTERNAL_PATH = " ${FFTW_EXTERNAL_PATH})
	message(STATUS "MAKE =" ${MAKE})
	
	include(ExternalProject)
	externalproject_add(own_fftw
	URL ${FFTW_FFTW3_TAR_FILE}
	URL_MD5 8ccbf6a5ea78a16dbc3e1306e234cc5c
	DOWNLOAD_DIR ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY}
	SOURCE_DIR ${FFTW_FFTW3_LIB_DIR}
	CONFIGURE_COMMAND <SOURCE_DIR>/configure ${ext_conf_flags_fft}
	INSTALL_DIR ${FFTW_EXTERNAL_PATH}
	BINARY_DIR ${FFTW_EXTERNAL_PATH}
	BUILD_COMMAND ${MAKE}
	LOG_INSTALL)
else()
	message(STATUS "########################################################")
	message(STATUS "########################################################")
	message(STATUS "###        FFTW LIBRARIES  FOUND ")
	message(STATUS "###        FFTW LIBRARIES DIR = " ${FFTW_FFTW3_LIB_DIR})	
	message(STATUS "###        OWN_FFTW_STATIC    = " ${OWN_FFTW_STATIC})
	message(STATUS "###        FFTW_INCLUDE       = " ${FFTW_INCLUDE})
	message(STATUS "########################################################")
	message(STATUS "########################################################")	
endif()


INCLUDE_DIRECTORIES(
   ${SRC_DIR}
)

if(EXISTS ${OWN_FFTW_STATIC} AND EXISTS ${FFTW_INCLUDE})
#	ADD_EXECUTABLE(${BIN_DIR}scorem_app_assessParticles ${SRC_DIR}scorem_app_assessParticles.cxx)
#	TARGET_LINK_LIBRARIES(${BIN_DIR}scorem_app_assessParticles ${OWN_FFTW_STATIC})

	ADD_EXECUTABLE(${BIN_DIR}emprove_app_starProcess ${SRC_DIR}emprove_app_starProcess.cxx)
	TARGET_LINK_LIBRARIES(${BIN_DIR}emprove_app_starProcess ${OWN_FFTW_STATIC})

	ADD_EXECUTABLE(${BIN_DIR}emprove_app_meanMinMax ${SRC_DIR}emprove_app_meanMinMax.cxx)
	TARGET_LINK_LIBRARIES(${BIN_DIR}emprove_app_meanMinMax )

	ADD_EXECUTABLE(${BIN_DIR}emprove_app_starSimulation ${SRC_DIR}emprove_app_starSimulation.cxx)
	TARGET_LINK_LIBRARIES(${BIN_DIR}emprove_app_starSimulation )
	
#	IF(UNIX)
#		SET (CMAKE_INSTALL_PREFIX /usr/local/bin) #in order to support the make install option
#		INSTALL(TARGETS ${BIN_DIR}scorem_app_assessParticles RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
#	ENDIF(UNIX)
endif()


################################################
#
#  Install into the system
#
################################################
#IF(UNIX)
#	SET (CMAKE_INSTALL_PREFIX /usr/local/bin) #in order to support the make install option
#	INSTALL(TARGETS ${BIN_DIR}scorem_app_assessParticles RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
#ENDIF(UNIX)

