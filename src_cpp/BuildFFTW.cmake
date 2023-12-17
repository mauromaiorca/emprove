
IF(UNIX AND NOT APPLE)
set(FFTW_EXTERNAL_PATH "${EXTERNAL_LIBS}/fftw")
ENDIF(UNIX AND NOT APPLE)
IF(APPLE)
set(FFTW_EXTERNAL_PATH "${EXTERNAL_LIBS}/fftw")
ENDIF(APPLE)
message(STATUS " FFTW_EXTERNAL_PATH >>>  " ${FFTW_EXTERNAL_PATH})


if (NOT DEFINED TARGET_X86)
    try_compile(TARGET_X86 ${CMAKE_BINARY_DIR}
                "${EXTERNAL_LIBS}/cmake/TestX86.c")
endif()

find_path(   OWN_FFTW_INCLUDES NAMES fftw3.h PATHS ${FFTW_EXTERNAL_PATH}/include NO_DEFAULT_PATH) 
find_library(OWN_FFTW_DOUBLE   NAMES fftw3   PATHS ${FFTW_EXTERNAL_PATH}/lib     NO_DEFAULT_PATH)

if(OWN_FFTW_INCLUDES AND OWN_FFTW_DOUBLE)
	set(FFTW_FOUND TRUE)
	set(BUILD_OWN_FFTW FALSE)
else()
	message(STATUS "########################################################")
	message(STATUS "########################################################")
	message(STATUS "###        FFTW LIBRARIES NOT FOUND                  ###")
	message(STATUS "###        FFTW WILL BE DOWNLOADED AND               ###")
	message(STATUS "###        BUILT DURING COMPILE-TIME.                ###")
	message(STATUS "###                                                  ###")
	message(STATUS "###        INTERNET CONNECTION IS REQUIRED           ###")
	message(STATUS "########################################################")
	message(STATUS "########################################################")	
	set(FFTW_FOUND FALSE)

	set(ext_conf_flags_fft --enable-shared --prefix=${FFTW_EXTERNAL_PATH})
	if(TARGET_X86)
		set(ext_conf_flags_fft ${ext_conf_flags_fft} --enable-avx)
	endif()
		
	#	set(OWN_FFTW_SINGLE ${FFTW_EXTERNAL_PATH}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3f${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(OWN_FFTW_DOUBLE ${FFTW_EXTERNAL_PATH}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(OWN_FFTW_STATIC ${FFTW_EXTERNAL_PATH}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})

	set(OWN_FFTW_INCLUDES "${FFTW_EXTERNAL_PATH}/include" )
	set(FFTW_PATH ${FFTW_PATH} ${FFTW_EXTERNAL_PATH})

	set(FFTW_EXTERNAL_LIBS_TAR_DIRECTORY  ${FFTW_EXTERNAL_PATH})
	set(FFTW_EXTERNAL_LIBS_EXTRACT_TARGET ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY})

	set(FFTW_FFTW3_TAR_FILE ftp://ftp.fftw.org/pub/fftw/fftw-3.3.10.tar.gz)
	set(FFTW_FFTW3_LIB_DIR ${FFTW_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3)

	message(STATUS "######  Download  ##################")


	include(ExternalProject)
	externalproject_add(own_fftw
	URL ${FFTW_FFTW3_TAR_FILE}
	URL_MD5 8ccbf6a5ea78a16dbc3e1306e234cc5c
	DOWNLOAD_DIR ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY}
	SOURCE_DIR ${FFTW_FFTW3_LIB_DIR}
	CONFIGURE_COMMAND <SOURCE_DIR>/configure ${ext_conf_flags_fft}
	INSTALL_DIR ${FFTW_EXTERNAL_PATH}/fftw3
	BINARY_DIR ${FFTW_EXTERNAL_PATH}/fftw3
	BUILD_COMMAND ${MAKE}
	LOG_INSTALL)
	#	message(STATUS "######  DONE  ##################")
endif()

if (FFTW_DOUBLE_REQUIRED)
	set(FFTW_LIBRARIES ${OWN_FFTW_DOUBLE} ${FFTW_LIBRARIES})
endif()

if (FFTW_INCLUDES)
	set(FFTW_INCLUDES ${OWN_FFTW_INCLUDES} ${FFTW_INCLUDES})
else()
	set(FFTW_INCLUDES ${OWN_FFTW_INCLUDES})
endif()

