#Set different variables according to the detected compiler and processor
#based on $CMAKE_CXX_COMPILER_ID it sets the following variables:
# WARNINGS, EXTRA_WARNINGS, EXTRA, OPTIM, ARCH, DEBUG, _VERSION, PROFILING
# It can also edit CMAKE_SHARED_LINKER_FLAGS and CMAKE_EXE_LINKER_FLAGS

INCLUDE("${CMAKE_SOURCE_DIR}/tools/cmake/BuildVersion.cmake")
BuildVersion()

MACRO (SET_COMPILER_OPTIONS)
	###########################################################
	IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
		#SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited configs"  FORCE)
		SET(WARNINGS "/W4 /D_CRT_SECURE_NO_WARNINGS /EHsc") #Za: strict ansi EHsc: handle c++ exceptions
		#SET(EXTRA_WARNINGS "/Wp64") #/Wall
		SET(OPTIM "/O2 /DNDEBUG /MD /DNOSAFECHECKS")
		SET(ARCH_OPTIM "/arch:SSE2")
		SET(ARCH_SAFE "")
		SET(DEBUG "/Z7 /Od /D__DEBUG /MDd")
		SET(_VERSION "/D_VERSION=${_versionString}")
		IF(BUILD_SHARED_LIBS)
			ADD_DEFINITIONS(/DMIO_DLL)
		ENDIF(BUILD_SHARED_LIBS)
		IF(GUI_EXCEPTIONS)
			SET(MSG_BOX "/DMESG_BOX")
		ENDIF(GUI_EXCEPTIONS)
		
	###########################################################
	ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
		SET(WARNINGS "-Wall -Wno-long-long  -Wswitch")
		SET(DEEP_WARNINGS "-Wshadow -Wpointer-arith -Wconversion -Winline -Wdisabled-optimization") #-Wfloat-equal -Wpadded
		SET(EXTRA_WARNINGS "-Wextra -pedantic ${DEEP_WARNINGS}")
		SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS")
		IF("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64" OR "${CMAKE_SYSTEM_PROCESSOR}" MATCHES "AMD64")
			SET(ARCH_SAFE  "-march=nocona -mtune=nocona")
		ENDIF()
		SET(DEBUG "-g3 -O0 -D__DEBUG")
		SET(_VERSION "-D_VERSION=${_versionString}")
		
	###########################################################
	ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")
		SET(WARNINGS "-Wall -Wno-long-long  -Wswitch")
		SET(DEEP_WARNINGS "-Wshadow -Wpointer-arith -Wconversion -Winline -Wdisabled-optimization") #-Wfloat-equal -Wpadded
		SET(EXTRA_WARNINGS "-Wextra -pedantic ${DEEP_WARNINGS}")
		SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS")
		IF("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64" OR "${CMAKE_SYSTEM_PROCESSOR}" MATCHES "AMD64")
			SET(ARCH_SAFE  "-march=nocona -mtune=nocona")
		ENDIF()
		SET(DEBUG "-g3 -O0 -D__DEBUG")
		SET(_VERSION "-D_VERSION=${_versionString}")
		
	###########################################################
	ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		#we consider that all other compilers support "-" options and silently ignore what they don't know
		SET(WARNINGS "-Wall -Wno-long-long  -Wswitch")
		SET(DEEP_WARNINGS "-Wunused-value -Wshadow -Wpointer-arith -Wconversion -Winline -Wdisabled-optimization -Wctor-dtor-privacy") #-Wfloat-equal -Wpadded
		SET(EXTRA_WARNINGS "-Wextra -pedantic -Weffc++ ${DEEP_WARNINGS}")
		SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS")
		IF("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64" OR "${CMAKE_SYSTEM_PROCESSOR}" MATCHES "AMD64")
			SET(ARCH_SAFE  "-march=nocona -mtune=nocona")
		ENDIF()
		SET(DEBUG "-g3 -O0 -D__DEBUG")
		SET(_VERSION "-D_VERSION=${_versionString}")
		#IF(BUILD_SHARED_LIBS)
		#	ADD_DEFINITIONS(-DMIO_DLL)
		#ENDIF(BUILD_SHARED_LIBS)
		
		SET(PROFILING "-pg -fprofile-arcs") #add ${PROFILING} to the CFLAGS when necessary
		SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -Wunsafe-loop-optimizations -Wvector-operation-performance -Wwrite-strings")
		IF(NOT ANDROID)
			SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -ansi")
			IF(WIN32) #for gcc on windows
				SET(CMAKE_SHARED_LINKER_FLAGS "--enable-auto-import")
				SET(CMAKE_EXE_LINKER_FLAGS "--enable-auto-import")
			ENDIF(WIN32)
		ENDIF(NOT ANDROID)
		EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
		IF(GCC_VERSION VERSION_GREATER 4.2 OR GCC_VERSION VERSION_EQUAL 4.2)
			SET(ARCH_OPTIM  "-march=native -mtune=native")
		ENDIF()
		IF(GCC_VERSION VERSION_GREATER 4.8 OR GCC_VERSION VERSION_EQUAL 4.8)
			IF(NOT WIN32)
				SET(OPTIM "${OPTIM} -flto") #for gcc>4.5, but first implementations were slow, so it is safe to enforce 4.8
			ENDIF(NOT WIN32)
			# if set to ON, all binaries depending on the library have to be compiled the same way.
			#Then, do an "export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.4" and run with "ASAN_OPTIONS=symbolize=1"
			SET(LEAKS_CHECK OFF CACHE BOOL "Set to ON to dynamically check for memory corruption (and do the same for applications linked with MeteoIO)")
			IF (LEAKS_CHECK)
				SET(EXTRA "${EXTRA} -fsanitize=address -fno-omit-frame-pointer")
			ENDIF(LEAKS_CHECK)
		ENDIF()
	
	###########################################################
	ELSEIF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
		SET(WARNINGS_OFF "-Wno-long-long -Wno-date-time -Wno-float-equal -Wno-documentation -Wno-documentation-unknown-command -Wno-old-style-cast -Wno-padded -Wno-missing-noreturn -Wno-weak-vtables -Wno-switch-enum -Wno-covered-switch-default -Wno-global-constructors -Wno-exit-time-destructors -Wno-unknown-pragmas")
		SET(WARNINGS "-Wall -Wswitch -Weverything ${WARNINGS_OFF}") #obviously, we should try to fix the warnings! Keeping in mind that some of these W are half buggy...
		SET(DEEP_WARNINGS "-Wunused-value -Wshadow -Wpointer-arith -Wconversion -Winline -Wdisabled-optimization -Wctor-dtor-privacy") #-Rpass=.* for static analysis
		SET(EXTRA_WARNINGS "-Wextra -pedantic -Weffc++ ${DEEP_WARNINGS}")
		SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS -flto")
		IF("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64" OR "${CMAKE_SYSTEM_PROCESSOR}" MATCHES "AMD64")
			SET(ARCH_SAFE  "-march=nocona -mtune=nocona")
		ENDIF()
		SET(DEBUG "-g3 -O0 -D__DEBUG")
		SET(_VERSION "-D_VERSION=${_versionString}")
		
		SET(PROFILING "-pg") #add ${PROFILING} to the CFLAGS when necessary
		SET(EXTRA "${EXTRA} -fcolor-diagnostics") #-fapple-pragma-pack does not seems necessary; -ftrapv should be replaced by sanitize=integer
		SET(LEAKS_CHECK OFF CACHE BOOL "Set to ON to dynamically check for memory corruption (and do the same for applications linked with MeteoIO)")
			IF (LEAKS_CHECK)
				SET(EXTRA "${EXTRA} -ftrapv -fno-omit-frame-pointer") #-fsanitize=address,undefined,integer,undefined-trap but this is currently not supported by Apple
			ENDIF(LEAKS_CHECK)
		SET(ARCH_OPTIM  "-march=native")
	ENDIF()
	
	###########################################################
	#targets providing SETs of compiler options
	IF(NOT DEST)
		SET(DEST "safe" CACHE STRING "Choose safe or optimized" FORCE)
	ENDIF(NOT DEST)

	IF (DEST STREQUAL "safe")
		SET(ARCH  "${ARCH_SAFE}")
	ENDIF(DEST STREQUAL "safe")

	IF(DEST STREQUAL "optimized")
		SET(ARCH  "${ARCH_OPTIM}")
	ENDIF(DEST STREQUAL "optimized")

	#show exception messages in a graphical message box
	SET(GUI_EXCEPTIONS OFF CACHE BOOL "Show a message box with exceptions texts ON or OFF")

ENDMACRO (SET_COMPILER_OPTIONS)
