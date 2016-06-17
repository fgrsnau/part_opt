CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)
SET(CMAKE_CONFIGURATION_TYPES "RelWithDebInfo;Debug" CACHE STRING "RelWithDebInfo;Debug" FORCE)
SET(CMAKE_BUILD_TYPE "RelWithDebInfo")

project (part_opt)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

set (VERSION_MAJOR 1)
set (VERSION_MINOR 3)

SET (SRC ../../code)
SET (MY_LIB_PATH ${CMAKE_SOURCE_DIR}/../../lib/)
SET (MY_BIN_PATH ${CMAKE_SOURCE_DIR}/../../bin/${PROJECT_NAME}/)
SET (MEX_INSTALL_PATH ${CMAKE_SOURCE_DIR}/../../matlab/${PROJECT_NAME}/)

# let these be placed wherewe cmake puts object files, will install then to the proper location
#set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../../tmp/${PROJECT_NAME})
#set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../../tmp/${PROJECT_NAME})
#set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../../tmp/${PROJECT_NAME})


set(BUILD_SHARED_LIBS OFF)

#configure_file (
#  "cfg.h.in"
#  "${PROJECT_BINARY_DIR}/cfg.h"
#  )

# C++11 support.
include(cmake_cpp11/EnableCPP11.cmake)

# Multithreading using OpenMP
# #OPTION(WITH_OPENMP "Enable multi-threading (requires OpenMP)" ON)
SET(WITH_OPENMP ON)
OPTION(WITH_OPENGM "Include test of OpenGM interface" OFF)
# OPTION(MAXFLOW_OPENGM "Use maxflow patch from OpenGM" OFF)
SET(MAXFLOW_OPENGM OFF)
OPTION(WITH_MATLAB "Build mex files for MATLAB" ON)

# Library dependencies
SET (LIBRARY_DEPENDENCIES)

#--------------------------------------------------------------
if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
     message(STATUS "Target is 64 bits")
     if (WIN32)
         set(WINXXBITS Win64)
	 set(MEXEXT mexw64)
     else()
	set(MEXEXT mexa64)
     endif(WIN32)
 else("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
     message(STATUS "Target is 32 bits")
     if (WIN32)
         set(WINXXBITS Win32)
	 set(MEXEXT mexw32)
     else()
	set(MEXEXT mexa32)
     endif(WIN32)
endif("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
#--------------------------------------------------------------


#--------------------------------------------------------------
# MATLAB
#--------------------------------------------------------------
if(WITH_MATLAB)
  message(STATUS "build with MatLab support")
# now matlab does not officially support some compilers, opt to ignore this

include(FindPackageHandleStandardArgs)
include(OpenCVFindMatlab.cmake)

message(STATUS "MATLAB_ROOT_DIR" ${MATLAB_ROOT_DIR})

#find_package(MATLAB REQUIRED)

message(STATUS "MATLAB_LIBRARY_DIRS:" ${MATLAB_LIBRARY_DIRS})
message(STATUS "MATLAB_LIBRARIES:" ${MATLAB_LIBRARIES})


#  message(STATUS "Matlab include dir: ${MATLAB_INCLUDE_DIR}")
#  message(STATUS "Matlab libmex: ${MATLAB_MEX_LIBRARY}")
#  message(STATUS "Matlab libeng: ${MATLAB_ENG_LIBRARY}")
#  message(STATUS "Matlab libmx:  ${MATLAB_MX_LIBRARY}")

#if(WIN32)
#	if(WINXXBITS STREQUAL Win32)
#		SET(MATLAB "C:/Program Files (x86)/MATLAB/R2009a")
#		SET(MATLAB_LIBRARY_PATH ${MATLAB}/extern/lib/win32/microsoft)
#	else()
#		SET(MATLAB "C:/Program Files/MATLAB/R2012b")
#		SET(MATLAB_LIBRARY_PATH ${MATLAB}/extern/lib/win64/microsoft)
#	endif()
#else()
#	#SET(MATLAB "/usr/local/matlab8")
#	SET(MATLAB "/usr/local/MATLAB/R2011a")
#	SET(MATLAB_LIBRARY_PATH ${MATLAB}/bin/glnxa64)
#endif()
#	find_path(MATLAB_INCLUDE_DIR NAMES mex.h PATHS ${MATLAB}/extern/include/)
#	find_library(MATLAB_MEX_LIBRARY NAMES mex PATHS ${MATLAB_LIBRARY_PATH})
#	find_library(MATLAB_MX_LIBRARY NAMES mx PATHS ${MATLAB_LIBRARY_PATH})
#	find_library(MATLAB_ENG_LIBRARY NAMES eng PATHS ${MATLAB_LIBRARY_PATH})
#

endif(WITH_MATLAB)


#--------------------------------------------------------------
# OPENGM
#--------------------------------------------------------------
if(WITH_OPENGM)
#   message(STATUS "build with OPENGM support")
#   find_package(opengm REQUIRED CONFIG)
    SET(OPENGM_INCLUDE_DIR "Z:/work/opengm/include" CACHE STRING "Z:/work/opengm/include")
    include_directories(${OPENGM_INCLUDE_DIR})
    add_definitions(-DWITH_OPENGM)
endif()

if(WITH_MATLAB)

endif()

#--------------------------------------------------------------
# HDF5
#--------------------------------------------------------------
if(WITH_OPENGM)
   message(STATUS "build with HDF5 support")
   if(WIN32)
     # FindHDF5 is broken on Windows
     # it wrongly discovers the dlls instead of the lib files
     # see: www.cmake.org/Bug/bug_relationship_graph.php?bug_id=14111
     # therefore we set the variable manually
     find_library( HDF5_IMPORT_LIB NAMES hdf5 hdf5dll )
     find_path(HDF5_INCLUDE_DIR hdf5.h)
     set(HDF5_LIBRARIES ${HDF5_IMPORT_LIB})
     message(STATUS "HDF5: autodiscovery is broken on WIN32; using global search paths")
     message(STATUS "HDF5: found ${HDF5_LIBRARIES}")
     message(STATUS "HDF5_LIBRARIES:" ${HDF5_LIBRARIES})
   else()
     IF(APPLE)
	find_library(HDF5_IMPORT_LIB NAMES hdf5)
	find_path(HDF5_INCLUDE_DIR hdf5.h)
	set(HDF5_LIBRARIES ${HDF5_IMPORT_LIB})
     ELSE()	
	find_path(HDF5_INCLUDE_DIR hdf5.h)
	message(STATUS "HDF5_INCLUDE_DIR: ${HDF5_INCLUDE_DIR}")
	find_package(HDF5 REQUIRED)
     ENDIF()
	message(STATUS "HDF5_LIBRARIES:" ${HDF5_LIBRARIES})
   endif()
   
   get_filename_component(HDF_LIB_DIR ${HDF5_LIBRARIES} PATH)
   message(STATUS "HDF_LIB_DIR:" ${HDF_LIB_DIR})
   LINK_DIRECTORIES(${HDF_LIB_DIR})

if(WIN32)
     find_library( HDF5_HL_IMPORT_LIB NAMES hdf5_hl hdf5_hldll )
     LIST(APPEND HDF5_LIBRARIES ${HDF5_HL_IMPORT_LIB} )
     message(STATUS "HDF5_LIBRARIES:" ${HDF5_LIBRARIES})
endif()
IF(APPLE)
  SET(HDF5_LIBSZ "libszip.a" CACHE STRING "only lib name: libszip.a")
  SET(HDF5_LIBZ "libz.a" CACHE STRING "only lib name: libz.a")
    LIST(APPEND HDF5_LIBRARIES ${HDF_LIB_DIR}/${HDF5_LIBSZ} ${HDF_LIB_DIR}/${HDF5_LIBZ})
   message(STATUS "HDF5_LIBRARIES:" ${HDF5_LIBRARIES})
ENDIF(APPLE)

   include_directories(${HDF5_INCLUDE_DIR})
   add_definitions(${HDF5_CPPFLAGS})
   add_definitions(-DWITH_HDF5)
   LIST(APPEND LIBRARY_DEPENDENCIES ${HDF5_LIBRARIES})
else()
   message(STATUS "build without HDF5")
endif()

if((WITH_OPENGM) AND (BLA))
    message(STATUS "build with HDF5 support")
    find_package(HDF5 REQUIRED)
    include_directories(${HDF5_INCLUDE_DIR})
    add_definitions(-DWITH_HDF5)
    SET(HDF5_LIB "${HDF5_LIBRARIES}" CACHE STRING "${HDF5_LIBRARIES}")
    message(STATUS "HDF5_LIBRARIES=" ${HDF5_LIBRARIES})
    message(STATUS "HDF5_LIB=" ${HDF5_LIB})
# does not work for me, HDF5_LIBRARIES is wrong and HDF5_LIBRARY_DIRS is empty
# have to hack it manually
IF(APPLE)
  SET(HDF5_LIBSZ "libszip.a" CACHE STRING "only lib name: libszip.a")
  SET(HDF5_LIBZ "libz.a" CACHE STRING "only lib name: libz.a")
    get_filename_component(HDF_LIB_DIR ${HDF5_LIB} PATH)
    message(STATUS "HDF5_LIB_DIR=" ${HDF_LIB_DIR})
    SET(U_HDF_LIBS ${HDF5_LIB})
    LIST(APPEND U_HDF_LIBS ${HDF_LIB_DIR}/${HDF5_LIBSZ} ${HDF_LIB_DIR}/${HDF5_LIBZ})
    message(STATUS "U_HDF_LIBS=" ${U_HDF_LIBS})
ENDIF(APPLE)
else()
   #message(STATUS "build without OPENGM")
endif()
#--------------------------------------------------------------
# MaxFlow
#--------------------------------------------------------------

if((WITH_OPENGM) AND (MAXFLOW_OPENGM))
	message(STATUS "build with external inference algorithm MaxFlow")
	SET(MAXFLOW_PATCHEDSRCDIR "${OPENGM_INCLUDE_DIR}/../src/external/MaxFlow-v3.02.src-patched/" CACHE STRING "MAXFLOW patched source code directory")
	add_definitions(-DWITH_MAXFLOW)
	include_directories(${MAXFLOW_PATCHEDSRCDIR})
endif((WITH_OPENGM) AND (MAXFLOW_OPENGM))

#---------------compiler setings-------------------------------------------
#add_definitions(-DNO_MEM_MANAGER)

## For Clang-based tools.
##SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF (WITH_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    MESSAGE("-- Found OpenMP.")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    ADD_DEFINITIONS(-DUSE_OPENMP)
  ELSE(OPENMP_FOUND)
    MESSAGE("-- Can't find OpenMP. Continuing without it.")
  ENDIF(OPENMP_FOUND)
ELSE (WITH_OPENMP)
  MESSAGE("-- OpenMP is disabled. Continuing without it.")
ENDIF (WITH_OPENMP)

IF (OPENMP_FOUND)
  IF (NOT MSVC)
    LIST(APPEND LIBRARY_DEPENDENCIES gomp)
  ENDIF (NOT MSVC)
ENDIF (OPENMP_FOUND)

# LIST(APPEND LIBRARY_DEPENDENCIES maxflow3)


#
# MSVC settings
#
IF (MSVC)
  # Disable deprecation warning for standard functions.
  ADD_DEFINITIONS("/wd4996")
  # To get rid of annoying min and max macros if windows.h
  ADD_DEFINITIONS(-DNOMINMAX=1)
  ADD_DEFINITIONS(-D__STDC_UTF_16__)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
  ADD_DEFINITIONS(-D_STL_SECURE_NO_WARNINGS)
  # ADD_DEFINITIONS(-D_HAS_CHAR16_T_LANGUAGE_SUPPORT) # VS does not have it
  # Macro redefinitions
  ADD_DEFINITIONS("/wd4005")
ELSE()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -fpermissive -fPIC")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w -fpermissive -O3 -fPIC")
  ###SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
 
ENDIF (MSVC)

IF(CMAKE_COMPILER_IS_GNUCC)
  LIST(APPEND LIBRARY_DEPENDENCIES rt)
ENDIF()


include_directories("${PROJECT_BINARY_DIR}")
include_directories(${SRC})

# function to find header files
function(find_headers out_header_list_name) # example: find_headers(HH ${SRCLIST})
	LIST(REMOVE_AT ARGV 0) # remove the first argument = out_header_list_name
	SET(HH)
	FOREACH(file ${ARGV})
		#message(${file})
		get_filename_component(ffile ${file} REALPATH) # full absolute path
		get_filename_component(dir ${ffile} PATH) # just path
		#message(${dir})
		set(H)
		set(HXX)
		set(HPP)
		file(GLOB_RECURSE H "${dir}/*.h")
		file(GLOB_RECURSE HXX "${dir}/*.hxx")
		file(GLOB_RECURSE HPP "${dir}/*.hpp")
		LIST(APPEND HH ${H} ${HXX} ${HPP})
	ENDFOREACH()
	LIST(REMOVE_DUPLICATES HH)
	#message(${HH})
	set(${out_header_list_name} ${HH} PARENT_SCOPE)
endfunction(find_headers)

# function to generate nice project files for IDEs
function(source_group_path srcdir)
get_filename_component(SRCDIR ${srcdir} REALPATH) # make sure relative pathes does not get confused
	LIST(REMOVE_AT ARGV 0) # remove the first argument = srcdir
	FOREACH(file ${ARGV})
		get_filename_component(ffile ${file} REALPATH)
		FILE(RELATIVE_PATH pth ${SRCDIR} ${ffile})
		get_filename_component(rdir ${pth} PATH)
		string(REPLACE "/" "\\" rrdir ${rdir})
		get_filename_component(name ${pth} NAME)
        	#MESSAGE(${rrdir} : ${name})
		source_group(${rrdir} FILES ${file})
		#ALLILUJA!
	ENDFOREACH()
endfunction(source_group_path)



message(STATUS "CMAKE_CXX_FLAGS_RELEASE=" ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS "CMAKE_EXE_LINKER_FLAGS=" ${CMAKE_EXE_LINKER_FLAGS})
message(STATUS "LIBRARY_DEPENDENCIES=" ${LIBRARY_DEPENDENCIES})

#___________________________________________
SET (MYLIBS)
#_______________libraries___________________
function(my_add_library TARGET CPP) # example: my_add_library(mylib ${CPP})
	LIST(REMOVE_AT ARGV 0)
	SET(CPP ${ARGV}) # this is the hack around to pass in a list
	message(STATUS "my_add_library(" ${TARGET} ",[" ${CPP} "])")

	find_headers(HH ${CPP})
	source_group_path(${SRC} ${CPP})
	source_group_path(${SRC} ${HH})
	add_library(${TARGET} ${CPP} ${HH})

	target_link_libraries(${TARGET} ${MYLIBS} ${LIBRARY_DEPENDENCIES})
	# copy that into lib folder
	add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${MY_LIB_PATH})
	add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${TARGET}> ${MY_LIB_PATH})
	# add to list
	SET(MYLIBS ${TARGET} ${MYLIBS} PARENT_SCOPE) # insert in front, which is important for linking
endfunction()


#____________exttype lib________________ 
SET(CPP)
list(APPEND CPP 
  ${SRC}/dynamic/block_allocator.cpp
  ${SRC}/dynamic/options.cpp
  ${SRC}/data/dataset.cpp
  ${SRC}/exttype/pvect.cpp
  ${SRC}/geom/vectn.cpp
  ${SRC}/files/xfs.cpp
  ${SRC}/debug/logs.cpp
  ${SRC}/debug/performance.cpp
  ${SRC}/streams/xstringstream.cpp
  ${SRC}/streams/file_stream.cpp
  ${SRC}/streams/binary_stream.cpp
  ${SRC}/streams/text_stream.cpp
)

my_add_library(exttype ${CPP})

#__________maxflow lib__________________ 
SET(CPP)
IF((WITH_OPENGM) AND (MAXFLOW_OPENGM))
	list(APPEND CPP ${MAXFLOW_PATCHEDSRCDIR}/graph.cpp ${MAXFLOW_PATCHEDSRCDIR}/maxflow.cpp)
ELSE((WITH_OPENGM) AND (MAXFLOW_OPENGM))
	list(APPEND CPP ${SRC}/maxflow/graph.cpp ${SRC}/maxflow/maxflow.cpp)
ENDIF((WITH_OPENGM) AND (MAXFLOW_OPENGM))

my_add_library(maxflow_BK ${CPP})

#___________part_opt lib________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/graph/mgraph.cpp
${SRC}/optim/trws/stream_graph.cpp
${SRC}/optim/part_opt/dee.cpp
${SRC}/optim/part_opt/dee2.cpp
${SRC}/optim/part_opt/energy.cpp
${SRC}/optim/part_opt/part_opt_TRWS.cpp
${SRC}/optim/part_opt/part_opt_interface.cpp
${SRC}/optim/part_opt/energy.cpp
${SRC}/optim/part_opt/msg_alg.cpp
${SRC}/optim/part_opt/trws_machine.cpp
${SRC}/optim/part_opt/vectorizers.cpp
)

my_add_library(part_opt ${CPP})

#_____________________________________________
#___________________TESTS_____________________
function(my_add_executable TARGET CPP) # example: my_add_executable(mylib ${CPP})
	LIST(REMOVE_AT ARGV 0)
	SET(CPP ${ARGV}) # this is the hack around to pass in a list
	message(STATUS "my_add_executable(" ${TARGET} ",[" ${CPP} "])")

	#find_headers(HH ${CPP})
	#source_group_path(${SRC} ${CPP})
	#source_group_path(${SRC} ${HH})

	#list(GET ${CPP} 0 CPP0)
	#get_filename_component(TPATH ${CPP0} PATH)
	#source_group_path(${TPATH} FILES ${CPP})

	source_group("\\" FILES ${CPP})
	add_executable(${TARGET} ${CPP} ${HH})

	message(STATUS "libs =" ${MYLIBS})
	target_link_libraries(${TARGET} ${MYLIBS} ${LIBRARY_DEPENDENCIES})
	# copy that into bin folder
	add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${MY_BIN_PATH})
	add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${TARGET}> ${MY_BIN_PATH})
endfunction()

#______________random______________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/part_opt/test_random.cpp
)

my_add_executable(test_random ${CPP})

##-----------------------------------------------------------------------------

SET(CPP
${SRC}/optim/qbranching/qbranching.cpp
${SRC}/optim/qbranching/test.cpp
)
my_add_executable(qbranching ${CPP})
                                 
#______________opengm______________________
if(WITH_OPENGM)
SET(CPP)
list(APPEND CPP
${SRC}/optim/part_opt/test_part_opt_opengm.cpp
)

my_add_executable(test_part_opt_opengm ${CPP})

ENDIF(WITH_OPENGM)

#____________________________________________________
#___________________matlab mexes_____________________
IF(WITH_MATLAB)

set(BUILD_SHARED_LIBS ON)
include_directories(${MATLAB_INCLUDE_DIRS})

LIST(APPEND LIBRARY_DEPENDENCIES ${MATLAB_LIBRARIES}) # ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} ${MATLAB_ENG_LIBRARY})
add_definitions(-DMX_COMPAT_32)
add_definitions(-DMATLAB_MEX_FILE)

if(WIN32)
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /export:mexFunction")
	LIST(APPEND LIBRARY_DEPENDENCIES gdi32 user32)
else()
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pthread -shared -Wl,--no-undefined -fopenmp")
	LIST(APPEND LIBRARY_DEPENDENCIES mx mex mat m eng /usr/lib/x86_64-linux-gnu/libstdc++.so.6)
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath-link,${MATLAB_LIBRARY_PATH}")
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath-link,${MATLAB}/extern/lib/glnxa64")
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath-link,${MATLAB}/bin/glnxa64")
	SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath-link,/usr/lib/x86_64-linux-gnu ${CMAKE_SHARED_LINKER_FLAGS}")
endif()

message(STATUS "MATLAB build settings")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE=" ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS "CMAKE_SHARED_LINKER_FLAGS=" ${CMAKE_SHARED_LINKER_FLAGS})
message(STATUS "LIBRARY_DEPENDENCIES=" ${LIBRARY_DEPENDENCIES})


function(add_mex TARGET CPP) # example: add_mex(my_mex ${CPP})
	LIST(REMOVE_AT ARGV 0)
	SET(CPP ${ARGV}) # this is the hack around to pass in a list
	list(APPEND CPP ${SRC}/mex/mex_io.cpp)
	message(STATUS "add_mex(" ${TARGET} ",[" ${CPP} "])")

	#list(GET ${CPP} 0 CPP0)
	#get_filename_component(TPATH ${CPP0} PATH)
	#source_group_path(${TPATH} FILES ${CPP})
	
	source_group("\\" FILES ${CPP})
        add_library(${TARGET} SHARED ${CPP})

	set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${TARGET} IMPORT_SUFFIX ".${MATLAB_MEXEXT}" IMPORT_PREFIX "" PREFIX "")
	target_link_libraries(${TARGET} part_opt exttype maxflow_BK ${LIBRARY_DEPENDENCIES})
	add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${MEX_INSTALL_PATH})
	add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${TARGET}> ${MEX_INSTALL_PATH}/${TARGET}.${MEXEXT})
endfunction()


#____________________opengm_read_mex____________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/part_opt/opengm_read_mex.cpp
)
add_mex(opengm_read_mex ${CPP})

#____________________opengm_read_mex____________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/part_opt/opengm_read_diff_mex.cpp
)
add_mex(opengm_read_diff_mex ${CPP})


#____________________opengm_write_mex____________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/part_opt/opengm_write_mex.cpp
)
add_mex(opengm_write_mex ${CPP})

#____________________part_opt_TRWS_mex___________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/part_opt/part_opt_TRWS_mex.cpp
${SRC}/optim/part_opt/opengm_read.cpp
)
add_mex(part_opt_TRWS_mex ${CPP})

#____________________part_opt_IK_mex__________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/part_opt/part_opt_IK_mex.cpp
${SRC}/optim/part_opt/emaxflow.cpp
)
add_mex(part_opt_IK_mex ${CPP})

#____________________part_opt_MQPBO_mex__________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/part_opt/part_opt_MQPBO_mex.cpp
${SRC}/qpbo-v1.3/QPBO.cpp
${SRC}/qpbo-v1.3/QPBO_maxflow.cpp
${SRC}/qpbo-v1.3/QPBO_extra.cpp
${SRC}/qpbo-v1.3/QPBO_postprocessing.cpp
)
add_mex(part_opt_MQPBO_mex ${CPP})

#____________________part_opt_DEE_mex__________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/part_opt/part_opt_DEE_mex.cpp
)
add_mex(part_opt_DEE_mex ${CPP})

#____________________part_opt_DEE2_mex_________________
SET(CPP)
list(APPEND CPP 
${SRC}/optim/part_opt/part_opt_DEE2_mex.cpp
)
add_mex(part_opt_DEE2_mex ${CPP})

ENDIF(WITH_MATLAB)