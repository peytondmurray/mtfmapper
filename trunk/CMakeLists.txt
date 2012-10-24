cmake_minimum_required(VERSION 2.4)

# set the extension used for C++ files here
# this is used to find all sources using a GLOB
set(cpp_ext .cc)

# Force cmake to default to a "Release" build, which implies full
# optimisation
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

project(mtfmapper CXX C)
set(mtfmapper_VERSION_MAJOR 0)
set(mtfmapper_VERSION_MINOR 4)
set(genrectangle_VERSION_MAJOR 1)
set(genrectangle_VERSION_MINOR 0)


# Turn all compiler warnings on when building with GCC
if(CMAKE_COMPILER_IS_GNUCXX)
  message(STATUS "Using a g++ compiler variant --- additional compiler warnings turned on")
  add_definitions(-Wall -W)
endif(CMAKE_COMPILER_IS_GNUCXX)

# Turn all compiler warnings on when building with the Intel compiler
if(CMAKE_CXX_COMPILER MATCHES "icpc")
  message(STATUS "Using an Intel compiler variant --- additional compiler warnings turned on")
  #-wd981 disables a warning about unspecified evaluation order of operands,
  # which appears to be spurious most of the time ...
  #-wd810 disables a warning about a loss of precision. This warning appears
  # to be spurious.
  add_definitions(-Wall -wd981 -wd810)
endif(CMAKE_CXX_COMPILER MATCHES "icpc")

# must check this on new MSVC
# Add stdint.h if we are compiling on MSVC
#if(CMAKE_GENERATOR MATCHES "Visual Studio")
#  INCLUDE_DIRECTORIES(${scl_SOURCE_DIR}/include)
#  add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
#endif(CMAKE_GENERATOR MATCHES "Visual Studio")

if (MSVC)
    add_definitions(/arch:SSE2)
endif(MSVC)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_helpers)
SET(OpenCV_BACKWARD_COMPA true)

# if we are not on windows, try the following for locating OpenCV
# this works well on Ubuntu
# and seems to work on Gentoo
if (NOT MSVC)
    if (EXISTS "/usr/share/OpenCV/OpenCVConfig.cmake")
        include("/usr/share/OpenCV/OpenCVConfig.cmake")
    elseif(EXIST "/usr/share/opencv/OpenCVConfig.cmake")
        include("/usr/share/opencv/OpenCVConfig.cmake")
    endif (EXISTS "/usr/share/OpenCV/OpenCVConfig.cmake")
endif(NOT MSVC)

# external libs
if (NOT OpenCV_VERSION)
    # if neither of the above worked, try the find script
    find_package(OpenCV REQUIRED)
endif (NOT OpenCV_VERSION)
find_package(FFTW REQUIRED)
find_package(TBB REQUIRED)
find_package(TCLAP REQUIRED)
find_package(Qt4)

include_directories("${FFTW_INCLUDE_DIR}")
include_directories("${OpenCV_INCLUDE_DIR}")
include_directories("${TCLAP_INCLUDE_DIR}")
include_directories("${TBB_INCLUDE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_BINARY_DIR}/bin")

#optional
if (BUILD_DOC)
  find_package(LATEX)
  # documenation in LaTeX
  if (LATEX_COMPILER)
    include(cmake_helpers/UseLATEX.cmake)
    add_subdirectory(doc)
  endif (LATEX_COMPILER)
endif (BUILD_DOC)

# man pages for unix platforms
if (UNIX)
  if(NOT MANDIR)
      set(MANDIR "/usr/man")
  endif(NOT MANDIR)
  set(HTML_HELP_DIR ${DOCDIR})
  add_subdirectory(manpages)
endif(UNIX)

if (MSVC)
  if (NOT DOCDIR)
      set(HTML_HELP_DIR "./")
  else(NOT DOCDIR)
      set(HTML_HELP_DIR ${DOCDIR})
  endif(NOT DOCDIR)  
endif (MSVC)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/config.h.in"
  "${PROJECT_BINARY_DIR}/bin/config.h"
  )

add_definitions(-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64)
add_subdirectory(bin)

# copy across the html documentation
if (HTML_HELP_DIR)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION ${HTML_HELP_DIR})
endif(HTML_HELP_DIR)



