# - Try to find MOAB
# Once done this will define
#
#  MOAB_FOUND - system has MOAB
#  MOAB_INCLUDE_DIRS - the MOAB include directory
#  MOAB_LIBRARIES - Link these to use MOAB
#  MOAB_DEFINITIONS - Compiler switches required for using MOAB
#
#  Copyright (c) 2010 Roman Putanowicz <putanowr@l5.pk.edu.pl>
#
#  Redistribution and use is allowed according to the terms of the New
#  BSD license.
#  For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#


find_path( MOAB_CMAKE_CONFIG NAMES MOABConfig.cmake
  HINTS ${MOAB_ROOT}
  PATHS ENV LD_LIBRARY_PATH
  PATH_SUFFIXES lib Lib
  NO_DEFAULT_PATH
  )

INCLUDE ( ${MOAB_CMAKE_CONFIG}/MOABConfig.cmake )



