#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([geo-demo], [1.0], [])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_SRCDIR([tests/test-geo.cc])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])

# Checks for programs.
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AC_PROG_CXX
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
m4_ifdef([LT_INIT], [LT_INIT])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
PKG_CHECK_MODULES([GTEST], [gtest])
PKG_CHECK_MODULES([GTEST_MAIN], [gtest_main])
PKG_CHECK_MODULES([GEOS], [geos])
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile
                 tests/Makefile])
AC_OUTPUT
