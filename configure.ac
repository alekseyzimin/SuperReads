AC_PREREQ([2.63])
AC_INIT([SuperReads], [0.0.1], [gmarcais@umd.edu])
AM_INIT_AUTOMAKE([SuperReads], [0.0.1])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for Jellyfish module
PKG_CHECK_MODULES([JELLYFISH], [jellyfish-1.1])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_REALLOC
AC_CHECK_FUNCS([alarm memset])

# Debugging flag
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [Enable debugging]))
AS_IF([test "x$enable_debug" = "x"], AC_DEFINE([NDEBUG], [1], [Defined if no debugging]))

# Check the version of strerror_r
AC_FUNC_STRERROR_R

AC_OUTPUT([Makefile src/Makefile])

