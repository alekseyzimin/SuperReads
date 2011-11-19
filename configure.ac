AC_PREREQ([2.63])
AC_INIT([SuperReads], [0.0.1], [gmarcais@umd.edu])
AM_INIT_AUTOMAKE([SuperReads], [0.0.1])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_PROG_LIBTOOL

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_CPP
AC_PROG_MKDIR_P
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([rt], [clock_gettime])

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h fcntl.h limits.h netdb.h stddef.h stdint.h stdlib.h string.h strings.h sys/time.h sys/timeb.h unistd.h wchar.h wctype.h])

# Checks for Jellyfish module
PKG_CHECK_MODULES([JELLYFISH], [jellyfish-1.1])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_FUNC_STRERROR_R
AC_FUNC_STRTOD
AC_CHECK_FUNCS([alarm clock_gettime dup2 getcwd getdelim getpagesize gettimeofday isascii memset mkdir munmap regcomp rmdir socket strcasecmp strchr strdup strerror strrchr strstr strtol strtoull])

# Debugging flag
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [Enable debugging]))
AS_IF([test "x$enable_debug" = "x"], AC_DEFINE([NDEBUG], [1], [Defined if no debugging]))

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 unittests/Makefile])
AC_OUTPUT
