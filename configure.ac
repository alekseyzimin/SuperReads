AC_PREREQ([2.63])
AC_INIT([SuperReads], [0.0.1], [alekseyz@ipst.umd.edu], [SuperReads])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_PROG_LIBTOOL

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O # To compile C program with subdir-objects
AC_PROG_INSTALL
AC_PROG_CPP
AC_PROG_MKDIR_P
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.
AC_CHECK_LIB([pthread], [pthread_create])
#AC_CHECK_LIB([rt], [clock_gettime])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [Define to 1 if you have the `clock_gettime' function])])

# Checks for header files.
dnl AC_CHECK_HEADERS([arpa/inet.h fcntl.h limits.h netdb.h stddef.h stdint.h stdlib.h string.h strings.h sys/time.h sys/timeb.h sys/mman.h unistd.h wchar.h wctype.h])

# Checks for Jellyfish module
PKG_CHECK_MODULES([JELLYFISH], [jellyfish-1.1])

# Check for yaggo
AC_ARG_VAR([YAGGO], [Yaggo switch parser generator])
AS_IF([test "x$YAGGO" = "x"], [AC_PATH_PROG([YAGGO], [yaggo], [false])])

# Checks for typedefs, structures, and compiler characteristics.
dnl AC_HEADER_STDBOOL
dnl AC_C_INLINE
dnl AC_TYPE_INT32_T
dnl AC_TYPE_INT64_T
dnl AC_TYPE_OFF_T
dnl AC_TYPE_PID_T
dnl AC_C_RESTRICT
dnl AC_TYPE_SIZE_T
dnl AC_TYPE_SSIZE_T
dnl AC_TYPE_UINT16_T
dnl AC_TYPE_UINT32_T
dnl AC_TYPE_UINT64_T
dnl AC_TYPE_UINT8_T
dnl AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
dnl AC_FUNC_ERROR_AT_LINE
dnl AC_FUNC_FORK
dnl AC_FUNC_MALLOC
AC_FUNC_MMAP
dnl AC_FUNC_REALLOC
dnl AC_FUNC_STRERROR_R
dnl AC_FUNC_STRTOD
dnl AC_CHECK_FUNCS([alarm clock_gettime dup2 getcwd getdelim getpagesize gettimeofday isascii memset mkdir munmap regcomp rmdir socket strcasecmp strchr strdup strerror strrchr strstr strtol strtoull mmap munmap mremap])
AC_CHECK_FUNCS([mremap])

# Debugging flag
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [Enable debugging]))
AS_IF([test "x$enable_debug" = "x"], AC_DEFINE([NDEBUG], [1], [Defined if no debugging]))

AC_CONFIG_FILES([Makefile doxygen.cfg])

AC_OUTPUT
