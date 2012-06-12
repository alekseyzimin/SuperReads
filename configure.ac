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
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [Define to 1 if you have the `clock_gettime' function])])

# Checks for header files.

# Checks for Jellyfish module
PKG_CHECK_MODULES([JELLYFISH], [jellyfish-1.1])

# Check for yaggo
AC_ARG_VAR([YAGGO], [Yaggo switch parser generator])
AS_IF([test "x$YAGGO" = "x"], [AC_PATH_PROG([YAGGO], [yaggo], [false])])

# Checks for typedefs, structures, and compiler characteristics.

# Extra compiler debugging
AC_ARG_VAR([DEV_CXXFLAGS], [Extra compiler flags for development])
if test "x$DEV_CXXFLAGS" = "x"; then
DEV_CXXFLAGS=-Werror
# Check for -Wunused-but-set-variable
my_save_cflags="$CFLAGS"
CFLAGS=-Wunused-but-set-variable
AC_MSG_CHECKING([whether CXX supports -Wunused-but-set-variable])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [DEV_CXXFLAGS="$DEV_CXXFLAGS -Wunused-but-set-variable"],
    [AC_MSG_RESULT([no])])
CFLAGS="$my_save_cflags"
fi
AC_SUBST([DEV_CXXFLAGS])

# Checks for library functions.
AC_FUNC_MMAP
AC_CHECK_FUNCS([mremap])

# Debugging flag
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [Enable debugging]))
AS_IF([test "x$enable_debug" = "x"], AC_DEFINE([NDEBUG], [1], [Defined if no debugging]))

AC_CONFIG_FILES([Makefile doxygen.cfg])

AC_OUTPUT
