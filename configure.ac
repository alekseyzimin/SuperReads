AC_PREREQ([2.63])
AC_INIT([SuperReads], [0.3.0], [alekseyz@ipst.umd.edu], [SuperReads])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_LANG(C++)
LT_INIT

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
PKG_CHECK_MODULES([JELLYFISH2_0], [jellyfish-2.0])

# Use valgrind to check memory allocation with mmap
AC_ARG_ENABLE([valgrind],
              [AS_HELP_STRING([--enable-valgrind], [Instrument mmap memory allocation with valgrind])])
AS_IF([test "x$enable_valgrind" = xyes],
      [AC_DEFINE([HAVE_VALGRIND], [1], [Define is using Valgrind])]
      [PKG_CHECK_MODULES([VALGRIND], [valgrind >= 1.8.0])])

# Check for yaggo
AC_ARG_VAR([YAGGO], [Yaggo switch parser generator])
AS_IF([test "x$YAGGO" = "x"], [AC_PATH_PROG([YAGGO], [yaggo], [false])])

# --enable-all-static
AC_ARG_ENABLE([all-static],
              [AC_HELP_STRING([--enable-all-static], [create statically linked executable])])
STATIC_FLAGS=
AS_IF([test "x$enable_all_static" = xyes],
      [AC_SUBST([STATIC_FLAGS], [-all-static])])

# Extra compiler debugging
AC_ARG_VAR([DEV_CXXFLAGS], [Extra compiler flags for development])
if test "x$DEV_CXXFLAGS" = "x"; then
DEV_CXXFLAGS="-Werror -D_FORTIFY_SOURCE=2"
# Check for -Wunused-but-set-variable
my_save_cxxflags="$CXXFLAGS"
CXXFLAGS=-Wunused-but-set-variable
AC_MSG_CHECKING([whether CXX supports -Wunused-but-set-variable])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [DEV_CXXFLAGS="$DEV_CXXFLAGS -Wunused-but-set-variable"],
    [AC_MSG_RESULT([no])])
CXXFLAGS="$my_save_cxxflags"
fi
AC_SUBST([DEV_CXXFLAGS])

# Checks for library functions.
AC_FUNC_MMAP
AC_CHECK_FUNCS([mremap])

# Check for builtin functions
AC_MSG_CHECKING([for __builtin_prefetch])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[int x = 0;]],
                                [[__builtin_prefetch(&x);]])],
               [AC_MSG_RESULT([yes])]
               [AC_DEFINE([HAVE_BUILTIN_PREFETCH], [1], [Defined if __builtin_prefetch is supported])],
               [AC_MSG_RESULT([no])])

# Debugging flag
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [Enable debugging]))
AS_IF([test "x$enable_debug" = "x"], AC_DEFINE([NDEBUG], [1], [Defined if no debugging]))

AC_CONFIG_FILES([Makefile doxygen.cfg])

AC_OUTPUT
