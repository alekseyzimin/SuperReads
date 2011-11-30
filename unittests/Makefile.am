TESTS = all_tests
check_PROGRAMS = all_tests
all_tests_LDADD = libgtest_main.la

src_dir = $(top_srcdir)/src
all_tests_SOURCES =

all_tests_SOURCES += test_charb.cc test_exp_buffer.cc test_charbuf.cc

AM_CPPFLAGS = -I$(srcdir)/gtest/include -I$(srcdir)/gtest -I$(top_srcdir)
#AM_CXXFLAGS = @PTHREAD_CFLAGS@
#LIBS += @PTHREAD_LIBS@

##############################
# Gtest build.
# Do not edit.
##############################
# Build rules for libraries.
check_LTLIBRARIES = libgtest.la libgtest_main.la

libgtest_la_SOURCES = gtest/src/gtest-all.cc
libgtest_main_la_SOURCES = gtest/src/gtest_main.cc
libgtest_main_la_LIBADD = libgtest.la

# gtest source files that we don't compile directly.  They are
# #included by gtest-all.cc.
GTEST_SRC = gtest/src/gtest-death-test.cc gtest/src/gtest-filepath.cc	\
  gtest/src/gtest-internal-inl.h gtest/src/gtest-port.cc		\
  gtest/src/gtest-printers.cc gtest/src/gtest-test-part.cc		\
  gtest/src/gtest-typed-test.cc gtest/src/gtest.cc

EXTRA_DIST = $(GTEST_SRC)

# Headers, not installed
noinst_HEADERS = gtest/include/gtest/gtest-death-test.h			\
  gtest/include/gtest/gtest-message.h					\
  gtest/include/gtest/gtest-param-test.h				\
  gtest/include/gtest/gtest-printers.h gtest/include/gtest/gtest-spi.h	\
  gtest/include/gtest/gtest-test-part.h					\
  gtest/include/gtest/gtest-typed-test.h gtest/include/gtest/gtest.h	\
  gtest/include/gtest/gtest_pred_impl.h					\
  gtest/include/gtest/gtest_prod.h					\
  gtest/include/gtest/internal/gtest-death-test-internal.h		\
  gtest/include/gtest/internal/gtest-filepath.h				\
  gtest/include/gtest/internal/gtest-internal.h				\
  gtest/include/gtest/internal/gtest-linked_ptr.h			\
  gtest/include/gtest/internal/gtest-param-util-generated.h		\
  gtest/include/gtest/internal/gtest-param-util.h			\
  gtest/include/gtest/internal/gtest-port.h				\
  gtest/include/gtest/internal/gtest-string.h				\
  gtest/include/gtest/internal/gtest-tuple.h				\
  gtest/include/gtest/internal/gtest-type-util.h

