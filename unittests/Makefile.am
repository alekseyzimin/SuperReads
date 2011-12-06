TESTS = all_tests
check_PROGRAMS = all_tests
all_tests_LDADD = libgtest_main.la

src_dir = $(top_srcdir)/src
all_tests_SOURCES =

all_tests_SOURCES += test_charb.cc test_exp_buffer.cc test_charbuf.cc

AM_CPPFLAGS = -I$(srcdir)/gtest/include -I$(srcdir)/gtest -I$(top_srcdir)
#AM_CXXFLAGS = @PTHREAD_CFLAGS@
#LIBS += @PTHREAD_LIBS@

