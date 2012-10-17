AM_CXXFLAGS += $(AM_DEV_CXXFLAGS) -g

# Count lines of code
.PHONY: cloc
cloc:
	cloc --force-lang="Ruby,yaggo" --force-lang="make,am" --force-lang="make,mk" \
	  --exclude-dir="gtest" --ignored=cloc_ignored_src_files \
	  $(top_srcdir)/src $(top_srcdir)/src2 $(top_srcdir)/lib $(top_srcdir)/include $(top_srcdir)/unittests \
	  $(top_srcdir)/Makefile.am $(top_srcdir)/*.mk
