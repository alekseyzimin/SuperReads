AM_CXXFLAGS += $(AM_DEV_CXXFLAGS) -g

# Count lines of code
.PHONY: cloc
cloc:
	cloc --force-lang="Ruby,yaggo" --force-lang="make,am" --force-lang="make,mk" \
	  --exclude-dir="gtest" --ignored=cloc_ignored_src_files \
	  $(srcdir)/src $(srcdir)/src2 $(srcdir)/lib $(srcdir)/include $(srcdir)/unittests \
	  $(srcdir)/Makefile.am $(srcdir)/*.mk


# Make a dependency on yaggo the software
$(YAGGO_SOURCES): $(YAGGO)
