# Enforce -Werror for development and remove for distribution
CWERROR = -Werror
dist-hook:
	for i in am in; do \
	  sed -e 's/^CWERROR = -Werror/CWERROR =/' $(distdir)/Makefile.$${i} > $(distdir)/Makefile.$${i}_new; \
	  touch -r $(distdir)/Makefile.$${i} $(distdir)/Makefile.$${i}_new; \
	  mv -f $(distdir)/Makefile.$${i}_new $(distdir)/Makefile.$${i}; \
	done


AM_CPPFLAGS = -Wall $(JELLYFISH_CFLAGS) -I$(top_srcdir)	\
-I$(top_srcdir)/include $(CWERROR)
AM_LDFLAGS = -lpthread $(JELLYFISH_LIBS)

bin_PROGRAMS = # Add executables here

dist_bin_SCRIPTS = # Add scripts here

noinst_HEADERS = # Add headers here

# Add source description here such as:
# my_prog_SOURCES = my_prog.cc something_else.cc
