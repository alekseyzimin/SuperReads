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

# Later put back to joinKUnitigs
bin_PROGRAMS = createFastaSuperReadSequences    \
	findMatchesBetweenKUnitigsAndReads \
	createKUnitigMaxOverlaps \
	eliminateBadSuperReadsUsingList \
	getSuperReadInsertCountsFromReadPlacementFile \
	joinKUnitigs_v3 \
	joinKUnitigs_v2

dist_bin_SCRIPTS = getLengthStatisticsForKUnitigsFile.perl createSuperReadsForDirectory.perl

noinst_HEADERS = charb.hpp diskBasedUnitigger.h \
	priorityHeap.h redBlackTreesInsertOnly.h \
	_redBlackTreesInsertOnly.h
# Add source description here such as:
# my_prog_SOURCES = my_prog.cc something_else.cc
createFastaSuperReadSequences_SOURCES = createFastaSuperReadSequences.cc
findMatchesBetweenKUnitigsAndReads_SOURCES =    \
findMatchesBetweenKUnitigsAndReads.cc
createKUnitigMaxOverlaps_SOURCES = createKUnitigMaxOverlaps.cc
# joinKUnitigs_SOURCES = joinKUnitigs.cc redBlackTreesInsertOnly.cc priorityHeap.cc
joinKUnitigs_v2_SOURCES = joinKUnitigs_v2.c _redBlackTreesInsertOnly.c priorityHeap.c
eliminateBadSuperReadsUsingList_SOURCES = eliminateBadSuperReadsUsingList.cc
getSuperReadInsertCountsFromReadPlacementFile_SOURCES = getSuperReadInsertCountsFromReadPlacementFile.cc
joinKUnitigs_v3_SOURCES = joinKUnitigs_v3.cc redBlackTreesInsertOnly.c priorityHeap.c
