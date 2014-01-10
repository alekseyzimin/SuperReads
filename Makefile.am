ACLOCAL_AMFLAGS = -I m4
CXXFLAGS_ALL = -I$(top_srcdir)/include -I$(top_srcdir) -Wall -std=c++0x $(VALGRIND_CFLAGS)
AM_CXXFLAGS = $(CXXFLAGS_ALL) $(JELLYFISH2_0_CFLAGS)
LDADD_ALL = -lrt -lpthread libutils.a
LDADD = $(LDADD_ALL) $(JELLYFISH2_0_LIBS)
AM_LDFLAGS = $(STATIC_FLAGS)
AM_DEV_CXXFLAGS = $(DEV_CXXFLAGS)

noinst_HEADERS = $(YAGGO_SOURCES)
bin_PROGRAMS =
dist_bin_SCRIPTS =
data_DATA =
BUILT_SOURCES = $(YAGGO_SOURCES)
CLEANFILES =
DISTCLEANFILES = $(BUILT_SOURCES)

EXTRA_DIST =

##############################
# jflib
##############################
JFLIB_I = include/jflib
noinst_HEADERS += $(JFLIB_I)/multiplexed_io.hpp $(JFLIB_I)/pool.hpp	\
                  $(JFLIB_I)/multiplexed_parser.hpp

##############################
# utils lib
##############################
UTILS_I = include
noinst_HEADERS += $(UTILS_I)/charb.hpp $(UTILS_I)/charbuf.hpp		\
                  $(UTILS_I)/exp_buffer.hpp				\
                  $(UTILS_I)/gzip_stream.hpp $(UTILS_I)/misc.hpp	\
                  $(UTILS_I)/heap.hpp $(UTILS_I)/reallocators.hpp	\
                  $(UTILS_I)/skip_list_set.hpp				\
                  $(UTILS_I)/skip_list_map.hpp				\
                  $(UTILS_I)/gcc_builtins.hpp $(UTILS_I)/rb_tree.hpp	\
                  $(UTILS_I)/rb_node.hpp				\
                  $(UTILS_I)/binary_search_tree_iterator.hpp		\
                  $(UTILS_I)/skip_list_common.hpp			\
                  $(UTILS_I)/multi_thread_skip_list_set.hpp		\
                  $(UTILS_I)/multi_thread_skip_list_map.hpp		\
                  $(UTILS_I)/tmpstream.hpp $(UTILS_I)/exp_vector.hpp	\
                  $(UTILS_I)/thread_pool.hpp $(UTILS_I)/fork_exec.hpp

noinst_LIBRARIES = libutils.a
libutils_a_SOURCES = lib/fork_exec.cpp
# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO " $@;
.yaggo.hpp:
	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -o $@ $<
YAGGO_SOURCES = # Append all file to be built by yaggo

##############################
# src pipeline
##############################
bin_PROGRAMS += error_corrected2frg reverse_complement reduce_sr	\
                sorted_merge create_k_unitigs_large_k			\
                create_k_unitigs_large_k2 homo_trim bloom_query		\
                filter_overlap_file rename_filter_fastq semaphore

dist_bin_SCRIPTS += src/getNumBasesPerReadInFastaFile.perl		\
                    src/create_end_pairs.pl src/close_gaps.sh		\
                    src/compute_sr_cov.pl src/create_sr_frg.pl		\
                    src2/compute_sr_cov.revisedForGCContig.pl		\
                    src/extractreads.pl src/extract_unjoined_pairs.pl	\
                    src/fasta2frg.pl src/filter_library.sh		\
                    src/fix_unitigs.sh					\
                    src/recompute_astat_superreads.sh src/masurca	\
                    src/runSRCA.pl src/extractreads_not.pl		\
                    src/run_ECR.sh src/filter_alt.pl src/parallel	\
                    src/rename_filter_fastq.pl				\
                    src/filter_redundancy.pl

noinst_HEADERS += src/aligned_simple_array.hpp				\
                  src/diskBasedUnitigger.h src/fibonacci_coding.hpp	\
                  src/sr_names.hpp src/bloom_filter.hpp			\
                  src/MurmurHash3.h src/bloom_hash.hpp			\
                  src/bloom_counter2.hpp src/read_parser.hpp		\
                  src/create_k_unitigs_common.hpp

homo_trim_SOURCES = src/homo_trim.cc
error_corrected2frg_SOURCES = src/error_corrected2frg.cc 
YAGGO_SOURCES              += src/error_corrected2frg_cmdline.hpp
reverse_complement_SOURCES = src/reverse_complement.cc
#create_k_unitigs_SOURCES = src/create_k_unitigs.cpp
#YAGGO_SOURCES           += src/create_k_unitigs_cmdline.hpp
reduce_sr_SOURCES = src/reduce_sr.cc
YAGGO_SOURCES    += src/reduce_sr_cmdline.hpp
sorted_merge_SOURCES = src/sorted_merge.cc
YAGGO_SOURCES       += src/sorted_merge_cmdline.hpp
create_k_unitigs_large_k_SOURCES = src/create_k_unitigs_large_k.cc
YAGGO_SOURCES += src/create_k_unitigs_large_k_cmdline.hpp
create_k_unitigs_large_k2_SOURCES = src/create_k_unitigs_large_k2.cc
YAGGO_SOURCES += src/create_k_unitigs_large_k2_cmdline.hpp
bloom_query_SOURCES = src/bloom_query.cc src/MurmurHash3.cpp
YAGGO_SOURCES += src/bloom_query_cmdline.hpp
filter_overlap_file_SOURCES = src/filter_overlap_file.cc
YAGGO_SOURCES += src/filter_overlap_file_cmdline.hpp
rename_filter_fastq_SOURCES = src/rename_filter_fastq.cc
YAGGO_SOURCES += src/rename_filter_fastq_cmdline.hpp

##############################
# src2 pipeline
##############################
bin_PROGRAMS += createFastaSuperReadSequences				\
                findMatchesBetweenKUnitigsAndReads			\
                createKUnitigMaxOverlaps				\
                eliminateBadSuperReadsUsingList				\
                getSuperReadInsertCountsFromReadPlacementFile		\
                getSuperReadInsertCountsFromReadPlacementFileTwoPasses	\
                joinKUnitigs_v3 outputRecordsNotOnList killBadKUnitigs \
                makeAdjustmentFactorsForNumReadsForAStatBasedOnGC \
		create_sr_frg

dist_bin_SCRIPTS += src2/getLengthStatisticsForKUnitigsFile.perl	\
                    src2/createSuperReadsForDirectory.perl src2/add_missing_mates.pl \
                    src2/outputAlekseysJellyfishReductionFile.perl \
                    src2/translateReduceFile.perl src/mergeSuperReadsUniquely.pl

killBadKUnitigs_SOURCES = src2/killBadKUnitigs.cc
createFastaSuperReadSequences_SOURCES = src2/createFastaSuperReadSequences.cc

findMatchesBetweenKUnitigsAndReads_SOURCES = src/read_parser.cc				\
                                             src2/findMatchesBetweenKUnitigsAndReads.cc	\
                                             src/gcc_builtins.cc
create_sr_frg_SOURCES = src2/create_sr_frg.cc
YAGGO_SOURCES += src2/findMatchesBetweenKUnitigsAndReads_cmdline.hpp

createKUnitigMaxOverlaps_SOURCES = src2/createKUnitigMaxOverlaps.cc
eliminateBadSuperReadsUsingList_SOURCES = src2/eliminateBadSuperReadsUsingList.cc
getSuperReadInsertCountsFromReadPlacementFile_SOURCES =	\
src2/getSuperReadInsertCountsFromReadPlacementFile.cc	\
src/fibonacci_coding.cc src/MurmurHash3.cpp
YAGGO_SOURCES += src2/getSuperReadInsertCountsFromReadPlacementFile_cmdline.hpp
getSuperReadInsertCountsFromReadPlacementFileTwoPasses_SOURCES =	\
src2/getSuperReadInsertCountsFromReadPlacementFileTwoPasses.cc		\
src/fibonacci_coding.cc src/MurmurHash3.cpp
YAGGO_SOURCES += src2/getSuperReadInsertCountsFromReadPlacementFileTwoPasses_cmdline.hpp
joinKUnitigs_v3_SOURCES = src2/joinKUnitigs_v3.cc
YAGGO_SOURCES          += src2/joinKUnitigs_v3_cmdline.hpp
outputRecordsNotOnList_SOURCES = src2/outputRecordsNotOnList.cc
makeAdjustmentFactorsForNumReadsForAStatBasedOnGC_SOURCES = src2/makeAdjustmentFactorsForNumReadsForAStatBasedOnGC.cc

semaphore_SOURCES = src/semaphore.cc
YAGGO_SOURCES += src/semaphore_cmdline.hpp

##############################
# gap closing
##############################
dist_bin_SCRIPTS += src2/closeGaps.perl src2/getEndSequencesOfContigs.perl \
                    src2/create_end_pairs.perl src2/outputJoinedPairs.perl \
                    src2/getSequenceForClosedGaps.perl

##############################
# local gap closing
##############################
dist_bin_SCRIPTS += src2/closeGapsLocally.perl \
	            src2/closeGaps.oneDirectory.perl \
	            src2/getSequenceForLocallyClosedGaps.perl \
                    src2/getMeanAndStdevForGapsByGapNumUsingCeleraTerminatorDirectory.perl \
                    src2/createSuperReadSequenceAndPlacementFileFromCombined.perl \
		    src2/getMeanAndStdevForGapsByGapNumUsingCeleraAsmFile.perl
bin_PROGRAMS += collectReadSequencesForLocalGapClosing \
	        runByDirectory
collectReadSequencesForLocalGapClosing_SOURCES = src2/collectReadSequencesForLocalGapClosing.cc
YAGGO_SOURCES += src2/collectReadSequencesForLocalGapClosing_cmdline.hpp
runByDirectory_SOURCES = src2/runByDirectory.cc
YAGGO_SOURCES += src2/runByDirectory_cmdline.hpp

##############################
# reversed jumping read library elimination
##############################
dist_bin_SCRIPTS += src2/reportReadsToExclude.perl src2/outputMatedReadsAsReverseComplement.perl \
	src2/findReversePointingJumpingReads.perl src2/closeGaps.oneDirectory.fromMinKmerLen.perl
bin_PROGRAMS += putReadsIntoGroupsBasedOnSuperReads
putReadsIntoGroupsBasedOnSuperReads_SOURCES = src2/putReadsIntoGroupsBasedOnSuperReads.cc
YAGGO_SOURCES += src2/putReadsIntoGroupsBasedOnSuperReads_cmdline.hpp

##############################
# reversed jumping read library elimination
##############################
bin_PROGRAMS += extendSuperReadsForUniqueKmerNeighbors extendSuperReadsBasedOnUniqueExtensions
extendSuperReadsForUniqueKmerNeighbors_SOURCES = src2/extendSuperReadsForUniqueKmerNeighbors.cc
extendSuperReadsBasedOnUniqueExtensions_SOURCES = src2/extendSuperReadsBasedOnUniqueExtensions.cc
YAGGO_SOURCES += src2/extendSuperReadsForUniqueKmerNeighbors_cmdline.hpp \
	src2/extendSuperReadsBasedOnUniqueExtensions_cmdline.hpp
##############################
# Celera post-processing
##############################
bin_PROGRAMS += addSurrogatesToFrgCtgFile

dist_bin_SCRIPTS += src2/addSurrogatesToFrgctg.perl src2/getUnitigTypeFromAsmFile.perl

addSurrogatesToFrgCtgFile_SOURCES = src2/addSurrogatesToFrgCtgFile.cc

##############################
# unittests
##############################
TESTS = all_tests
EXTRA_DIST += $(TESTS)
check_PROGRAMS = all_tests

all_tests_SOURCES = unittests/test_main.cc unittests/test_charb.cc	\
                    unittests/test_exp_buffer.cc			\
                    unittests/test_charbuf.cc unittests/test_heap.cc	\
                    unittests/test_sr_names.cc				\
                    unittests/test_bloom_filter.cc			\
                    unittests/test_bloom_counter2.cc			\
                    unittests/test_skip_list_set.cc			\
                    unittests/test_skip_list_map.cc			\
                    unittests/test_rbtree.cc				\
                    unittests/test_multi_thread_skip_list_set.cc	\
                    unittests/test_multi_thread_skip_list_map.cc	\
                    unittests/test_read_parser.cc			\
                    unittests/test_hash.cc				\
                    unittests/test_exp_vector.cc			\
                    unittests/test_thread_pool.cc			\
                    unittests/test_fork_exec.cc
all_tests_SOURCES += src/fibonacci_coding.cc src/MurmurHash3.cpp	\
                     src/gcc_builtins.cc src/read_parser.cc		\
                     unittests/misc.cc
all_tests_CXXFLAGS = $(AM_CXXFLAGS) $(JELLYFISH2_0_CFLAGS)	\
-I$(top_srcdir)/unittests/gtest/include				\
-I$(top_srcdir)unittests/gtest -I$(top_srcdir)/include		\
-I$(top_srcdir)
all_tests_LDFLAGS = $(AM_LDFLAGS) $(JELLYFISH2_0_LIBS)
all_tests_LDADD = libutils.a libgtest.a
YAGGO_SOURCES += unittests/test_main_cmdline.hpp
noinst_HEADERS += unittests/test_main_cmdline.hpp unittests/misc.hpp

include gtest.mk

##############################
# Doxygen
##############################
#include doxygen.am
MOSTLYCLEANFILE = doxygen-doc
EXTRA_DIST += doxygen.cfg
.PHONY: doxygen-doc
doxygen-doc: doxygen.cfg
	rm -rf $@
	doxygen $<

##############################
# Extra configuration for development tree
##############################
-include $(top_srcdir)/development.mk
-include $(top_srcdir)/personal.mk

##############################
# Info
##############################
print-%:
	@echo $($*)

##############################
# data generation
##############################
#dist_data_DATA = data/adapter_0
#CLEANFILES += data/adapter_0
#SUFFIXES = _0 .fa
#.fa_0:
#	$(AM_V_GEN)mkdir -p $(@D); jellyfish count -m 24 -s 5k -C -t 10 -r -o $$(dirname $@)/$$(basename $@ _0) $<

