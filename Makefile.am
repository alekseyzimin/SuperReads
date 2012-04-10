ACLOCAL_AMFLAGS = -I m4
AM_CXXFLAGS = -I$(top_srcdir)/include -I$(top_srcdir) -Wall $(JELLYFISH_CFLAGS) -std=c++0x
AM_LDFLAGS = -lpthread $(JELLYFISH_LIBS)

noinst_HEADERS = $(BUILT_SOURCES)
bin_PROGRAMS =
dist_bin_SCRIPTS =
BUILT_SOURCES =
CLEANFILES = $(BUILT_SOURCES)

EXTRA_DIST =

##############################
# jflib
##############################
JFLIB_I = include/jflib
noinst_HEADERS += $(JFLIB_I)/circular_buffer.hpp			\
                  $(JFLIB_I)/compare_and_swap.hpp			\
                  $(JFLIB_I)/locks_pthread.hpp				\
                  $(JFLIB_I)/multiplexed_io.hpp $(JFLIB_I)/pool.hpp	\
                  $(JFLIB_I)/atomic_field.hpp				\
                  $(JFLIB_I)/multiplexed_parser.hpp

##############################
# utils lib
##############################
UTILS_I = include
noinst_HEADERS += $(UTILS_I)/err.hpp $(UTILS_I)/thread_exec.hpp		\
                  $(UTILS_I)/charb.hpp $(UTILS_I)/charbuf.hpp		\
                  $(UTILS_I)/exp_buffer.hpp $(UTILS_I)/gzip_stream.hpp	\
                  $(UTILS_I)/yaggo.hpp $(UTILS_I)/misc.hpp		\
                  $(UTILS_I)/heap.hpp $(UTILS_I)/reallocators.hpp	\
                  $(UTILS_I)/skip_list_set.hpp				\
                  $(UTILS_I)/skip_list_map.hpp				\
                  $(UTILS_I)/gcc_builtins.hpp $(UTILS_I)/rb_tree.hpp	\
                  $(UTILS_I)/rb_node.hpp				\
                  $(UTILS_I)/binary_search_tree_iterator.hpp		\
                  $(UTILS_I)/skip_list_common.hpp			\
                  $(UTILS_I)/multi_thread_skip_list_set.hpp		\
                  $(UTILS_I)/multi_thread_skip_list_map.hpp		\
                  $(UTILS_I)/tmpstream.hpp $(UTILS_I)/exp_vector.hpp

noinst_LIBRARIES = libutils.a
libutils_a_SOURCES = lib/thread_exec.cc lib/err.cc lib/yaggo.cpp
LDADD = libutils.a

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO " $@;
.yaggo.hpp:
	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -o $@ $<

##############################
# src pipeline
##############################
bin_PROGRAMS += create_k_unitigs error_correct_reads			\
                error_corrected2frg homo_trim reverse_complement	\
                reduce_sr sorted_merge combine_jf_dbs

dist_bin_SCRIPTS += src/getNumBasesPerReadInFastaFile.perl		\
                    src/create_end_pairs.pl		\
                    src/close_gaps.sh src/compute_sr_cov.pl		\
                    src/create_sr_frg.pl src/extractreads.pl		\
                    src/extract_unjoined_pairs.pl src/fasta2frg.pl	\
                    src/filter_library.sh src/fix_unitigs.sh		\
                    src/recompute_astat_superreads.sh src/runSRCA.pl	\
                    src/extractreads_not.pl
noinst_HEADERS += src/aligned_simple_array.hpp			\
                  src/diskBasedUnitigger.h			\
                  src/error_correct_reads.hpp src/err_log.hpp	\
                  src/kmer.hpp src/fibonacci_coding.hpp		\
                  src/sr_names.hpp src/bloom_filter.hpp		\
                  src/MurmurHash3.h src/bloom_hash.hpp		\
                  src/bloom_counter2.hpp src/mer_dna.hpp	\
                  src/read_parser.hpp 
error_corrected2frg_SOURCES = src/error_corrected2frg.cc 
BUILT_SOURCES              += src/error_corrected2frg_cmdline.hpp
homo_trim_SOURCES = src/homo_trim.cc 
reverse_complement_SOURCES = src/reverse_complement.cc
create_k_unitigs_SOURCES = src/create_k_unitigs.cpp
BUILT_SOURCES           += src/create_k_unitigs_cmdline.hpp
error_correct_reads_SOURCES = src/error_correct_reads.cc src/kmer.cc	\
                              src/err_log.cc
BUILT_SOURCES              += src/error_correct_reads_cmdline.hpp
reduce_sr_SOURCES = src/reduce_sr.cc
BUILT_SOURCES    += src/reduce_sr.hpp
sorted_merge_SOURCES = src/sorted_merge.cc
BUILT_SOURCES       += src/sorted_merge.hpp
combine_jf_dbs_SOURCES = src/combine_jf_dbs.cc
BUILT_SOURCES         += src/combine_jf_dbs.hpp
#create_k_unitigs_large_k_SOURCES = src/create_k_unitigs_large_k.cc	\
#                                   src/read_parser.cc src/mer_dna.cc

##############################
# src2 pipeline
##############################
bin_PROGRAMS += createFastaSuperReadSequences				\
                findMatchesBetweenKUnitigsAndReads			\
                createKUnitigMaxOverlaps				\
                eliminateBadSuperReadsUsingList				\
                getSuperReadInsertCountsFromReadPlacementFile		\
                getSuperReadInsertCountsFromReadPlacementFileTwoPasses	\
                joinKUnitigs_v3 outputRecordsNotOnList

dist_bin_SCRIPTS += src2/getLengthStatisticsForKUnitigsFile.perl	\
                    src2/createSuperReadsForDirectory.perl src2/add_missing_mates.pl

createFastaSuperReadSequences_SOURCES = src2/createFastaSuperReadSequences.cc
findMatchesBetweenKUnitigsAndReads_SOURCES =    \
src2/findMatchesBetweenKUnitigsAndReads.cc src/read_parser.cc src/mer_dna.cc
BUILT_SOURCES += src2/findMatchesBetweenKUnitigsAndReads.hpp
createKUnitigMaxOverlaps_SOURCES = src2/createKUnitigMaxOverlaps.cc
eliminateBadSuperReadsUsingList_SOURCES = src2/eliminateBadSuperReadsUsingList.cc
getSuperReadInsertCountsFromReadPlacementFile_SOURCES =	\
src2/getSuperReadInsertCountsFromReadPlacementFile.cc	\
src/fibonacci_coding.cc src/MurmurHash3.cpp
BUILT_SOURCES += src2/getSuperReadInsertCountsFromReadPlacementFile.hpp
getSuperReadInsertCountsFromReadPlacementFileTwoPasses_SOURCES =	\
src2/getSuperReadInsertCountsFromReadPlacementFileTwoPasses.cc		\
src/fibonacci_coding.cc src/MurmurHash3.cpp
BUILT_SOURCES += src2/getSuperReadInsertCountsFromReadPlacementFileTwoPasses.hpp
joinKUnitigs_v3_SOURCES = src2/joinKUnitigs_v3.cc
BUILT_SOURCES          += src2/joinKUnitigs_v3.hpp
outputRecordsNotOnList_SOURCES = src2/outputRecordsNotOnList.cc


##############################
# unittests
##############################
TESTS = all_tests
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
                    unittests/test_mer_dna.cc				\
                    unittests/test_read_parser.cc			\
                    unittests/test_hash.cc				\
                    unittests/test_exp_vector.cc
all_tests_SOURCES += src/fibonacci_coding.cc src/MurmurHash3.cpp	\
                     src/gcc_builtins.cc src/mer_dna.cc			\
                     src/read_parser.cc unittests/misc.cc
all_tests_CXXFLAGS = $(AM_CXXFLAGS)					\
-I$(top_srcdir)/unittests/gtest/include -I$(top_srcdir)unittests/gtest	\
-I$(top_srcdir)/include -I$(top_srcdir)
all_tests_LDADD = libgtest.la
noinst_HEADERS += unittests/test_main.hpp unittests/misc.hpp

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
