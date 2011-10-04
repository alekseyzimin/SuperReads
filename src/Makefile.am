AM_CPPFLAGS = -Wall -Werror $(JELLYFISH_CFLAGS) -I$(top_srcdir)
AM_LDFLAGS = -lpthread $(JELLYFISH_LIBS)

bin_PROGRAMS = create_k_unitigs error_correct_reads		\
        findMatchesBetweenKUnitigsAndReads			\
        groupConsensusRecordsForRead				\
        extractCoveringKUnitigsFromMinimalCoveringKUnitigSets	\
        createFastaSuperReadSequences				\
        getRecordsForTandomRepeatKUnitigs myUniq		\
        changeReadNumsToReadNamesAtBeginOfLine			\
        convertReadNamesToReadNumbers				\
        getMaxReadNumbersByPrefixForFastaFile			\
        mergeShortMatesIntoSuperReads error_corrected2frg filter_illumina_data_quals homo_trim reverse_complement \
	createKUnitigMaxOverlaps getMateStatistics_pt1 \
	joinKUnitigs getReadStartsOffsetsAndOrientationsInSuperReads_1stPass

dist_bin_SCRIPTS = createSuperReadsForDirectory.perl			\
getNumBasesInKUnitigsFile.perl mergeShortMatesIntoSuperReads.perl	\
eliminateInfrequentlyOccurringSuperReadsUsingList.perl			\
getNumReadsFromReadPrefixCountsFile.perl				\
outputSuperReadSeqForJumpLibrary.perl					\
extractMinimalKUnitigSetToCoverReads.perl				\
getReadStartsOffsetsAndOrientationsInSuperReads.1stPass.perl		\
reduceRecordsForKillingKUnitigsToConnectMates.perl			\
findAndKillSuperReadsWithKUnitigOverlapsGeTheMerLength.perl		\
killUnwantedDataFromFirstSuperReadInfo.perl				\
reportFinalReadPlacementsInSuperReads.perl				\
findMatchRecordsForGivenKUnitig.perl					\
listInfrequentlyOccurringSuperReads.perl reportSuperReadGroups.perl	\
findSuperReadGroups.perl						\
makeListOfSuperReadsAndCounts.fromSuperReadGroupsFile.perl		\
tryToKillKUnitigsOfLength61.perl close_gaps.sh compute_sr_cov.pl create_sr_frg.pl extractreads.pl extract_unjoined_pairs.pl fasta2frg.pl filter_library.sh fix_unitigs.sh recompute_astat_superreads.sh runSRCA.pl \
postSuperReadPipelineCommandsForJoiningMates.perl \
createSuperReadGroupsFromMatePairJoinedMates.perl \
createCombinedGoodJoinsFile.perl \
extractMatePairConnectionFromResults.perl \
convertReadMateInfoIntoKUnitigMateInfoForMateJoining.perl \
getStdevOfLibrary.perl getMateStatistics_pt2.perl \
findEquivalentSuperReads.perl \
mergePostMateMergeAndPriorSuperReadGroupsByReadFiles.perl

noinst_HEADERS = aligned_simple_array.hpp misc.hpp kmer.hpp		\
diskBasedUnitigger.h create_k_unitigs_cmdline.hpp			\
error_correct_reads_cmdline.hpp error_correct_reads.hpp err_log.hpp	\
kmer.hpp priorityHeap.h redBlackTreesInsertOnly.h			\
create_k_unitigs_cmdline.hpp yaggo.hpp

error_corrected2frg_SOURCES = error_corrected2frg.c 
filter_illumina_data_quals_SOURCES = filter_illumina_data_quals.c 
homo_trim_SOURCES = homo_trim.c 
reverse_complement_SOURCES = reverse_complement.c

create_k_unitigs_SOURCES = create_k_unitigs.cpp yaggo.cpp

error_correct_reads_SOURCES = error_correct_reads.cc kmer.cc	\
err_log.cc yaggo.cpp

findMatchesBetweenKUnitigsAndReads_SOURCES =	\
findMatchesBetweenKUnitigsAndReads.cc
groupConsensusRecordsForRead_SOURCES = groupConsensusRecordsForRead.c
extractCoveringKUnitigsFromMinimalCoveringKUnitigSets_SOURCES =	\
extractCoveringKUnitigsFromMinimalCoveringKUnitigSets.c
createFastaSuperReadSequences_SOURCES = createFastaSuperReadSequences.c
getRecordsForTandomRepeatKUnitigs_SOURCES = getRecordsForTandomRepeatKUnitigs.c
myUniq_SOURCES = myUniq.c
changeReadNumsToReadNamesAtBeginOfLine_SOURCES = changeReadNumsToReadNamesAtBeginOfLine.c
mergeShortMatesIntoSuperReads_SOURCES = mergeShortMatesIntoSuperReads.cc
convertReadNamesToReadNumbers_SOURCES = convertReadNamesToReadNumbers.c
getMaxReadNumbersByPrefixForFastaFile_SOURCES = getMaxReadNumbersByPrefixForFastaFile.c
getMateStatistics_pt1_SOURCES = getMateStatistics_pt1.cc
createKUnitigMaxOverlaps_SOURCES = createKUnitigMaxOverlaps.cc
joinKUnitigs_SOURCES = joinKUnitigs.c redBlackTreesInsertOnly.c priorityHeap.c
getReadStartsOffsetsAndOrientationsInSuperReads_1stPass_SOURCES = getReadStartsOffsetsAndOrientationsInSuperReads_1stPass.cc
