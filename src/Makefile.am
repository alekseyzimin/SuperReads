AM_CPPFLAGS = -Wall -Werror $(JELLYFISH_CFLAGS) -I$(top_srcdir)
AM_LDFLAGS = -lpthread $(JELLYFISH_LIBS)

# To add  error_correct_reads
bin_PROGRAMS = create_k_unitigs findMatchesBetweenKUnitigsAndReads	\
        groupConsensusRecordsForRead					\
        extractCoveringKUnitigsFromMinimalCoveringKUnitigSets		\
        createFastaSuperReadSequences					\
        getRecordsForTandomRepeatKUnitigs myUniq			\
        changeReadNumsToReadNamesAtBeginOfLine				\
        convertReadNamesToReadNumbers					\
        getMaxReadNumbersByPrefixForFastaFile				\
        mergeShortMatesIntoSuperReads

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
tryToKillKUnitigsOfLength61.perl

noinst_HEADERS = aligned_simple_array.hpp misc.hpp kmer.hpp	\
diskBasedUnitigger.h fasta_read_parser.hpp			\
create_k_unitigs_cmdline.hpp

create_k_unitigs_SOURCES = create_k_unitigs.cpp	\
create_k_unitigs_cmdline.cpp


#error_correct_reads_SOURCES = error_correct_reads.cc		\
#error_correct_reads_cmdline.cpp kmer.hpp kmer.cc err_log.hpp
#error_correct_reads_SOURCES += error_correct_reads_cmdline.hpp	\
#error_correct_reads.hpp err_log.hpp

findMatchesBetweenKUnitigsAndReads_SOURCES =			\
findMatchesBetweenKUnitigsAndReads.cc fasta_read_parser.cc
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

