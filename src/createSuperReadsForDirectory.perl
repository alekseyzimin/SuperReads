#!/usr/bin/perl
#
# This exec takes a (set of) input read files (in fasta format) and outputs
# the set of super-reads for these reads (in fasta format).
#
# The args are the input fasta files as well as (optionally) the directory
# where you want the work to occur. If the directory is not specified, the
# work is done in the current directory.
# If the work directory doesn't exist then it is created.
# 
# The flags are as follows:
# -l merLen : the length of the k-mer to use for the calculations (31)
# -s tableSize : the size of the table when running jellyfish (2,000,000,000)
# -t numProcessors : the number of processors to run jellyfish and create_k_unitigs (16)
# -M minCoverageToContinueAUnitig : same as the parameter into create_k_unitigs (3)
# -m minCoverageToCreateAFork : same as the parameter into create_k_unitigs (2)
# -kunitigsfile filename : a user-given k-unitigs file; otherwise we calculate
# -mkudisr numBaseDiffs : max base diffs between overlapping k-unitigs in super-reads (0)
# -join-mates : join mate pairs
# -join-shooting : join mate pairs using the shooting method
# -force-join : force the mates from all libraries to be joined
#                                (used when -join-shooting is set)
# -default-mean # : the value used for the mean for -force-join (default 300)
# -default-stdev # : the value used for the stdev for -force-join
#                                 (default: 10% of the default mean)
# -minreadsinsuperread minReads : super-reads containing fewer than numReads
#                                reads will be eliminated (2)
# -noclean : don't clean up the files afterwards
# -debug : output the longer form of output when comparing reads to k-unitigs
# -mikedebug : don't kill off intermediate results
# -jumplibraryreads : we are generating for jump-library reads; a k-unitigs
#                                 file must be specified
# -elim-dupls : eliminate duplicate reads ahead of time (not implemented yet).
# -h : help 
use File::Basename;
use Cwd;
$exeDir = dirname ($0);
$pwd = cwd;
if ($exeDir !~ /^\//) {
    $exeDir = "$pwd/$exeDir"; }

&processArgs;

$maxHashFillFactor = .8;

$successFile = "$workingDirectory/superReads.success";
unlink ($successFile) if (-e $successFile);

$minKUnitigLen = $merLen+1;
if (! -d $workingDirectory) {
    $cmd = "mkdir $workingDirectory";
    print "$cmd\n"; system ($cmd); }
$jellyfishDataPrefix = "$workingDirectory/organismMerCounts";
$jellyfishHashFile = $jellyfishDataPrefix . "_0";
$kUnitigFastaSequencePrefix = "$workingDirectory/guillaumeKUnitigsAtLeast32bases";
$totalKUnitigFastaSequence = "${kUnitigFastaSequencePrefix}_all.fasta";
$totalKUnitigFastaSequenceComplete = $totalKUnitigFastaSequence;
if ($totalKUnitigFastaSequenceComplete !~ /^\//) {
    $totalKUnitigFastaSequenceComplete = "$pwd/$totalKUnitigFastaSequenceComplete"; }
$jellyfishKUnitigDataPrefix = "$workingDirectory/organismMerCountsForKUnitigs";
$jellyfishKUnitigHashFile = $jellyfishKUnitigDataPrefix . "_0";
$kUnitigLengthsFile = "$workingDirectory/kUnitigLengths.txt";
$numKUnitigsFile = "$workingDirectory/numKUnitigs.txt";
$totReadFile = "$workingDirectory/inputReads.fasta";
$myProgOutput1 = "$workingDirectory/testOutput.txt";
$myProgOutput0_1 = "$workingDirectory/numReadsPerPrefix.txt";
$numReadsFile = "$workingDirectory/numReads.txt";

$myProgOutput1prefix = "$workingDirectory/newTestOutput";
$myProgOutput1_1 = "$workingDirectory/testOutput.nucmerLinesOnly.txt";
$myProgOutput1_1prefix = "$workingDirectory/newTestOutput.nucmerLinesOnly";
$myProgOutput2 = "$workingDirectory/arrangedCoordsResultsByRead.txt";
$myProgOutput3 = "$workingDirectory/arrangedCoordsResultsByRead.reducedKUnitigMatches.txt";
$myProgOutput4 = "$workingDirectory/extractedCoveringKUnitigSummariesForReads.txt";
$myProgOutput5 = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.txt";
$myProgOutput6good = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.overlyLongOverlapsEliminated.txt";
$myProgOutput6bad = "$workingDirectory/superReadsWKUnitigOvlGeKmerLen.numReads.superReadName.txt";
$myProgOutput7 = "$workingDirectory/superReadGroups.onePerLine.withReadInfoIncluded.txt";
$myProgOutput8 = "$workingDirectory/superReadSequences.passList.txt";
$sequenceCreationErrorFile1 = "$workingDirectory/createFastaSuperReadSequences.errors.part1.txt";
$sequenceCreationErrorFile2 = "$workingDirectory/createFastaSuperReadSequences.errors.part2.txt";
$sequenceCreationErrorFileCombined = "$workingDirectory/createFastaSuperReadSequences.errors.txt";
# The following line is an output file not specified on any command line
$errorOutput9 = "$workingDirectory/chimeric_read.txt";
$repetitiveKUnitigInfoFile = "$workingDirectory/multiCopyKUnitigs.kUnitig.rptLength.txt";
$myProgOutput10 = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.chimericReadsAndBadSuperReadsKilled.txt";
$myProgOutput14 = "$workingDirectory/superReadCounts.count.superRead.txt";
$myProgOutput15 = "$workingDirectory/infrequentlyOccurringSuperReads.txt";
$myProgOutput16 = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.preMateMerge.txt";
$myProgOutput18 = "$workingDirectory/readPlacementsInSuperReads.preMateMerge.read.superRead.offset.ori.txt";
$myProgOutput22 = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.postMateMerge.txt";
$myProgOutput23 = "$workingDirectory/readPlacementsInSuperReads.postMateMerge.read.superRead.offset.ori.usingReadNumbers.txt";
$myProgOutput23complete = $myProgOutput23;
if ($myProgOutput23complete !~ /^\//) {
    $myProgOutput23complete = "$pwd/$myProgOutput23complete"; }
$myProgOutput24 = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.overlapJoinedMates.txt";
$myProgOutput25 = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.postOverlapJoinedMates.txt";
$myProgOutput26 = "$workingDirectory/findEquivalentSuperReads.equivOutput.txt";
$myProgOutput27 = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.semiFinal.txt";
$myProgOutput28 = "$workingDirectory/superReadGroups.onePerLine.withReadInfoIncluded.semiFinal.txt";
$myProgOutput27 = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.final.txt"; # Re-define here
$myProgOutput28 = "$workingDirectory/superReadGroups.onePerLine.withReadInfoIncluded.final.txt"; # Re-define here
$finalSuperReadSequenceFile = "$workingDirectory/superReadSequences.fasta";
$finalReadPlacementFileUsingReadNumbers = "$workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.usingReadNumbers.txt";
$finalReadPlacementFile = "$workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt";

if ($jumpLibraryReads) {
    $myProgOutput1_1 = "$workingDirectory/testOutput.nucmerLinesOnly.jumpLibrary.txt";
    $myProgOutput1_1prefix = "$workingDirectory/newTestOutput.nucmerLinesOnly.jumpLibrary";
    $myProgOutput2 = "$workingDirectory/arrangedCoordsResultsByRead.jumpLibrary.txt";
    $myProgOutput3 = "$workingDirectory/arrangedCoordsResultsByRead.reducedKUnitigMatches.jumpLibrary.txt";
    $myProgOutput4 = "$workingDirectory/extractedCoveringKUnitigSummariesForReads.jumpLibrary.txt";
    $myProgOutput5 = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.jumpLibrary.txt";
    $myProgOutput6good = "$workingDirectory/superReadGroupsForEachReadWhichHasAGroup.overlyLongOverlapsEliminated.jumpLibrary.txt";
    $myProgOutput6bad = "$workingDirectory/superReadsWKUnitigOvlGeKmerLen.numReads.superReadName.jumpLibrary.txt";
    $myProgOutput7 = "$workingDirectory/superReadGroups.onePerLine.withReadInfoIncluded.jumpLibrary.txt";
    $myProgOutput8 = "$workingDirectory/superReadSequences.jumpLibrary.bothSidesExtended.fasta";
    $myProgOutput12usingReadNumbers = "$workingDirectory/readPlacementsInSuperReads.forJumpLibraryWBothSidesExtended.read.superRead.offset.ori.usingReadNumbers.txt";
    $myProgOutput12 = "$workingDirectory/readPlacementsInSuperReads.forJumpLibraryWBothSidesExtended.read.superRead.offset.ori.txt";
    $finalSuperReadSequenceFile = "$workingDirectory/superReadSequences.jumpLibrary.fasta";
}

# goto startHere;

$cmd = "$exeDir/getMaxReadNumbersByPrefixForFastaFile @fastaFiles > $myProgOutput0_1";
&runCommandAndExitIfBad ($cmd, $myProgOutput0_1, 1);

$cmd = "$exeDir/convertReadNamesToReadNumbers @fastaFiles $myProgOutput0_1 > $totReadFile";
&runCommandAndExitIfBad ($cmd, $totReadFile, 1);

$cmd = "$exeDir/getNumReadsFromReadPrefixCountsFile.perl $myProgOutput0_1 > $numReadsFile";
&runCommandAndExitIfBad ($cmd, $numReadsFile, 1);


if (! $kUnitigsFile) {

  redoJellyfish:
    $cmd = "time jellyfish count -m $merLen -r -o $jellyfishDataPrefix -c 6 -p 126--both-strands -s $tableSize -t $numProcessors @fastaFiles";
    &runCommandAndExitIfBad ($cmd, "", 0);
    
    $tableResizeFactor = &returnTableResizeAmount ($jellyfishDataPrefix, $jellyfishHashFile);
    if ($tableResizeFactor > 1) {
	$tableSize *= 2;
	print "Resizing the table to $tableSize for jellyfish\n";
	goto redoJellyfish; }
    
    &runCommandAndExitIfBad ("", $jellyfishHashFile, 1);
    
    $cmd = "create_k_unitigs -C -m $minCoverageToCreateAFork -M $minCoverageToContinueAUnitig -l $minKUnitigLen -o $kUnitigFastaSequencePrefix -t $numProcessors $jellyfishHashFile";
    &runCommandAndExitIfBad ($cmd, "${kUnitigFastaSequencePrefix}_${i}.fa", 1); 

    for ($i=0; 1; $i++) {
	$fn = "${kUnitigFastaSequencePrefix}_${i}.fa";
	last unless (-e $fn);
	if (! $mikedebug) {
	    $fn2 = "${kUnitigFastaSequencePrefix}_${i}.counts";
	    unlink ($fn2); }
	push (@guillaumeKUnitigFiles, $fn); }

    $cmd = "cat @guillaumeKUnitigFiles > $totalKUnitigFastaSequence";
    print "$cmd\n"; system ($cmd); }
else { # A k-unitigs file was passed on the command line
    if ($kUnitigsFile ne $totalKUnitigFastaSequence) {
	if ($kUnitigsFile !~ /^\//) {
	    $kUnitigsFile = "$pwd/$kUnitigsFile"; }
	$cmd = "ln -s $kUnitigsFile $totalKUnitigFastaSequence";
	print "$cmd\n"; system ($cmd); } }
if (! $mikedebug) { &killFiles ($jellyfishHashFile, @guillaumeKUnitigFiles); }
    
$cmd = "cat $totalKUnitigFastaSequence | $exeDir/getNumBasesInKUnitigsFile.perl > $kUnitigLengthsFile";
&runCommandAndExitIfBad ($cmd, $kUnitigLengthsFile, 1);

$numKUnitigs = &getNumLines ($kUnitigLengthsFile);
$cmd = "echo $numKUnitigs > $numKUnitigsFile";
&runCommandAndExitIfBad ($cmd, $numKUnitigsFile, 1);

$minSizeNeededForTable = &reportMinJellyfishTableSizeForKUnitigs;

redoKUnitigsJellyfish:
$cmd = "jellyfish count -m $merLen -r -o $jellyfishKUnitigDataPrefix -c 6 -p 126 --both-strands -s $minSizeNeededForTable -t $numProcessors $totalKUnitigFastaSequence";
&runCommandAndExitIfBad ($cmd, $jellyfishKUnitigHashFile, 1);
    
$tableResizeFactor = &returnTableResizeAmount ($jellyfishKUnitigDataPrefix, $jellyfishKUnitigHashFile);
if ($tableResizeFactor > 1) {
    $tableSize *= 2;
    print "Resizing the table to $tableSize for the k-unitig jellyfish run\n";
    goto redoKUnitigsJellyfish; }

if ($debug) {
    $cmd = "$exeDir/findMatchesBetweenKUnitigsAndReads -l $jellyfishKUnitigHashFile -t $numProcessors -p $myProgOutput1prefix $totalKUnitigFastaSequence $numKUnitigsFile $totReadFile";
    &runCommandAndExitIfBad ($cmd, "", 0);
    $cmd = "grep -h \"^myNucmerLine\" ${myProgOutput1prefix}_* | $exeDir/myUniq > $myProgOutput1_1";
    &runCommandAndExitIfBad ($cmd, $myProgOutput1_1, 1); }
else {
    $cmd = "$exeDir/findMatchesBetweenKUnitigsAndReads $jellyfishKUnitigHashFile -t $numProcessors -p $myProgOutput1_1prefix $totalKUnitigFastaSequence $numKUnitigsFile $totReadFile";
    &runCommandAndExitIfBad ($cmd, "", 0);
    $cmd = "cat ${myProgOutput1_1prefix}_* | $exeDir/myUniq > $myProgOutput1_1";
    &runCommandAndExitIfBad ($cmd, $myProgOutput1_1, 1); }
$cmd = "\\rm ${myProgOutput1_1prefix}_*"; print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($jellyfishKUnitigHashFile, $totReadFile); }
    
$cmd = "$exeDir/groupConsensusRecordsForRead $myProgOutput1_1 > $myProgOutput2";
&runCommandAndExitIfBad ($cmd, $myProgOutput2, 1);

$cmd = "$exeDir/extractMinimalKUnitigSetToCoverReads.perl $myProgOutput2 > $myProgOutput3";
&runCommandAndExitIfBad ($cmd, $myProgOutput3, 1);
if (! $mikedebug) { &killFiles ($myProgOutput2); }

$cmd = "$exeDir/extractCoveringKUnitigsFromMinimalCoveringKUnitigSets $myProgOutput3 > $myProgOutput4";
&runCommandAndExitIfBad ($cmd, $myProgOutput4, 1);

$cmd = "cat $myProgOutput4 | $exeDir/findSuperReadGroups.perl > $myProgOutput5";
&runCommandAndExitIfBad ($cmd, $myProgOutput5, 1);
if (! $mikedebug) { &killFiles ($myProgOutput4); }

$cmd = "$exeDir/findAndKillSuperReadsWithKUnitigOverlapsGeTheMerLength.perl $myProgOutput5 $merLen $myProgOutput6good $myProgOutput6bad";
&runCommandAndExitIfBad ($cmd, $myProgOutput6good, 1);
if (! $mikedebug) { &killFiles ($myProgOutput5); }

$cmd = "cat $myProgOutput6good | $exeDir/reportSuperReadGroups.perl > $myProgOutput7";
&runCommandAndExitIfBad ($cmd, $myProgOutput7, 1);

if ($jumpLibraryReads) {
    goto jumpLibraryCalculations; }

$cmd = "$exeDir/createFastaSuperReadSequences $workingDirectory $myProgOutput7 -seqdiffmax $seqDiffMax -error-filename $sequenceCreationErrorFile1 -nosequence > $myProgOutput8";
&runCommandAndExitIfBad ($cmd, $myProgOutput8, 0);
if (! $mikedebug) { &killFiles ($myProgOutput7); }

# The following generate the following files in $workingDirectory:
# 1) chimeric_read.txt: A list of reads with k-unitigs which overlap themselves
#      by too much (unless all matches have a small relative offset)
# 2) multiCopyKUnitigs.kUnitig.rptLength.txt: A list of k-unitigs which
#      have multiple copies overlapping a given read; don't use for merging
#      mates into super-reads.
$cmd = "cat $myProgOutput1_1 | $exeDir/getRecordsForTandomRepeatKUnitigs | $exeDir/reduceRecordsForKillingKUnitigsToConnectMates.perl $workingDirectory";
&runCommandAndExitIfBad ($cmd, "", 0);
if (! $mikedebug) { if (! $joinShooting) { &killFiles ($myProgOutput1_1); } }

# =====================================================================
# PUTTING IN NEW STUFF HERE (3/15/11)
$cmd = "cat $myProgOutput6good | $exeDir/killUnwantedDataFromFirstSuperReadInfo.perl $workingDirectory > $myProgOutput10";
&runCommandAndExitIfBad ($cmd, $myProgOutput10, 1);
if (! $mikedebug) { &killFiles ($myProgOutput6good); }

$cmd = "cat $myProgOutput10 | $exeDir/makeListOfSuperReadsAndCounts.fromSuperReadGroupsFile.perl > $myProgOutput14";
&runCommandAndExitIfBad ($cmd, $myProgOutput14, 1);

$maxReadCountToEliminateSuperRead = $minReadsInSuperRead-1;
$cmd = "cat $myProgOutput14 | $exeDir/listInfrequentlyOccurringSuperReads.perl $maxReadCountToEliminateSuperRead > $myProgOutput15";
&runCommandAndExitIfBad ($cmd, $myProgOutput15, 0);
if (! $mikedebug) { &killFiles ($myProgOutput14); }

$cmd = "$exeDir/eliminateInfrequentlyOccurringSuperReadsUsingList.perl $myProgOutput10 $myProgOutput15 > $myProgOutput16";
&runCommandAndExitIfBad ($cmd, $myProgOutput16, 1);
if (! $mikedebug) { &killFiles ($myProgOutput10); }

$cmd = "$exeDir/getReadStartsOffsetsAndOrientationsInSuperReads_1stPass $myProgOutput3 $myProgOutput16 $kUnitigLengthsFile $sequenceCreationErrorFile1 > $myProgOutput18";
&runCommandAndExitIfBad ($cmd, $myProgOutput18, 1);

if ($joinMates) {
    $cmd = "$exeDir/mergeShortMatesIntoSuperReads $myProgOutput16 $numKUnitigsFile $errorOutput9 $repetitiveKUnitigInfoFile $kUnitigLengthsFile $numReadsFile > $myProgOutput22"; }
else {
    $myProgOutput16complete = $myProgOutput16;
    if ($myProgOutput16complete !~ /^\//) {
	$myProgOutput16complete = "$pwd/$myProgOutput16complete"; }
    
    $cmd = "ln -s $myProgOutput16complete $myProgOutput22"; }
&runCommandAndExitIfBad ($cmd, $myProgOutput22, 1);

# $cmd = "$exeDir/reportFinalReadPlacementsInSuperReads.perl $kUnitigLengthsFile $myProgOutput18 $myProgOutput22 > $myProgOutput23";
$cmd = "$exeDir/getReadStartsOffsetsAndOrientationsInSuperReads_1stPass $myProgOutput3 $myProgOutput22 $kUnitigLengthsFile $sequenceCreationErrorFile1 > $myProgOutput23";
&runCommandAndExitIfBad ($cmd, $myProgOutput23, 1);

if ($joinShooting) {
    chdir ($workingDirectory);
    $cmd = "$exeDir/postSuperReadPipelineCommandsForJoiningMates.perl $forceJoin $defaultMean $defaultStdev -l $merLen -kunitig-files-prefix $totalKUnitigFastaSequenceComplete -read-placements-file $myProgOutput23complete";
    if ($mikedebug) {
	$cmd .= " -mikedebug"; }
    &runCommandAndExitIfBad ($cmd, "", 0);
    chdir ($pwd);
    $cmd = "$exeDir/mergePostMateMergeAndPriorSuperReadGroupsByReadFiles.perl $myProgOutput24 $myProgOutput22 > $myProgOutput25";
    &runCommandAndExitIfBad ($cmd, $myProgOutput25, 1);
    if (! $mikedebug) { &killFiles ($myProgOutput22, $myProgOutput24, $myProgOutput1_1); }
}
else {
    $myProgOutput22complete = $myProgOutput22;
    if ($myProgOutput22complete !~ /^\//) {
	$myProgOutput22complete = "$pwd/$myProgOutput22complete"; }
    $cmd = "ln -s $myProgOutput22complete $myProgOutput25";
    &runCommandAndExitIfBad ($cmd, $myProgOutput22complete, 1);
}
    
$cmd = "$exeDir/findEquivalentSuperReads.perl $myProgOutput25 -equiv-file $myProgOutput26 -out-file $myProgOutput27";
&runCommandAndExitIfBad ($cmd, $myProgOutput27, 1);
if (! $mikedebug) { &killFiles ($myProgOutput26, $myProgOutput16, $myProgOutput22, $myProgOutput24, $myProgOutput25); }

$cmd = "cat $myProgOutput27 | $exeDir/reportSuperReadGroups.perl > $myProgOutput28";
&runCommandAndExitIfBad ($cmd, $myProgOutput28, 1);

$cmd = "$exeDir/createFastaSuperReadSequences $workingDirectory $myProgOutput28 -seqdiffmax $seqDiffMax -error-filename $sequenceCreationErrorFile2 > $finalSuperReadSequenceFile";
&runCommandAndExitIfBad ($cmd, $finalSuperReadSequenceFile, 1);

$cmd = "cat $sequenceCreationErrorFile1 $sequenceCreationErrorFile2 > $sequenceCreationErrorFileCombined";
&runCommandAndExitIfBad ($cmd, $sequenceCreationErrorFileCombined, 0);
if (! $mikedebug) { &killFiles ($sequenceCreationErrorFile1, $sequenceCreationErrorFile2); }

$cmd = "$exeDir/getReadStartsOffsetsAndOrientationsInSuperReads_1stPass $myProgOutput3 $myProgOutput27 $kUnitigLengthsFile $sequenceCreationErrorFileCombined > $finalReadPlacementFileUsingReadNumbers";
&runCommandAndExitIfBad ($cmd, $finalReadPlacementFileUsingReadNumbers, 1);
if (! $mikedebug) { &killFiles ($myProgOutput3, $myProgOutput18, $myProgOutput23); }

$cmd = "$exeDir/changeReadNumsToReadNamesAtBeginOfLine $myProgOutput0_1 $finalReadPlacementFileUsingReadNumbers > $finalReadPlacementFile";
&runCommandAndExitIfBad ($cmd, $finalReadPlacementFile, 1);
if (! $mikedebug) { &killFiles ($myProgOutput0_1, $myProgOutput8); }

#
#
#
#$cmd = "time cat $myProgOutput16 | $exeDir/reportSuperReadGroups.perl > $myProgOutput92";
#print "$cmd\n"; system ($cmd);
#if (! $mikedebug) {
#    if (-f $myProgOutput91) {
#	&killFiles ($myProgOutput16); }
#}
#
#$cmd = "time $exeDir/createFastaSuperReadSequences $workingDirectory $myProgOutput92 -seqdiffmax $seqDiffMax > $finalSuperReadSequenceFile";
#print "$cmd\n"; system ($cmd);
#
#$cmd = "time $exeDir/reportFinalReadPlacementsInSuperReads.perl $kUnitigLengthsFile $myProgOutput18 $myProgOutput91 > $finalReadPlacementFileUsingReadNumbers";
#print "$cmd\n"; system ($cmd);
#
$cmd = "touch $successFile";
system ($cmd);

exit (0);

jumpLibraryCalculations:
$cmd = "$exeDir/createFastaSuperReadSequences $workingDirectory $myProgOutput7 -seqdiffmax $seqDiffMax  -error-filename $sequenceCreationErrorFile1 > $myProgOutput8";
&runCommandAndExitIfBad ($cmd, $myProgOutput8, 1);
if (! $mikedebug) { &killFiles ($myProgOutput7); }

$cmd = "$exeDir/getReadStartsOffsetsAndOrientationsInSuperReads_1stPass $myProgOutput3 $myProgOutput6good $kUnitigLengthsFile $sequenceCreationErrorFile1 > $myProgOutput12usingReadNumbers";
&runCommandAndExitIfBad ($cmd, $myProgOutput12usingReadNumbers, 1);
if (! $mikedebug) { &killFiles ($myProgOutput3, $myProgOutput6good); }

$cmd = "$exeDir/changeReadNumsToReadNamesAtBeginOfLine $myProgOutput0_1 $myProgOutput12usingReadNumbers > $myProgOutput12";
&runCommandAndExitIfBad ($cmd, $myProgOutput12, 1);

$cmd = "$exeDir/outputSuperReadSeqForJumpLibrary.perl $myProgOutput8 $myProgOutput12 > $finalSuperReadSequenceFile";
&runCommandAndExitIfBad ($cmd, $finalSuperReadSequenceFile, 1);
if (! $mikedebug) { &killFiles ($myProgOutput12usingReadNumbers); }

$cmd = "touch $successFile";
system ($cmd);

sub processArgs
{
    $tableSize = 2000000000;
    $merLen = 31;
    $numProcessors = 16;
    $minCoverageToCreateAFork = 2;
    $minCoverageToContinueAUnitig = 3;
    $minReadsInSuperRead = 2;
    $joinMates = 0;
    $joinShooting = 0;
    $seqDiffMax = 0;
    $elimDupls = 0;
    $clean = 1;
    $help = 0;
    if ($#ARGV < 0) {
	$help = 1; }
    for ($i=0; $i<=$#ARGV; $i++) {
	if ($ARGV[$i] eq "-l") {
	    ++$i;
	    $merLen = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-s") {
	    ++$i;
	    $tableSize = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-t") {
	    ++$i;
	    $numProcessors = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-M") {
	    ++$i;
	    $minCoverageToContinueAUnitig = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-kunitigsfile") {
	    ++$i;
	    $kUnitigsFile = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-m") {
	    ++$i;
	    $minCoverageToCreateAFork = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-mkudisr") {
	    ++$i;
	    $seqDiffMax = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-minreadsinsuperread") {
	    ++$i;
	    $minReadsInSuperRead = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-default-mean") {
	    ++$i;
	    $defaultMean = "-default-mean " . $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-default-stdev") {
	    ++$i;
	    $defaultStdev = "-default-stdev " . $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-noclean") {
	    $clean = 0;
	    next; }
	elsif ($ARGV[$i] eq "-elim-dupls") {
	    $elimDupls = 1;
	    next; }
	elsif ($ARGV[$i] eq "-jumplibraryreads") {
	    $jumpLibraryReads = 1;
	    next; }
	elsif ($ARGV[$i] eq "-join-mates") {
	    $joinMates = 1;
	    next; }
	elsif ($ARGV[$i] eq "-join-shooting") {
	    $joinShooting = 1;
	    next; }
	elsif ($ARGV[$i] eq "-force-join") {
	    $forceJoin = "-force-join";
	    next; }
	elsif ($ARGV[$i] eq "-h") {
	    $help = 1;
	    next; }
	elsif ($ARGV[$i] eq "-debug") {
	    $debug = 1;
	    next; }
	elsif ($ARGV[$i] eq "-mikedebug") {
	    $mikedebug = 1;
	    next; }
	elsif (-f $ARGV[$i]) {
	    push (@fastaFiles, $ARGV[$i]);
	    next; }
	else {
	    if (! $workingDirectory) {
		$workingDirectory = $ARGV[$i]; }
	    else {
		print "Working directory was specified as \"$workingDirectory\", now specified as ",$ARGV[$i],".\nThis message could also occur if an input read file doesn't exist.\n";
		$help = 1; } } }
    if ($#fastaFiles < 0) {
	$help = 1; }
    if ($jumpLibraryReads && (! $kUnitigsFile)) {
	print "A k-unitigs file must be supplied with using jump library reads\n";
	$help = 1; }
    if ($help) {
	&giveUsageAndExit; }
}

sub giveUsageAndExit
{
    open (FILE, $0);
    $line = <FILE>;
    while ($line = <FILE>) {
	chomp ($line);
	last unless ($line =~ /^\#/);
	substr ($line, 0, 2) = "";
	print "$line\n"; }
    exit (0);
}
	
sub returnTableResizeAmount
{
    my ($dataPrefix, $hashFile) = @_;

    for ($i=1; 1; $i++) {
	$fn = $dataPrefix . "_$i";
	last unless (-e $fn);
	unlink ($fn); }
    unlink ($hashFile) if ($i > 1);
    return ($i);
}

sub getNumLines
{
    my ($fn) = @_;
    my ($cmd, $line, @flds);

    $cmd = "wc -l $fn |";
    open (CRAZYFILE, $cmd);
    $line = <CRAZYFILE>;
    close (CRAZYFILE);
    @flds = split (" ", $line);
    return ($flds[0]);
}

sub reportMinJellyfishTableSizeForKUnitigs
{
    my ($line, $numKMersInKUnitigs, $numKUnitigs, $minSizeNeededForTable, @flds);

    open (FILE, $kUnitigLengthsFile);
    while ($line = <FILE>) {
	chomp ($line);
	@flds = split (" ", $line);
	$numKMersInKUnitigs += $flds[1]; }
    close (FILE);
    
    open (FILE, $numKUnitigsFile);
    $numKUnitigs = <FILE>;
    chomp ($numKUnitigs);
    close (FILE);
    $numKMersInKUnitigs -= ($numKUnitigs * ($merLen-1));
    $minSizeNeededForTable = int ($numKMersInKUnitigs/$maxHashFillFactor + .5);
    return ($minSizeNeededForTable);
}

sub killFiles
{
    my (@filesToKill) = @_;
    my ($file);

    for (@filesToKill) {
	$file = $_;
	if (-e $file) {
	    unlink ($file); } }
}

# If localCmd is set, it captures the return code and makes
# sure it ran properly. Otherwise it exits.
# If one sets the fileName then we assume it must exist
# With minSize one can set the minimum output file size
sub runCommandAndExitIfBad
{
    my ($localCmd, $fileName, $minSize) = @_;
    my ($retCode, $exitValue, $sz);

    if ($localCmd =~ /\S/) {
        print "$localCmd\n";
        system ("time $localCmd");
        $retCode = $?;
        if ($retCode == -1) {
            print STDERR "failed to execute: $!\n";
            exit ($retCode); }
        elsif ($retCode & 127) {
            printf STDERR "child died with signal %d, %s coredump\n",
            ($retCode & 127), ($retCode & 128) ? 'with' : 'without';
            exit ($retCode); }
        else {
            $exitValue = $retCode >> 8;
            if ($exitValue == 255) { $exitValue = -1; }
            if ($exitValue != 0) {
                printf STDERR "child exited with value %d\n", $exitValue;
                print STDERR "Command \"$localCmd\" failed. Bye!\n";
                exit ($exitValue); }
        }
    }
    return unless ($fileName =~ /\S/);
    if (! -e $fileName) {
        print STDERR "Output file \"$fileName\" doesn't exist. Bye!\n";
        exit (1); }
    $sz = -s $fileName;
    if ($sz < $minSize) {
        print STDERR "Output file \"$fileName\" is of size $sz, must be at least of size $minSize. Bye!\n";
        exit (1); }
}

