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

&processArgs;

$maxHashFillFactor = .8;

$minKUnitigLen = $merLen+1;
if (! -d $workingDirectory) {
    $cmd = "mkdir $workingDirectory";
    print "$cmd\n"; system ($cmd); }
$jellyfishDataPrefix = "$workingDirectory/organismMerCounts";
$jellyfishHashFile = $jellyfishDataPrefix . "_0";
$kUnitigFastaSequencePrefix = "$workingDirectory/guillaumeKUnitigsAtLeast32bases";
$totalKUnitigFastaSequence = "${kUnitigFastaSequencePrefix}_all.fasta";
$kUnitigFastaSequencePrefixComplete = $kUnitigFastaSequencePrefix;
if ($kUnitigFastaSequencePrefixComplete !~ /^\//) {
    $kUnitigFastaSequencePrefixComplete = "$pwd/$kUnitigFastaSequencePrefixComplete"; }
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

$cmd = "time $exeDir/getMaxReadNumbersByPrefixForFastaFile @fastaFiles > $myProgOutput0_1";
print "$cmd\n"; system ($cmd);

$cmd = "time $exeDir/convertReadNamesToReadNumbers @fastaFiles $myProgOutput0_1 > $totReadFile";
print "$cmd\n"; system ($cmd);

$cmd = "time $exeDir/getNumReadsFromReadPrefixCountsFile.perl $myProgOutput0_1 > $numReadsFile";
print "$cmd\n"; system ($cmd);

if (! $kUnitigsFile) {

  redoJellyfish:
    $cmd = "time jellyfish count -m $merLen -r -o $jellyfishDataPrefix -c 6 -p 126--both-strands -s $tableSize -t $numProcessors @fastaFiles";
    print "$cmd\n"; system ($cmd);
    
    $tableResizeFactor = &returnTableResizeAmount ($jellyfishDataPrefix, $jellyfishHashFile);
    if ($tableResizeFactor > 1) {
	$tableSize *= 2;
	print "Resizing the table to $tableSize for jellyfish\n";
	goto redoJellyfish; }
    
    print "Making sure the hash table is in memory...\n";
    $cmd = "cat $jellyfishHashFile > /dev/null";
    print "$cmd\n"; system ($cmd);
    
    $cmd = "time create_k_unitigs -C -m $minCoverageToCreateAFork -M $minCoverageToContinueAUnitig -l $minKUnitigLen -o $kUnitigFastaSequencePrefix -t $numProcessors $jellyfishHashFile";
    print "$cmd\n"; system ($cmd);

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
print "$cmd\n"; system ($cmd);

$numKUnitigs = &getNumLines ($kUnitigLengthsFile);
$cmd = "echo $numKUnitigs > $numKUnitigsFile";
print "$cmd\n"; system ($cmd);

$minSizeNeededForTable = &reportMinJellyfishTableSizeForKUnitigs;

redoKUnitigsJellyfish:
$cmd = "time jellyfish count -m $merLen -r -o $jellyfishKUnitigDataPrefix -c 6 -p 126 --both-strands -s $minSizeNeededForTable -t $numProcessors $totalKUnitigFastaSequence";
print "$cmd\n"; system ($cmd);
    
$tableResizeFactor = &returnTableResizeAmount ($jellyfishKUnitigDataPrefix, $jellyfishKUnitigHashFile);
if ($tableResizeFactor > 1) {
    $tableSize *= 2;
    print "Resizing the table to $tableSize for the k-unitig jellyfish run\n";
    goto redoKUnitigsJellyfish; }

print "Making sure the k-unitig hash table is in memory...\n";
$cmd = "cat $jellyfishKUnitigHashFile > /dev/null";
print "$cmd\n"; system ($cmd);

if ($debug) {
    $cmd = "time $exeDir/findMatchesBetweenKUnitigsAndReads -l $jellyfishKUnitigHashFile -t $numProcessors -p $myProgOutput1prefix $totalKUnitigFastaSequence $numKUnitigsFile $totReadFile";
    print "$cmd\n"; system ($cmd);
    $cmd = "grep -h \"^myNucmerLine\" ${myProgOutput1prefix}_* | $exeDir/myUniq > $myProgOutput1_1";
    print "$cmd\n"; system ($cmd); }
else {
    $cmd = "time $exeDir/findMatchesBetweenKUnitigsAndReads $jellyfishKUnitigHashFile -t $numProcessors -p $myProgOutput1_1prefix $totalKUnitigFastaSequence $numKUnitigsFile $totReadFile";
    print "$cmd\n"; system ($cmd);
    $cmd = "time cat ${myProgOutput1_1prefix}_* | $exeDir/myUniq > $myProgOutput1_1";
    print "$cmd\n"; system ($cmd); }
$cmd = "\\rm ${myProgOutput1_1prefix}_*"; print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($jellyfishKUnitigHashFile, $totReadFile); }
    
$cmd = "time $exeDir/groupConsensusRecordsForRead $myProgOutput1_1 > $myProgOutput2";
print "$cmd\n"; system ($cmd);

$cmd = "time $exeDir/extractMinimalKUnitigSetToCoverReads.perl $myProgOutput2 > $myProgOutput3";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput2); }

$cmd = "time $exeDir/extractCoveringKUnitigsFromMinimalCoveringKUnitigSets $myProgOutput3 > $myProgOutput4";
print "$cmd\n"; system ($cmd);

$cmd = "time cat $myProgOutput4 | $exeDir/findSuperReadGroups.perl > $myProgOutput5";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput4); }

$cmd = "time $exeDir/findAndKillSuperReadsWithKUnitigOverlapsGeTheMerLength.perl $myProgOutput5 $merLen $myProgOutput6good $myProgOutput6bad";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput5); }

$cmd = "time cat $myProgOutput6good | $exeDir/reportSuperReadGroups.perl > $myProgOutput7";
print "$cmd\n"; system ($cmd);

if ($jumpLibraryReads) {
    goto jumpLibraryCalculations; }

$cmd = "$exeDir/createFastaSuperReadSequences $workingDirectory $myProgOutput7 -seqdiffmax $seqDiffMax -error-filename $sequenceCreationErrorFile1 -nosequence > $myProgOutput8";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput7); }

# The following generate the following files in $workingDirectory:
# 1) chimeric_read.txt: A list of reads with k-unitigs which overlap themselves
#      by too much (unless all matches have a small relative offset)
# 2) multiCopyKUnitigs.kUnitig.rptLength.txt: A list of k-unitigs which
#      have multiple copies overlapping a given read; don't use for merging
#      mates into super-reads.
$cmd = "cat $myProgOutput1_1 | $exeDir/getRecordsForTandomRepeatKUnitigs | $exeDir/reduceRecordsForKillingKUnitigsToConnectMates.perl $workingDirectory";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput1_1); }

# =====================================================================
# PUTTING IN NEW STUFF HERE (3/15/11)
$cmd = "time cat $myProgOutput6good | $exeDir/killUnwantedDataFromFirstSuperReadInfo.perl $workingDirectory > $myProgOutput10";
print "$cmd\n"; system ($cmd);

$cmd = "time cat $myProgOutput10 | $exeDir/makeListOfSuperReadsAndCounts.fromSuperReadGroupsFile.perl > $myProgOutput14";
print "$cmd\n"; system ($cmd);

$maxReadCountToEliminateSuperRead = $minReadsInSuperRead-1;
$cmd = "cat $myProgOutput14 | $exeDir/listInfrequentlyOccurringSuperReads.perl $maxReadCountToEliminateSuperRead > $myProgOutput15";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput14); }

$cmd = "time $exeDir/eliminateInfrequentlyOccurringSuperReadsUsingList.perl $myProgOutput10 $myProgOutput15 > $myProgOutput16";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput10); }

$cmd = "time $exeDir/getReadStartsOffsetsAndOrientationsInSuperReads_1stPass $myProgOutput3 $myProgOutput16 $kUnitigLengthsFile $sequenceCreationErrorFile1 > $myProgOutput18";
print "$cmd\n"; system ($cmd);

if ($joinMates) {
    $cmd = "time $exeDir/mergeShortMatesIntoSuperReads $myProgOutput16 $numKUnitigsFile $errorOutput9 $repetitiveKUnitigInfoFile $kUnitigLengthsFile $numReadsFile > $myProgOutput22"; }
else {
    $myProgOutput16complete = $myProgOutput16;
    if ($myProgOutput16complete !~ /^\//) {
	$myProgOutput16complete = "$pwd/$myProgOutput16complete"; }
    
    $cmd = "ln -s $myProgOutput16complete $myProgOutput22"; }
print "$cmd\n"; system ($cmd);

$cmd = "time $exeDir/reportFinalReadPlacementsInSuperReads.perl $kUnitigLengthsFile $myProgOutput18 $myProgOutput22 > $myProgOutput23";
print "$cmd\n"; system ($cmd);

if ($joinShooting) {
    chdir ($workingDirectory);
    $cmd = "$exeDir/postSuperReadPipelineCommandsForJoiningMates.perl $forceJoin $defaultMean $defaultStdev -l $merLen -kunitig-files-prefix $kUnitigFastaSequencePrefixComplete -read-placements-file $myProgOutput23complete";
    print "$cmd\n"; system ($cmd);
    chdir ($pwd);
    $cmd = "$exeDir/mergePostMateMergeAndPriorSuperReadGroupsByReadFiles.perl $myProgOutput24 $myProgOutput22 > $myProgOutput25";
    print "$cmd\n"; system ($cmd);
    if (! $mikedebug) { &killFiles ($myProgOutput22, $myProgOutput24); }
}
else {
    $myProgOutput22complete = $myProgOutput22;
    if ($myProgOutput22complete !~ /^\//) {
	$myProgOutput22complete = "$pwd/$myProgOutput22complete"; }
    $cmd = "ln -s $myProgOutput22complete $myProgOutput25";
    print "$cmd\n"; system ($cmd);
}
    
$cmd = "$exeDir/findEquivalentSuperReads.perl $myProgOutput25 -equiv-file $myProgOutput26 -out-file $myProgOutput27";
print "$cmd\n"; system ($cmd);

$cmd = "cat $myProgOutput27 | $exeDir/reportSuperReadGroups.perl > $myProgOutput28";
print "$cmd\n"; system ($cmd);

$cmd = "time $exeDir/createFastaSuperReadSequences $workingDirectory $myProgOutput28 -seqdiffmax $seqDiffMax -error-filename $sequenceCreationErrorFile2 > $finalSuperReadSequenceFile";
print "$cmd\n"; system ($cmd);

$cmd = "cat $sequenceCreationErrorFile1 $sequenceCreationErrorFile2 > $sequenceCreationErrorFileCombined";
print "$cmd\n"; system ($cmd);

$cmd = "$exeDir/getReadStartsOffsetsAndOrientationsInSuperReads_1stPass $myProgOutput3 $myProgOutput27 $kUnitigLengthsFile $sequenceCreationErrorFileCombined > $finalReadPlacementFileUsingReadNumbers";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput3); }

$cmd = "time $exeDir/changeReadNumsToReadNamesAtBeginOfLine $myProgOutput0_1 $finalReadPlacementFileUsingReadNumbers > $finalReadPlacementFile";
print "$cmd\n"; system ($cmd);

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
exit (0);

jumpLibraryCalculations:
$cmd = "$exeDir/createFastaSuperReadSequences $workingDirectory $myProgOutput7 -seqdiffmax $seqDiffMax  -error-filename $sequenceCreationErrorFile1 > $myProgOutput8";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput7); }

$cmd = "time $exeDir/getReadStartsOffsetsAndOrientationsInSuperReads_1stPass $myProgOutput3 $myProgOutput6good $kUnitigLengthsFile $sequenceCreationErrorFile1 > $myProgOutput12usingReadNumbers";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput3, $myProgOutput6good); }

$cmd = "time $exeDir/changeReadNumsToReadNamesAtBeginOfLine $myProgOutput0_1 $myProgOutput12usingReadNumbers > $myProgOutput12";
print "$cmd\n"; system ($cmd);

$cmd = "time $exeDir/outputSuperReadSeqForJumpLibrary.perl $myProgOutput8 $myProgOutput12 > $finalSuperReadSequenceFile";
print "$cmd\n"; system ($cmd);
if (! $mikedebug) { &killFiles ($myProgOutput12usingReadNumbers); }

exit(0);

# =====================================================================
# REDONE TO HERE

@guillaumeKUnitigCountFiles = @guillaumeKUnitigFiles;
for ($i=0; $i<=$#guillaumeKUnitigCountFiles; $i++) {
    $guillaumeKUnitigCountFiles[$i] =~ s/fa$/counts/; }
@filesToKill = (@guillaumeKUnitigFiles, $jellyfishKUnitigHashFile, @guillaumeKUnitigCountFiles, $totReadFile, $myProgOutput1, $myProgOutput2, $myProgOutput3, $myProgOutput4);
if ($clean) {
    print "Cleaning up...\n";
    for (@filesToKill) {
	$file = $_;
	unlink ($file); } }

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
	unlink ($file); }
}

