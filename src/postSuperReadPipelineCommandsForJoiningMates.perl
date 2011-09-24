#!/usr/bin/perl
use File::Basename;
use Cwd;
$exeDir = dirname ($0);
$pwd = cwd;

$program0 = "$exeDir/createKUnitigMaxOverlaps";
$program1 = "$exeDir/getMateStatistics_pt1";
$program2 = "$exeDir/getMateStatistics_pt2.perl";
$program3 = "$exeDir/getStdevOfLibrary.perl";
$program11 = "$exeDir/convertReadMateInfoIntoKUnitigMateInfoForMateJoining.perl";
$program12 = "$exeDir/joinKUnitigs";
$program13 = "$exeDir/extractMatePairConnectionFromResults.perl";
$program21 = "$exeDir/createCombinedGoodJoinsFile.perl";
$program22 = "$exeDir/createSuperReadGroupsFromMatePairJoinedMates.perl";
&processCommandLine;

sub processCommandLine
{
    my ($i);

    $kmerLen = 31;
    $kUnitigFilesPrefix = "guillaumeKUnitigsAtLeast32bases";
    $readPlacementsInSuperReadsFile = "readPlacementsInSuperReads.final.read.superRead.offset.ori.usingReadNumbers.txt";
    $forceJoin = 0;
    for ($i=0; $i<=$#ARGV; $i++) {
	if ($ARGV[$i] eq "-kmer-len") {
	    ++$i;
	    $kmerLen = $ARGV[$i];
	    next; }
	if ($ARGV[$i] eq "-l") {
	    ++$i;
	    $kmerLen = $ARGV[$i];
	    next; }
	if ($ARGV[$i] eq "-kunitig-file") {
	    ++$i;
	    $kUnitigFilesPrefix = $ARGV[$i]; # This works if a file is specified as well
	    next; }
	if ($ARGV[$i] eq "-kunitig-files-prefix") {
	    ++$i;
	    $kUnitigFilesPrefix = $ARGV[$i];
	    next; }
	if ($ARGV[$i] eq "-read-placements-file") {
	    ++$i;
	    $readPlacementsInSuperReadsFile = $ARGV[$i];
	    next; }
	if ($ARGV[$i] eq "-default-mean") {
	    ++$i;
	    $defaultMean = $ARGV[$i];
	    next; }
	if ($ARGV[$i] eq "-default-stdev") {
	    ++$i;
	    $defaultStdev = $ARGV[$i];
	    next; }
	if ($ARGV[$i] eq "-force-join") {
	    $forceJoin = 1;
	    next; }
	if ($ARGV[$i] eq "-mikedebug") {
	    $mikedebug = 1;
	    next; }
	
    }
}
if ($forceJoin) {
    if (! $defaultMean) {
	$defaultMean = 300; }
    if (! $defaultStdev) {
	$defaultStdev = $defaultMean * .10; }
}

# INPUT FILES
$readToKUnitigMatchFile = "testOutput.nucmerLinesOnly.txt";
$kUnitigLengthsFile = "kUnitigLengths.txt";
$numReadsPerPrefixFile = "numReadsPerPrefix.txt";
# FILES CREATED DURING THE PIPELINE
$numKUnitigsFile = "numKUnitigs.txt";
$completeInputReadDistInfo = "prelimInsertLenInfo.readNum_insertLen_distToEndOfSuperRead_numKUnisInSuperRead.txt";
push (@filesToKill, $completeInputReadDistInfo);
$readMatesToBeJoined1 = "readMatesToBeJoinedWMeanAndStdev.txt";
push (@filesToKill, $readMatesToBeJoined1);
$kUnitigPairsForReadMatesToBeJoined1 = "kUnitigPairsToBeJoinedDueToReadMates.wMeanAndStdev.txt";
push (@filesToKill, $kUnitigPairsForReadMatesToBeJoined1);
$matePairJoiningResultsFile1 = "kUnitigPairResults.txt";
push (@filesToKill, $matePairJoiningResultsFile1);
$prefixForOverlapsBetweenKUnitigs = "overlap"; # Prefix for the end-to-end ovls
$kmerLenMinus1 = $kmerLen - 1;
$wellJoinedMatePairFile1 = "goodJoinsFile.001.txt";
push (@filesToKill, $wellJoinedMatePairFile1);
$readMatesToBeJoinedWUniInfo1 = "failingMatesWUnisToUse.001.txt";
push (@filesToKill, $readMatesToBeJoinedWUniInfo1);
$kUnitigPairsForReadMatesToBeJoined2 = "kUnitigPairsToBeJoinedDueToReadMates.pass2.wMeanAndStdev.txt";
push (@filesToKill, $kUnitigPairsForReadMatesToBeJoined2);
$matePairJoiningResultsFile2 = "kUnitigPairResults.pass2.txt";
push (@filesToKill, $matePairJoiningResultsFile2);
$wellJoinedMatePairFile2 = "goodJoinsFile.002.txt";
push (@filesToKill, $wellJoinedMatePairFile2);
$readMatesToBeJoinedWUniInfo2 = "failingMatesWUnisToUse.002.txt";
push (@filesToKill, $readMatesToBeJoinedWUniInfo2);
$kUnitigPairsForReadMatesToBeJoined3 = "kUnitigPairsToBeJoinedDueToReadMates.pass3.wMeanAndStdev.txt";
push (@filesToKill, $kUnitigPairsForReadMatesToBeJoined3);
$matePairJoiningResultsFile3 = "kUnitigPairResults.pass3.txt";
push (@filesToKill, $matePairJoiningResultsFile3);
$wellJoinedMatePairFile3 = "goodJoinsFile.003.txt";
push (@filesToKill, $wellJoinedMatePairFile3);
$readMatesToBeJoinedWUniInfo3 = "failingMatesWUnisToUse.003.txt";
push (@filesToKill, $readMatesToBeJoinedWUniInfo3);
# The next two files are created by $program11 on passes 2 and 3 respectively
$mateJoinerPrefixFile = "mateJoinerPrefixFile.txt";
$mateJoinerSuffixFile = "mateJoinerSuffixFile.txt";
push (@filesToKill, $mateJoinerPrefixFile, $mateJoinerSuffixFile);
#
$combinedWellJoinedMatePairFile = "goodJoinsFile.combined.txt";
push (@filesToKill, $combinedWellJoinedMatePairFile);
$mateJoinedSuperReadGroupForEachReadFile = "superReadGroupsForEachReadWhichHasAGroup.overlapJoinedMates.txt";

# goto mainLoop;
# Pre-programs:

open (FILE, $kUnitigLengthsFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $kUniLen[$flds[0]] = $flds[1]; }
close (FILE);
for ($kUni=0; $kUni<=$#kUniLen; $kUni++) {
    if (! $kUniLen[$kUni]) {
	$kUniLen[$kUni] = 0; } }
$numKUnitigs = $#kUniLen + 1;
open (OUTFILE, ">$numKUnitigsFile");
print OUTFILE "$numKUnitigs\n";
close (OUTFILE);

$cmd = "$program0 $kUnitigFilesPrefix -kmervalue $kmerLen $prefixForOverlapsBetweenKUnitigs -largest-kunitig-number $numKUnitigs";
print "$cmd\n"; $cmd = "time $cmd\n"; system ($cmd);

$cmd = "$program1 $kUnitigLengthsFile $readPlacementsInSuperReadsFile > $completeInputReadDistInfo";
print "$cmd\n";
$cmd = "time $cmd\n"; system ($cmd);

$cmd = "$program2 $numReadsPerPrefixFile $completeInputReadDistInfo";
print "$cmd\n";
$cmd = "time $cmd\n"; system ($cmd);

$cmd = "cat $numReadsPerPrefixFile | wc -l > numDiffPrefixes.txt";
print "$cmd\n";
$cmd = "time $cmd\n"; system ($cmd);
push (@filesToKill, "numDiffPrefixes.txt");

open (FILE, "numDiffPrefixes.txt");
$numDiffPrefixes = <FILE>; close (FILE); chomp ($numDiffPrefixes);

for ($i=1; $i<=$numDiffPrefixes; $i++) {
    $infile = "insertSizes.lib${i}.txt";
    next unless (-e $infile);
    push (@filesToKill, $infile);
    $cmd = "cat $infile | $program3 -stdev 4 > lib${i}.meanAndStdev.txt";
    print "$cmd\n";
    $cmd = "time $cmd\n"; system ($cmd);
}

$kUnitigOverlapsFile = "${prefixForOverlapsBetweenKUnitigs}.overlaps";
open (FILE, $readPlacementsInSuperReadsFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $super[$flds[0]] = $flds[1];
}
close (FILE);
$firstReadForLib[0] = 0;
open (FILE, $numReadsPerPrefixFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    push (@firstReadForLib, $flds[2]+$flds[1]);
}
close (FILE);
for ($i=1; $i<=$numDiffPrefixes; $i++) {
    $infile = "lib${i}.meanAndStdev.txt";
    if (! -e $infile) {
	next unless ($forceJoin);
	$mean[$i-1] = $defaultMean;
	$stdev[$i-1] = $defaultStdev;
	next; }

    open (FILE, $infile);
    $line = <FILE>;
    if ($line =~ /^Initial:/) {
	while ($line = <FILE>) {
	    last if ($line =~ /^Final:/); }
	chomp ($line);
	@flds = split (" ", $line);
	$line = "$flds[-2] $flds[-1]"; # To make the same anal below
    }
    close (FILE);
    chomp ($line);
    @flds = split (" ", $line);
    $mean[$i-1] = int ($flds[0]+.5);
    $stdev[$i-1] = int ($flds[1]+.5);
}
    
open (OUTFILE, ">$readMatesToBeJoined1");
for ($i=0;$i<=$#super; $i+=2) {
    $ip1 = $i+1;
    next unless ($super[$i]);
    next unless ($super[$ip1]);
    next if ($super[$i] eq $super[$ip1]);
    for ($j=0; $j<$#firstReadForLib; $j++) {
	last if ($firstReadForLib[$j+1] > $i); }
    next unless ($mean[$j]);
    print OUTFILE $i,"F ",$ip1,"R $mean[$j] $stdev[$j]\n";
}
close (OUTFILE);

mainLoop:
for ($pass=1; $pass<=3; $pass++) {
    if ($pass == 1) { $readMatesToBeJoinedFile = $readMatesToBeJoined1; }
    elsif ($pass == 2) { $readMatesToBeJoinedFile = $readMatesToBeJoinedWUniInfo1; }
    elsif ($pass == 3) { $readMatesToBeJoinedFile = $readMatesToBeJoinedWUniInfo2; }
    $kUnitigPairsForReadMatesToBeJoinedFile = ${"kUnitigPairsForReadMatesToBeJoined$pass"};
    $matePairJoiningResultsFile = ${"matePairJoiningResultsFile$pass"};
    $wellJoinedMatePairFile = ${"wellJoinedMatePairFile$pass"};
    $readMatesToBeJoinedWUniInfoOutputFile = ${"readMatesToBeJoinedWUniInfo$pass"};
    $cmd = "$program11 $readMatesToBeJoinedFile $readToKUnitigMatchFile $pass > $kUnitigPairsForReadMatesToBeJoinedFile";
    print "$cmd\n";
    $cmd = "time $cmd\n"; system ($cmd);
    
    $cmd = "$program12 $kUnitigOverlapsFile $kUnitigPairsForReadMatesToBeJoinedFile -min-overlap-length $kmerLenMinus1 -report-paths > $matePairJoiningResultsFile";
    print "$cmd\n";
    $cmd = "time $cmd\n"; system ($cmd);
    
    $cmd = "$program13 $readMatesToBeJoinedFile $kUnitigPairsForReadMatesToBeJoinedFile $matePairJoiningResultsFile $wellJoinedMatePairFile $readMatesToBeJoinedWUniInfoOutputFile";
    print "$cmd\n";
    $cmd = "time $cmd\n"; system ($cmd);
}

$cmd = "$program21 $mateJoinerPrefixFile $mateJoinerSuffixFile $wellJoinedMatePairFile1 $wellJoinedMatePairFile2 $wellJoinedMatePairFile3 > $combinedWellJoinedMatePairFile";
print "$cmd\n";
$cmd = "time $cmd\n"; system ($cmd);

$cmd = "cat $combinedWellJoinedMatePairFile | $program22 > $mateJoinedSuperReadGroupForEachReadFile";
print "$cmd\n";
$cmd = "time $cmd\n"; system ($cmd);

if (! $mikedebug) {
    &killFiles (@filesToKill); }

sub killFiles
{
    my (@filesToKill) = @_;
    my ($file);

    for (@filesToKill) {
        $file = $_;
        unlink ($file); }
}

