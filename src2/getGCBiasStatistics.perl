#!/usr/bin/perl
use File::Basename;
$exeDir = dirname ($0);
$defaultReadLength = 100;
$minFromEnd = 200;
$intervalLen = 200;
&parseArgs;
# $assemblyDirectory = $ARGV[0]; # /genome3/raid/alekseyz/rhodobacter/assembly2.1.2
$reportedPlacementsFile = "tplacements.txt";
$intermResultsFile = "coverageCountsByGC.intermFile.txt";
$finalResultsFile = "adjustmentFactorsForGCPct.txt";
$cmd = "$exeDir/getSuperReadPlacements.perl -dir $assemblyDirectory > $reportedPlacementsFile";
print "$cmd\n";
system ($cmd);
$cmd = "$exeDir/getATBiasInCoverageForIllumina.v2 --sequence-file localUnitigSequenceFile.fasta --default-read-length $defaultReadLength --placement-file $reportedPlacementsFile --interval-len $intervalLen --min-from-end $minFromEnd > $intermResultsFile";
print "$cmd\n";
system ($cmd);
$cmd = "cat $intermResultsFile | $exeDir/getMeanAndStdevByGCCount.perl --interval-len $intervalLen > $finalResultsFile";
print "$cmd\n";
system ($cmd);

sub parseArgs
{
    for ($i=0; $i<=$#ARGV; ++$i) {
	if ($ARGV[$i] eq "--assemblyDirectory") {
	    ++$i;
	    $assemblyDirectory = $ARGV[$i];
	    next; }
	if ($ARGV[$i] eq "--min-from-end") {
	    ++$i;
	    $minFromEnd = $ARGV[$i];
	    next; }
	if ($ARGV[$i] eq "--default-read-length") {
	    ++$i;
	    $defaultReadLength = $ARGV[$i];
	    next; }
    }
}
