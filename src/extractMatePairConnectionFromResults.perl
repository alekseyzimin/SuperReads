#!/usr/bin/perl
$inputFileWithReadPairsToBeJoined = $ARGV[0];
$inputFileWithKUnitigPairsToBeJoined = $ARGV[1];
$matePairJoinerResultsFile = $ARGV[2];
$outputJoinFile = $ARGV[3];
$outputFileWithKUnitigsToKeep = $ARGV[4];

open (FILE1, $inputFileWithReadPairsToBeJoined);
open (FILE2, $inputFileWithKUnitigPairsToBeJoined);
open (INPUT_RESULTS_FILE, $matePairJoinerResultsFile);
open (OUTPUT_JOIN_FILE, ">$outputJoinFile");
open (OUTPUT_FILE_WITH_KUNIS_TO_KEEP, ">$outputFileWithKUnitigsToKeep");

$resultsLine = <INPUT_RESULTS_FILE>;
while ($tline = <INPUT_RESULTS_FILE>) {
    chomp ($tline);
    if ($tline !~ /^MATE/) {
	$resultsLine .= "$tline\n";
	next; }
    $readPairLine = <FILE1>; chomp ($readPairLine);
    @flds = split (" ", $readPairLine); $readPairLine = "@flds[0..3]";
    $kUnitigPairLine = <FILE2>; chomp ($kUnitigPairLine);
    # Now do the analysis for the resultsLine
    &reportResults ($resultsLine);
    # Now start over
    $resultsLine = "$tline\n";
}
$readPairLine = <FILE1>; chomp ($readPairLine);
@flds = split (" ", $readPairLine); $readPairLine = "@flds[0..3]";
$kUnitigPairLine = <FILE2>; chomp ($kUnitigPairLine);
# Now do the analysis for the resultsLine
&reportResults ($resultsLine);

close (FILE1);
close (FILE2);
close (INPUT_RESULTS_FILE);
close (OUTPUT_JOIN_FILE);
close (OUTPUT_FILE_WITH_KUNIS_TO_KEEP);

sub reportResults
{
    my ($resultsLines) = @_;
    my (@resultsLines, @kUnis, @oris, @numOvlsIn, @numOvlsOut);
    my ($resultsLine, @flds, $kUni, $ori, $numOvlsIn, $numOvlsOut);
    my ($firstBadLine, $lastBadLine, $i, $numPossibleLengths);

    chomp ($resultsLines);
    @resultsLines = split (/\n/, $resultsLines);
    return if ($#resultsLines < 3); # No connection possible
    @kUnis = @oris = @numOvlsIn = @numOvlsOut = ();
    $numPossibleLengths = 0;
    for (@resultsLines) {
	$resultsLine = $_;
	if ($resultsLine =~ /^\d+\s*$/) {
	    ++ $numPossibleLengths; }
	next unless ($resultsLine =~ /^uni/);
	@flds = split (" ", $resultsLine);
	($kUni) = ($flds[2] =~ /^(\d+)/);
	($ori) = ($flds[8] =~ /^(.)/);
	($numOvlsIn) = ($flds[17] =~ /^(\d+)/);
	($numOvlsOut) = ($flds[20] =~ /^(\d+)/);
	push (@kUnis, $kUni);
	push (@oris, $ori);
	push (@numOvlsIn, $numOvlsIn);
	push (@numOvlsOut, $numOvlsOut);
    }
    $firstBadLine = -1; $lastBadLine = -1;
    if ($numPossibleLengths > 1) {
	$firstBadLine = 0;
	$lastBadLine = $#numOvlsOut+1;
	goto outputABadPair; }
    for ($i=0; $i<=$#numOvlsOut; $i++) {
	if (($numOvlsOut[$i] > 1) || ($numOvlsIn[$i] > 1)) {
	    $firstBadLine = $i;
	    last; } }
    if ($firstBadLine < 0) { goto itsGood; }
    for ($i=$#numOvlsIn; $i>=0; $i--) {
	if (($numOvlsOut[$i] > 1) || ($numOvlsIn[$i] > 1)) {
	    $lastBadLine = $i;
	    last; } }
  outputABadPair:
    if ($lastBadLine <= $firstBadLine+1) { # Nothing acceptable
	print "NOTHING ACCEPTABLE (joins of too many possible lengths): $readPairLine ; $kUnitigPairLine\n";
	return; }
    print OUTPUT_FILE_WITH_KUNIS_TO_KEEP $readPairLine;
    print OUTPUT_FILE_WITH_KUNIS_TO_KEEP " $kUnis[0]";
    for ($i=$firstBadLine+1; $i<$lastBadLine; $i++) {
	print OUTPUT_FILE_WITH_KUNIS_TO_KEEP " $kUnis[$i]"; }
    print OUTPUT_FILE_WITH_KUNIS_TO_KEEP "\n";
    return;

    itsGood:
    print OUTPUT_JOIN_FILE $readPairLine;
    for ($i=0; $i<=$#kUnis; $i++) {
	if ($i > 0) { print OUTPUT_JOIN_FILE " 30"; }
	print OUTPUT_JOIN_FILE " $kUnis[$i] $oris[$i]"; }
    print OUTPUT_JOIN_FILE "\n";
}


