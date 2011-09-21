#!/usr/bin/perl
$readNumsPerPrefixFile = $ARGV[0];
$insertSizeFile = $ARGV[1];
$minInsertsToStudy = 400;
$numInsertsToStudy = 20000;
# Outputs "insertSizes.lib#.txt"; # starts from 1
open (FILE1, $readNumsPerPrefixFile);
open (FILE2, $insertSizeFile);
$line2 = <FILE2>;
chomp ($line2);
$outputFileNumber = 0;
while ($line = <FILE1>) {
    chomp ($line);
    ++$outputFileNumber;
    @flds = split (" ", $line);
    $first = $flds[2];
    $last = $first + $flds[1];
    $readPrefixLetters = $flds[0];
    @distsToEnd = ();
    @insertSizes = ();
    @numOfEachDistToEnd = ();
    @flds2 = split (" ", $line2);
    if ($flds2[1] > -1000) {
	push (@distsToEnd, $flds2[2]);
	push (@insertSizes, $flds2[1]);
	++$numOfEachDistToEnd[$flds2[2]]; }
    $atEof = 1;
    while ($line2 = <FILE2>) {
	$atEof = 0;
	@flds2 = split (" ", $line2);
	last if ($flds2[0] >= $last);
	if ($flds2[1] > -1000) {
	    push (@distsToEnd, $flds2[2]);
	    push (@insertSizes, $flds2[1]);
	    ++$numOfEachDistToEnd[$flds2[2]];
	}
    }
    $insertsAccumulated = 0;
    for ($i=$#numOfEachDistToEnd; $i>=0; $i--) {
	$insertsAccumulated += $numOfEachDistToEnd[$i];
	last if ($insertsAccumulated >= $minInsertsToStudy); }
    $minDistToConsider = $i;
    $numInsertsOverHalfLength = 0;
    for ($i=0; $i<=$#distsToEnd; $i++) {
	next unless ($distsToEnd[$i] >= $minDistToConsider);
	if ($insertSizes[$i] > $minDistToConsider / 2) {
	    ++ $numInsertsOverHalfLength; }
    }
    if ($numInsertsOverHalfLength > $insertsAccumulated * .1) {
	print STDERR "Read prefix $readPrefixLetters insert size estimation failed. numInsertsOverHalfLength = ${numInsertsOverHalfLength}; insertsAccumulated = ${insertsAccumulated}; minDistToConsider = ${minDistToConsider}. Skipping!\n";
	goto endOfLoop; }
    $insertsAccumulated = 0;
    for ($i=$#numOfEachDistToEnd; $i>=0; $i--) {
	$insertsAccumulated += $numOfEachDistToEnd[$i];
	last if ($insertsAccumulated >= $numInsertsToStudy); }
    if ($i <= 0) {
	print STDERR "Read prefix $readPrefixLetters insert size estimation failed. Not enough spanned inserts. Skipping!\n";
	goto endOfLoop; }
    $minDistToAnalyze = $i;
    open (OUTFILE, ">insertSizes.lib${outputFileNumber}.txt");
    for ($i=0; $i<=$#distsToEnd; $i++) {
	next unless ($distsToEnd[$i] >= $minDistToAnalyze);
	print OUTFILE "$insertSizes[$i]\n";
    }
    close OUTFILE;

  endOfLoop:
    last if ($atEof);
}
close (FILE1);
close (FILE2);
