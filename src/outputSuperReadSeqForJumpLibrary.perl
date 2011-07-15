#!/usr/bin/perl
$minReadSequenceLength = 64;
$maxReadLength = 2047;
$superReadSequenceFile = $ARGV[0];
$readPlacementFile = $ARGV[1];

open (FILE, $superReadSequenceFile);
while ($seqName = <FILE>) {
    chomp ($seqName);
    ($seqName) = ($seqName =~ /^>(.+)$/);
    $superReadSequence{$seqName} = <FILE>;
    chomp ($superReadSequence{$seqName});
}
close (FILE);
@superReadNames = keys %superReadSequence;
for (@superReadNames) {
    $superReadName = $_;
    $superReadLength{$superReadName} = length ($superReadSequence{$superReadName}); }

open (FILE, $readPlacementFile);
while ($line = <FILE>) {
    chomp ($line);
    ($readName, $superReadName, $offset, $ori) = split (" ", $line);
    # Keep super-read on the 5' end if needed to use the read as per
    # conversation with Jim on 5/2/11 at 3pm
    $adjustedOffset = $offset;
    $superReadLength = $superReadLength{$superReadName};
    if ($ori eq "F") {
	if ($superReadLength - $offset < $minReadSequenceLength) {
	    $adjustedOffset = $superReadLength - $minReadSequenceLength; }
	if ($adjustedOffset < 0) {
	    $adjustedOffset = 0; } }
    else {
	if ($offset < $minReadSequenceLength) {
	    $adjustedOffset = $minReadSequenceLength; }
	if ($adjustedOffset > $superReadLength) {
	    $adjustedOffset = $superReadLength; } }
    if ($ori eq "F") {
	$outputString = substr ($superReadSequence{$superReadName}, $adjustedOffset); }
    else {
	$outputString = substr ($superReadSequence{$superReadName}, 0, $adjustedOffset);
	$outputString = reverse ($outputString);
	$outputString =~ tr/acgtACGT/tgcaTGCA/; }
    if ($outputString) {
	print ">$readName\n";
	if (length ($outputString) > $maxReadLength) {
	    $outputString = substr ($outputString, 0, $maxReadLength); }
	print "$outputString\n";
    }
}


