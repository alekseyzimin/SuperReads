#!/usr/bin/perl
$reducedKUnitigMatchesFile = $ARGV[0];
$superReadGroupFile = $ARGV[1];
$kUnitigLengthsFile = $ARGV[2];
if ($#ARGV >= 3) {
    $superReadsEliminatedDueToBadKUnitigMatchFile = $ARGV[3]; }

if ($superReadsEliminatedDueToBadKUnitigMatchFile) {
    open (FILE, $superReadsEliminatedDueToBadKUnitigMatchFile);
    while ($line = <FILE>) {
	($superReadName) = ($line =~ /^(\S+)\s/);
	$isBadSuperRead{$superReadName} = 1;
	$line = <FILE>;
	$line = <FILE>; }
    close (FILE); }

open (FILE, $superReadGroupFile);
while ($line = <FILE>) {
    chomp ($line);
    ($readName, $origSuperReadNameInfo) = ($line =~ /readName = (\S+) :\s*(\S.+\S)\s*$/);
    $newSuperReadNameInfo = <FILE>;
    chomp ($newSuperReadNameInfo);
    $wasReversed = 0;
    if ($origSuperReadNameInfo ne $newSuperReadNameInfo) {
	$wasReversed = 1; }
    @newSuperReadNameInfoFlds = split (" ", $newSuperReadNameInfo);
    $superReadName = $newSuperReadNameInfoFlds[0] . $newSuperReadNameInfoFlds[1];
    for ($i=2; $i<$#newSuperReadNameInfoFlds; $i+=3) {
	$superReadName .= ("_" . $newSuperReadNameInfoFlds[$i] . "_" . $newSuperReadNameInfoFlds[$i+1] . $newSuperReadNameInfoFlds[$i+2]); }
    next if ($isBadSuperRead{$superReadName});
    $wasReversed[$readName] = $wasReversed;
    $newSuperReadNameInfo[$readName] = $newSuperReadNameInfo;
}
close (FILE);

open (FILE, $kUnitigLengthsFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $kUnitigLengths[$flds[0]] = $flds[1];
}
close (FILE);

    
open (FILE, $reducedKUnitigMatchesFile);
$isFirstLine = 1;
while ($line = <FILE>) {
    chomp ($line);
    if ($line =~ /readNum/) {
	$isFirstLine = 1;
	next; }
    if ($line !~ /\S/) {
	goto processTheRead; }
    @flds = split (" ", $line);
    $kUniNum = $flds[6];
    $kUniLen = $flds[4];
    $ahg = $flds[0]-1;
    $bhg = $kUniLen - $flds[1];
    if ($flds[2] < $flds[3]) {
	$ahg -= ($flds[2]-1);
	$bhg = ($flds[5] - $flds[3]) - $bhg;
    }
    else {
	$ahg -= ($flds[5]-$flds[2]);
	$bhg = ($flds[3]-1) - $bhg;
    }
#    print "ahg = $ahg ; bhg = $bhg ; line = $line\n";
    $kUnitigLengths[$kUniNum] = $kUniLen;
    if ($isFirstLine) {
	$isFirstLine = 0;
	$readName = $flds[7];
	$wasReversed = $wasReversed[$readName];
	@newSuperReadNameInfoFlds = split (" ", $newSuperReadNameInfo[$readName]);
	$superReadLength = $kUnitigLengths[$newSuperReadNameInfoFlds[0]];
	$superReadName = $newSuperReadNameInfoFlds[0] . $newSuperReadNameInfoFlds[1];
	for ($i=2; $i<$#newSuperReadNameInfoFlds; $i+=3) {
	    $superReadName .= ("_" . $newSuperReadNameInfoFlds[$i] . "_" . $newSuperReadNameInfoFlds[$i+1] . $newSuperReadNameInfoFlds[$i+2]);
	    $superReadLength += ($kUnitigLengths[$newSuperReadNameInfoFlds[$i+1]] - $newSuperReadNameInfoFlds[$i]);
	}
#	if ($wasReversed == 0) {
	    $ahgHold = $ahg;
	    $bhgHold = $bhg;
#	}
    }
#    if ($wasReversed != 0) {
#	$ahgHold = $ahg;
#	$bhgHold = $bhg; }
    next;

  processTheRead:
    next if ($isBadSuperRead);
    if ($wasReversed[$readName] == 0) {
	if ($newSuperReadNameInfoFlds[1] eq "F") {
	    $firstBaseOffset = $ahgHold; }
	else {
	    $firstBaseOffset = - $bhgHold; } }
    else {
	if ($newSuperReadNameInfoFlds[-1] eq "R") {
	    $firstBaseOffset = $superReadLength - $ahgHold; }
	else {
	    $firstBaseOffset = $superReadLength + $bhgHold; }
    }
#    print "wasReversed = $wasReversed[$readName] ; ahgHold = $ahgHold ; bhgHold = $bhgHold\n";
    if ($wasReversed[$readName]) {
	$readOriInSuperRead = "R"; }
    else {
	$readOriInSuperRead = "F"; }
    if ($superReadName) {
	$firstBaseOffset *= 1;
	print "$readName $superReadName $firstBaseOffset $readOriInSuperRead\n";
    }
}
close (FILE);

