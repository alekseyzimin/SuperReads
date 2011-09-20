#!/usr/bin/perl
$readMateFile = $ARGV[0];
$kUnitigVsReadNucmerFile = $ARGV[1];
$pass2file = "mateJoinerPrefixFile.txt";
$pass3file = "mateJoinerSuffixFile.txt";
$pass = 1;
if ($#ARGV >= 2) {
    $pass = $ARGV[2]; }

open (FILE, $readMateFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    ($read1, $ori1) = ($flds[0] =~ /^(.+)(.)$/);
    if ($ori1 =~ /\d/) { $read1 .= $ori1; $ori1 = "F"; }
    ($read2, $ori2) = ($flds[1] =~ /^(.+)(.)$/);
    if ($ori2 =~ /\d/) { $read2 .= $ori2; $ori2 = "R"; }
    push (@read1s, $read1);
    $readOri{$read1} = $ori1;
    push (@read2s, $read2);
    $readOri{$read2} = $ori2;
    $insertSize{$read1} = $flds[2];
    $insertStdDev{$read1} = $flds[3];
    $isNeeded{$read1} = 1;
    $isNeeded{$read2} = 2;
    if ($pass > 1) {
	$kUniSeq{$read1} = $kUniSeq{$read2} = " @flds[5..$#flds] ";
	if ($pass > 2) {
	    $firstReadsKUniSeq{$read1} = $flds[4]; }
    }
}
close (FILE);

if ($pass == 2) {
    open (SPCL_OUTFILE, ">$pass2file"); }
elsif ($pass == 3) {
    open (SPCL_OUTFILE, ">$pass3file"); }
# For the moment we will assume that all read1s are 'F' and read2s are 'R'
open (FILE, $kUnitigVsReadNucmerFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    next unless ($isNeeded{$flds[-1]});
    $read = $flds[-1];
    next if ($wasFound{$read} && ($isNeeded{$read} >= $pass));
    $kUni = $flds[-2];
    if ($read ne $readHold) {
	# Do the prior output for $readHold
	if ($pass > 1) {
	    &analyzeNucmerLines ($readHold); }
	# Now continue
	$readHold = $read;
	@heldNucmerLines = ();
    }
    push (@heldNucmerLines, $line);
    if (($pass > 1) && $wasFound{$read}) {
	if ($pass - $isNeeded{$read} == 1) {
	    next unless (index ($kUniSeq{$read}, " $kUni ") >= 0); }
	if ($pass - $isNeeded{$read} == 2) {
	    next unless ($kUni eq $firstReadsKUniSeq{$read}); }
    }
    $kUni{$read} = $kUni;
    if ($flds[6] < $flds[7]) {
	$lengthAdjustment{$read} = $flds[4] - $flds[6]; }
    else {
	$lengthAdjustment{$read} = ($flds[4]-1) - $flds[6]; }
    if ($flds[6] < $flds[7]) {
	$kUniOri{$read} = $readOri{$read}; }
    else {
	if ($isNeeded{$read} == 1) { # The first read of a pair
	    $kUniOri{$read} = "R"; }
	else { # The second read of a pair
	    $kUniOri{$read} = "F"; }
    }
    $wasFound{$read} = 1;
}
close (FILE);
close (SPCL_OUTFILE);

for ($i=0; $i<=$#read1s; $i++) {
    $read1 = $read1s[$i];
    $read2 = $read2s[$i];
    $kUni1 = $kUni{$read1};
    $kUni2 = $kUni{$read2};
    $insertSize = $insertSize{$read1};
    $insertStdDev = $insertStdDev{$read1};
    $ori1 = $kUniOri{$read1};
    $ori2 = $kUniOri{$read2};
    $insertSize += ($lengthAdjustment{$read1} + $lengthAdjustment{$read2});
    print "${kUni1}$ori1 ${kUni2}$ori2 $insertSize $insertStdDev\n";
}

sub analyzeNucmerLines
{
    my ($localRead) = @_;
    my (@oris, @minOffsets, @maxOffsets, $line, @flds, $readOffset1, $readOffset2);
    my ($good, $i, $index, @overlaps, @kUnitigNums, $tOri);

    return unless ($localRead =~ /\S/);
    @kUnitigNums = @oris = @minOffsets = @maxOffsets = @overlaps = ();
    for (@heldNucmerLines) {
	$line = $_;
	@flds = split (" ", $line);
	$readOffset1 = $flds[6];
	$readOffset2 = $flds[7];
	$kUnitigNum = $flds[-2];
	push (@kUnitigNums, $kUnitigNum);
	if ($readOffset1 < $readOffset2) {
	    push (@oris, "F");
	    push (@minOffsets, $readOffset1);
	    push (@maxOffsets, $readOffset2); }
	else {
	    push (@oris, "R");
	    push (@minOffsets, $readOffset2);
	    push (@maxOffsets, $readOffset1); }
    }
    $good = 1;
    for ($i=1; $i<=$#minOffsets; $i++) {
	$overlaps[$i-1] = ($maxOffsets[$i-1]+1) - $minOffsets[$i];
	# The following test should never apply for what we're doing
	if ($overlaps[$i-1] < 0) {
	    $good = 0; }
    }
    print SPCL_OUTFILE $localRead;
    if ($pass == 2) {
	for ($i=0; 1; $i++) {
	    last if ($kUnitigNums[$i] eq $kUni{$localRead});
	    print SPCL_OUTFILE " $kUnitigNums[$i] $oris[$i] $overlaps[$i]"; }
    }
    elsif ($pass == 3) {
	$isStarted = 0;
	for ($i=$#kUnitigNums; $i>=0; $i--) {
	    if ($kUnitigNums[$i] eq $kUni{$localRead}) {
		$isStarted = 1;
		next; }
	    if ($isStarted) {
		if ($oris[$i] eq "F") { $tOri = "R"; }
		else { $tOri = "F"; }
		print SPCL_OUTFILE " $overlaps[$i] $kUnitigNums[$i] $tOri";
	    }
	}
    }
    print SPCL_OUTFILE "\n";
}

