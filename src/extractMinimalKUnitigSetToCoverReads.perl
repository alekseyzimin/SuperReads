#!/usr/bin/perl
if ($#ARGV >= 0) {
    $inputIsFromAFile = 1;
    $inputFile = $ARGV[0]; 
    open (FILE, $inputFile); }

while (1) {
    if ($inputIsFromAFile) { last if (!($line = <FILE>)); }
    else { last if (!($line = <STDIN>)); }
    if ($line =~ /^readNum/) {
	$hdrLine = $line;
	next; }
    elsif ($line =~ /\S/) {
	chomp ($line);
	push (@lines, $line);
	next; }
    $lastLineIndex = $#lines;
    $numLines = $#lines+1;
    if ($numLines == 1) {
	print $hdrLine, $lines[0], "\n";
	print "\n";
	@lines = ();
	next; }

    @lines2 = ();
    $maxReadOffsetSeen = -1;
    for ($i=0; $i<$numLines; $i++) {
	($readOffset1, $readOffset2) = ($lines[$i] =~ /^\s*\S+\s+\S+\s+(\S+)\s+(\S+)\s/);
	if ($readOffset1 < $readOffset2) {
	    $minReadOffset = $readOffset1;
	    $maxReadOffset = $readOffset2; }
	else {
	    $minReadOffset = $readOffset2;
	    $maxReadOffset = $readOffset1; }
#	print $lines[$i]," min = $minReadOffset ; max = $maxReadOffset ; maxReadOffsetSeen = $maxReadOffsetSeen\n";
	next if ($maxReadOffset <= $maxReadOffsetSeen);
#	print "Got here\n";
	push (@lines2, $lines[$i]);
#	print "lines2 = @lines2\n";
	$maxReadOffsetSeen = $maxReadOffset; }
    @lines = @lines2;
#    print "Final lines2 = @lines\n";
    $lastLineIndex = $#lines;
    $numLines = $#lines+1;

    if ($numLines == 1) {
	print $hdrLine, $lines[0], "\n";
	print "\n";
	@lines = ();
	next; }

    ($readOffset1, $readOffset2) = ($lines[0] =~ /^\s*\S+\s+\S+\s+(\S+)\s+(\S+)\s/);
    $lastBaseContainedInFirstKUnitig = $readOffset1;
    if ($readOffset2 > $lastBaseContainedInFirstKUnitig) {
	$lastBaseContainedInFirstKUnitig = $readOffset2; }

    ($readOffset1, $readOffset2) = ($lines[$lastLineIndex] =~ /^\s*\S+\s+\S+\s+(\S+)\s+(\S+)\s/);
    $firstBaseContainedInLastKUnitig = $readOffset1;
    if ($readOffset2 < $firstBaseContainedInLastKUnitig) {
	$firstBaseContainedInLastKUnitig = $readOffset2; }

    if ($firstBaseContainedInLastKUnitig <= $lastBaseContainedInFirstKUnitig+1) {
	print $hdrLine, $lines[0], "\n", $lines[$lastLineIndex], "\n\n";
	@lines = ();
	next; }
    
    print $hdrLine, $lines[0], "\n";
    $minSpclKUnitigNum = 1000000000000;
    for ($i=1; $i<$lastLineIndex; $i++) {
	($readOffset1, $readOffset2) = ($lines[$i] =~ /^\s*\S+\s+\S+\s+(\S+)\s+(\S+)\s/);
	if ($readOffset1 < $readOffset2) {
	    $minReadOffset = $readOffset1;
	    $maxReadOffset = $readOffset2; }
	else {
	    $minReadOffset = $readOffset2;
	    $maxReadOffset = $readOffset1; }
	next if ($minReadOffset > $lastBaseContainedInFirstKUnitig + 1);
	next if ($maxReadOffset < $firstBaseContainedInLastKUnitig - 1);
#	print $lines[$i],"\n"; # For debugging
	($kUnitigNum) = ($lines[$i] =~ /\s(\S+)\s+\S+\s*$/);
	if ($kUnitigNum < $minSpclKUnitigNum) {
	    $lineIndexOfCoveringKUnitig = $i;
	    $minSpclKUnitigNum = $kUnitigNum; } }
#    print "Picked line = " if ($minSpclKUnitigNum < 1000000000000); # For debugging
    if ($minSpclKUnitigNum < 1000000000000) { 
	print "$lines[$lineIndexOfCoveringKUnitig]\n"; }
    else {
	for ($i=1; $i<$lastLineIndex; $i++) {
	    print $lines[$i], "\n"; } }
    print $lines[$lastLineIndex], "\n";
    print "\n";
    @lines = (); }
