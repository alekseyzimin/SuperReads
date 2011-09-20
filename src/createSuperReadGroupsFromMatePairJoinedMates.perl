#!/usr/bin/perl
while ($line = <STDIN>) {
    chomp ($line);
    ($read1, $read2, $rest) = ($line =~ /^\s*(\S+)\S\s+(\S+)\S\s+\d+\s+\d+\s+(\S.*\S)\s*$/);
    @flds = split (" ", $rest);
    @flds2 = ();
    for ($i=0; $i<=$#flds; $i+=3) {
	$flds2[$#flds-1-$i] = $flds[$i]; }
    for ($i=1; $i<=$#flds; $i+=3) {
	if ($flds[$i] eq "F") {
	    $flds2[$#flds+1-$i] = "R"; }
	else {
	    $flds2[$#flds+1-$i] = "F"; } }
    for ($i=2; $i<=$#flds; $i+=3) {
	$flds2[$#flds-$i] = $flds[$i]; }
#    print "@flds\n";
#    print "@flds2\n";
    $lowerValuedRecord = 0;
    for ($i=0; $i<=$#flds; $i++) {
	if ($i % 3 != 1) {
	    if ($flds[$i] < $flds2[$i]) {
		$lowerValuedRecord = 1; }
	    elsif ($flds[$i] > $flds2[$i]) {
		$lowerValuedRecord = 2; } }
	else {
	    if ($flds[$i] lt $flds2[$i]) {
		$lowerValuedRecord = 1; }
	    elsif ($flds[$i] gt $flds2[$i]) {
		$lowerValuedRecord = 2; } }
	last if ($lowerValuedRecord > 0); }
#    if ($lowerValuedRecord != 0) {
#	print "Chosen record: "; }
    if ($lowerValuedRecord == 1) {
	$superReadName = "@flds"; }
#    elsif ($lowerValuedRecord == 2) {
    else {
	$superReadName = "@flds2"; }
    print "readName = $read1 : @flds\n$superReadName\n";
    print "readName = $read2 : @flds2\n$superReadName\n";
}

