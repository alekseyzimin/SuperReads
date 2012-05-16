#!/usr/bin/env perl
# Pipe in readPositionsInSuperReads
while ($line = <STDIN>) {
    chomp ($line);
    @flds = split (" ", $line);
    ($val) = ($flds[0] =~ /^..(\d+)/);
    $val = int ($val/2);
    $lineForMatch = "$val $flds[1]";
    if ($lineForMatch eq $lineForMatchHold) {
	$val *= 2;
	print "cc$val\ncc", $val+1, "\n"; }
    $lineForMatchHold = $lineForMatch;
}

