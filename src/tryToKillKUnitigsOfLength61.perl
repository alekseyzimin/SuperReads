#!/usr/bin/perl
$dir = "/genome6/raid/alekseyz/varroa_mite/Illumina_data/work1";
$file1 = "$dir/kUnitigLengths.txt";
$file2 = "$dir/superReadCounts.count.superRead.txt";


open (FILE, $file1);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    if ($flds[1] == 61) {
	$isNeeded{$flds[0]} = 1; }
}
close (FILE);

open (FILE, $file2);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $lineIsNeeded = 0;
    for ($i=1; $i<=$#flds; $i+=3) {
	if ($isNeeded{$flds[$i]}) {
	    $lineIsNeeded = 1;
	    break; }
    }
    if ($lineIsNeeded) {
	push (@lines, $line); }
}
close (FILE);

for (@lines) {
    print "$_\n";
}
