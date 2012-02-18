#!/usr/bin/perl
# Cat the k-unitig fasta file through
$workingDir = ".";
$outputPrefix = "";
for ($i=0; $i<=$#ARGV; $i++) {
    $arg = $ARGV[$i];
    if ($arg eq "-output-prefix") {
	++$i;
	$outputPrefix = $ARGV[$i];
	if ($outputPrefix !~ /\.$/) { $outputPrefix .= "."; }
	next; }
    $workingDir = $arg;
}

$numKUnitigsFile = $outputPrefix . "numKUnitigs.txt";
$maxKUnitigNumberFile = $outputPrefix . "maxKUnitigNumber.txt";
$totBasesInKUnitigsFile = $outputPrefix . "totBasesInKUnitigs.txt";
$isFirstRead = 1;
while ($line = <STDIN>) {
    if ($line =~ /^>/) {
	if (! $isFirstRead) { $kUnitigLengths[$kUnitig] = $kUnitigLength; }
	$kUnitigLength = 0;
	$isFirstRead = 0;
	($kUnitig) = ($line =~ /^.(\S+)\s/);
    }
    else {
	$len = length ($line)-1;
	$kUnitigLength += $len;
    }
}
if (! $isFirstRead) { $kUnitigLengths[$kUnitig] = $kUnitigLength; }

for ($i=0; $i<=$#kUnitigLengths; $i++) {
    $length = $kUnitigLengths[$i];
    $totBasesInKUnitigs += $length;
    if (! $length) {
	$length = 0; }
    else {
	++$numKUnitigs; }
    print "$i $length\n";
}
open (OUTFILE, ">", "$workingDir/$numKUnitigsFile");
print OUTFILE "$numKUnitigs\n";
close (OUTFILE);
open (OUTFILE, ">", "$workingDir/$maxKUnitigNumberFile");
$arraySizeForKUnitigData = $#kUnitigLengths+1;
print OUTFILE "$arraySizeForKUnitigData\n";
close (OUTFILE);
open (OUTFILE, ">", "$workingDir/$totBasesInKUnitigsFile");
print OUTFILE "$totBasesInKUnitigs\n";
close (OUTFILE);
