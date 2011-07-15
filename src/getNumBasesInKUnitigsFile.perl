#!/usr/bin/perl
# Cat the k-unitig fasta file through
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
    if (! $length) { $length = 0; }
    print "$i $length\n";
}
