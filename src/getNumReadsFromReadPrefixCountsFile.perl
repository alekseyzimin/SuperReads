#!/usr/bin/perl
# The numReadsPerPrefix.txt file is the arg
open (FILE, $ARGV[0]);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $numReads += $flds[1];
}
close (FILE);
print "$numReads\n";
