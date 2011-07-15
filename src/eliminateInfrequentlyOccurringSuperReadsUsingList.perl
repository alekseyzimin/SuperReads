#!/usr/bin/perl
# This exec expects 2 args:
# The name of the input file containing reads and their containing super-reads
# The name of a file listing super-reads to exclude
# Output is to STDOUT
$readsPlusSuperReadsFile = $ARGV[0];
$fileOfSuperReadsToExclude = $ARGV[1];
open (FILE, $fileOfSuperReadsToExclude);
while ($superRead = <FILE>) {
    chomp ($superRead);
    $killIt{$superRead} = 1; }
close (FILE);

open (FILE, $readsPlusSuperReadsFile);
while ($line = <FILE>) {
    $superRead = <FILE>;
    chomp ($superRead);
    next if ($killIt{$superRead});
    print $line, $superRead, "\n";
}
close (FILE);

