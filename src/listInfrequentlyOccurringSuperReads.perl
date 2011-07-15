#!/usr/bin/perl
# This exec expects 1 arg: the threshold number of occurrences which should
# be the maximum on the list; an arg of 5 would list all super-reads with
# a count of 5 or less
# Input on STDIN and should be one line for each super-read:
#  The number of times it occurs followed by the name of the super-read
$maxCount = $ARGV[0];
while ($line = <STDIN>) {
    chomp ($line);
    ($curCount, $superRead) = ($line =~ /^(\d+)\s+(\S.*\S)\s*$/);
    if ($curCount <= $maxCount) {
	print "$superRead\n"; }
}
