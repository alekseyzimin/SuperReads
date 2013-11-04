#! /usr/bin/env perl
# This expects 1 arg: the name of an asm file
# Output is to STDOUT
# Unitig# [UNS], where 'U' is U-unitig, 'S' is surrogate, 'N' is degenerate
$fn = $ARGV[0];
$cmd = "grep \"\\\{UTG\" -A 10 $fn |";
open (FILE, $cmd);
while ($line = <FILE>) {
    if ($line =~ /^acc/) {
	($uniNum) = ($line =~ /^acc..(\d+)\D/); }
    if ($line =~ /^sta:/) {
	($code) = ($line =~ /^sta.(.)/);
	print "$uniNum $code\n"; }
}
close (FILE);

