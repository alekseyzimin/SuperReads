#!/usr/bin/perl
# The file names might be
# output.txt and reverseComplemented.jumpingReadFile.fa
$fileWithListOfConnectedDirectories = $ARGV[0];
$fileWithReadsCorrespondingToDirectories = $ARGV[1];
open (FILE, $fileWithListOfConnectedDirectories);
while ($line = <FILE>) {
    ($num) = ($line =~ /^(\d+)\s/);
    $kickOut[$num] = 1; }
close (FILE);

open (FILE, "grep \"^>\" $fileWithReadsCorrespondingToDirectories |");
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    ($val) = ($flds[0] =~ /^...(\d+)$/);
    $val = int ($val/2);
    if ($kickOut[$val]) {
	print "$flds[1]\n"; }
}
close (FILE);

