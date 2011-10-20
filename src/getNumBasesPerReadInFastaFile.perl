#!/usr/bin/perl
# Pass the name of a fasta file as an arg and this gives the number of
# bases in each read of the fasta file
if ($#ARGV != 0) {
    open (FILE, $0);
    while ($line = <FILE>) {
	last unless ($line =~ /^\#/);
	print $line;
    }
    close (FILE);
    exit;
}
$file = $ARGV[0];
$cmd = "zcat -f $file |";
open (FILE, $cmd);
$isFirstRead = 1;
while ($line = <FILE>) {
    if ($line =~ /^>/) {
	if (! $isFirstRead) { print "$readLen\n"; }
	$readLen = 0;
	$isFirstRead = 0;
    }
    else {
	$len = length ($line)-1;
	$readLen += $len;
    }
}
if (! $isFirstRead) { print "$readLen\n"; }
