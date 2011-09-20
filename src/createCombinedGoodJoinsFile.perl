#!/usr/bin/perl
$prefixFile = $ARGV[0];
$suffixFile = $ARGV[1];
@goodJoinsFiles = (@ARGV[2..4]);

open (FILE, $prefixFile);
while ($line = <FILE>) {
    chomp ($line);
    next unless ($line =~ /\S\s+\S/);
    ($read, $prefix) = ($line =~ /^(\d+)\s+(\S.+\S)\s*$/);
    $prefix[$read] = $prefix;
}
close (FILE);

open (FILE, "$suffixFile");
while ($line = <FILE>) {
    chomp ($line);
    next unless ($line =~ /\S\s+\S/);
    ($read, $suffix) = ($line =~ /^(\d+)\s+(\S.+\S)\s*$/);
    $suffix[$read] = $suffix;
}
close (FILE);

for (@goodJoinsFiles) {
    $goodJoinFile = $_;
    open (FILE, $goodJoinFile);
    while ($line = <FILE>) {
	chomp ($line);
	($begin, $rest) = ($line =~ /^(\S+\s+\S+\s+\S+\s+\S+\s+)(\S.+\S)\s*$/);
	($read1, $read2) = ($begin =~ /^(\d+)\S\s+(\d+)\S\s/);
	print $begin;
	print "$prefix[$read1] " if ($prefix[$read1]);
	print $rest;
	print " $suffix[$read2]" if ($suffix[$read2]);
	print "\n";
    }
}



	
