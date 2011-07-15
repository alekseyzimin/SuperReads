#!/usr/bin/perl
while ($line = <STDIN>) {
    $line = <STDIN>;
    if (! $count{$line}) {
	push (@superReads, $line); }
    ++$count{$line}; }
# @superReads = keys %count;
for (@superReads) {
    $superRead = $_;
    print $count{$superRead}, " ", $superRead;
}
