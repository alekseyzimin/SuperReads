# SuperRead pipeline
# Copyright (C) 2012  Genome group at University of Maryland.
# 
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



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



	
