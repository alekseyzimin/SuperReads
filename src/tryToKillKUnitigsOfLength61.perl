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
$dir = "/genome6/raid/alekseyz/varroa_mite/Illumina_data/work1";
$file1 = "$dir/kUnitigLengths.txt";
$file2 = "$dir/superReadCounts.count.superRead.txt";


open (FILE, $file1);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    if ($flds[1] == 61) {
	$isNeeded{$flds[0]} = 1; }
}
close (FILE);

open (FILE, $file2);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $lineIsNeeded = 0;
    for ($i=1; $i<=$#flds; $i+=3) {
	if ($isNeeded{$flds[$i]}) {
	    $lineIsNeeded = 1;
	    break; }
    }
    if ($lineIsNeeded) {
	push (@lines, $line); }
}
close (FILE);

for (@lines) {
    print "$_\n";
}
