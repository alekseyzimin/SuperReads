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
$kUni = $ARGV[-1];
$cmd = "grep --text -10 \" $kUni \" Illumina_data/work1/arrangedCoordsResultsByRead.txt |";
open (FILE, $cmd);
while ($line = <FILE>) {
    chomp ($line);
    if ($line =~ /^readNum/) {
	($readNum) = ($line =~ /^readNum = (\d+)\D/);
	$lines{$readNum} .= "$line\n";
	next; }
    if ($line !~ /\S/) {
	$lines{$readNum} .= "\n";
	next; }
    $lines{$readNum} .= "$line\n";
    ++$lineCount{$readNum};
    @flds = split (" ", $line);
    if ($flds[6] == $kUni) {
	$isNeeded{$readNum} = 1;
	if ($flds[2] < $flds[3]) { $minOffset = $flds[2]; $maxOffset = $flds[3]; }
	else { $minOffset = $flds[3]; $maxOffset = $flds[2]; }
	if (($minOffset == 1) && ($maxOffset == $flds[5])) {
	    $isCovered{$readNum} = 1; } }
}

@readNames = keys %lines;
@readNames = sort byNum @readNames;
for (@readNames) {
    $readName = $_;
    next unless ($isNeeded{$readName});
    print "readName = \"$readName\"; ";
    if ($isNeeded{$readName}) { print "isNeeded "; }
    print "lineCount = $lineCount{$readName}; ";
    print "isCovered = $isCovered{$readName}\n";
    print $lines{$readName};
}

sub byNum
{
    return ($a <=> $b);
}

