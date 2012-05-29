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
# Cat the superReadGroupsForEachReadWhichHasAGroup.txt file in and this
# outputs the superReadGroups with the read info for them to STDOUT
$isFirstLine = 1;
while ($line = <STDIN>) {
    chomp ($line);
    if ($isFirstLine) {
	@flds = split (" ", $line);
	if ($flds[0] eq "readNum") {
	    $hasReadNum = 1; }
	else {
	    $hasReadNum = 0; }
	$isFirstLine = 0; }
    if ($hasReadNum) {
	($readNum, $readName) = ($line =~ /^readNum = (\d+); readName = (\S+)\s/); }
    else {
	($readName) = ($line =~ /^readName = (\S+)\s/);
	$readNum = $readName; }
    $superReadGroup = <STDIN>;
    if (! $readsInSuperReadGroup{$superReadGroup}) {
	push (@superReadGroups, $superReadGroup); }
    else {
	$readsInSuperReadGroup{$superReadGroup} .= " ; "; }
    $readsInSuperReadGroup{$superReadGroup} .= "$readNum $readName";
}

for (@superReadGroups) {
    $superReadGroup = $_;
    $readsInSuperReadGroup = $readsInSuperReadGroup{$superReadGroup};
    chomp ($superReadGroup);
    print "$superReadGroup : $readsInSuperReadGroup\n"; }
