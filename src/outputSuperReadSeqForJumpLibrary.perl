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
$minReadSequenceLength = 64;
$maxReadLength = 2047;
$superReadSequenceFile = $ARGV[0];
$readPlacementFile = $ARGV[1];

open (FILE, $superReadSequenceFile);
while ($seqName = <FILE>) {
    chomp ($seqName);
    ($seqName) = ($seqName =~ /^>(.+)$/);
    $superReadSequence{$seqName} = <FILE>;
    chomp ($superReadSequence{$seqName});
}
close (FILE);
@superReadNames = keys %superReadSequence;
for (@superReadNames) {
    $superReadName = $_;
    $superReadLength{$superReadName} = length ($superReadSequence{$superReadName}); }

open (FILE, $readPlacementFile);
while ($line = <FILE>) {
    chomp ($line);
    ($readName, $superReadName, $offset, $ori) = split (" ", $line);
    # Keep super-read on the 5' end if needed to use the read as per
    # conversation with Jim on 5/2/11 at 3pm
    $adjustedOffset = $offset;
    $superReadLength = $superReadLength{$superReadName};
    if ($ori eq "F") {
	if ($superReadLength - $offset < $minReadSequenceLength) {
	    $adjustedOffset = $superReadLength - $minReadSequenceLength; }
	if ($adjustedOffset < 0) {
	    $adjustedOffset = 0; } }
    else {
	if ($offset < $minReadSequenceLength) {
	    $adjustedOffset = $minReadSequenceLength; }
	if ($adjustedOffset > $superReadLength) {
	    $adjustedOffset = $superReadLength; } }
    if ($ori eq "F") {
	$outputString = substr ($superReadSequence{$superReadName}, $adjustedOffset); }
    else {
	$outputString = substr ($superReadSequence{$superReadName}, 0, $adjustedOffset);
	$outputString = reverse ($outputString);
	$outputString =~ tr/acgtACGT/tgcaTGCA/; }
    if ($outputString) {
	print ">$readName\n";
	if (length ($outputString) > $maxReadLength) {
	    $outputString = substr ($outputString, 0, $maxReadLength); }
	print "$outputString\n";
    }
}


