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
while ($line = <STDIN>) {
    chomp ($line);
    ($read1, $read2, $rest) = ($line =~ /^\s*(\S+)\S\s+(\S+)\S\s+\d+\s+\d+\s+(\S.*\S)\s*$/);
    @flds = split (" ", $rest);
    @flds2 = ();
    for ($i=0; $i<=$#flds; $i+=3) {
	$flds2[$#flds-1-$i] = $flds[$i]; }
    for ($i=1; $i<=$#flds; $i+=3) {
	if ($flds[$i] eq "F") {
	    $flds2[$#flds+1-$i] = "R"; }
	else {
	    $flds2[$#flds+1-$i] = "F"; } }
    for ($i=2; $i<=$#flds; $i+=3) {
	$flds2[$#flds-$i] = $flds[$i]; }
#    print "@flds\n";
#    print "@flds2\n";
    $lowerValuedRecord = 0;
    for ($i=0; $i<=$#flds; $i++) {
	if ($i % 3 != 1) {
	    if ($flds[$i] < $flds2[$i]) {
		$lowerValuedRecord = 1; }
	    elsif ($flds[$i] > $flds2[$i]) {
		$lowerValuedRecord = 2; } }
	else {
	    if ($flds[$i] lt $flds2[$i]) {
		$lowerValuedRecord = 1; }
	    elsif ($flds[$i] gt $flds2[$i]) {
		$lowerValuedRecord = 2; } }
	last if ($lowerValuedRecord > 0); }
#    if ($lowerValuedRecord != 0) {
#	print "Chosen record: "; }
    if ($lowerValuedRecord == 1) {
	$superReadName = "@flds"; }
#    elsif ($lowerValuedRecord == 2) {
    else {
	$superReadName = "@flds2"; }
    print "readName = $read1 : @flds\n$superReadName\n";
    print "readName = $read2 : @flds2\n$superReadName\n";
}

