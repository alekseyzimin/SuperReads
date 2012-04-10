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
$kmerLen = 31;
$fileNum = 0;
for ($i=0; $i<=$#ARGV; $i++) {
    if ($ARGV[$i] eq "-l") {
	++$i;
	$kmerLen = $ARGV[$i];
	next; }
    if ($fileNum == 0) { $inputFileWithReadPairsToBeJoined = $ARGV[$i]; }
    elsif ($fileNum == 1) { $inputFileWithKUnitigPairsToBeJoined = $ARGV[$i]; }
    elsif ($fileNum == 2) { $matePairJoinerResultsFile = $ARGV[$i]; }
    elsif ($fileNum == 3) { $outputJoinFile = $ARGV[$i]; }
    elsif ($fileNum == 4) { $outputFileWithKUnitigsToKeep = $ARGV[$i]; }
    else {
	print STDERR "Exec $0 was called with too many args. Bye!\n";
	exit (1); }
    ++$fileNum;
}
$kmerLenMinus1 = $kmerLen - 1;

open (FILE1, $inputFileWithReadPairsToBeJoined);
open (FILE2, $inputFileWithKUnitigPairsToBeJoined);
open (INPUT_RESULTS_FILE, $matePairJoinerResultsFile);
open (OUTPUT_JOIN_FILE, ">$outputJoinFile");
open (OUTPUT_FILE_WITH_KUNIS_TO_KEEP, ">$outputFileWithKUnitigsToKeep");

$resultsLine = <INPUT_RESULTS_FILE>;
while ($tline = <INPUT_RESULTS_FILE>) {
    chomp ($tline);
    if ($tline !~ /^MATE/) {
	$resultsLine .= "$tline\n";
	next; }
    $readPairLine = <FILE1>; chomp ($readPairLine);
    @flds = split (" ", $readPairLine); $readPairLine = "@flds[0..3]";
    $kUnitigPairLine = <FILE2>; chomp ($kUnitigPairLine);
    # Now do the analysis for the resultsLine
    &reportResults ($resultsLine);
    # Now start over
    $resultsLine = "$tline\n";
}
$readPairLine = <FILE1>; chomp ($readPairLine);
@flds = split (" ", $readPairLine); $readPairLine = "@flds[0..3]";
$kUnitigPairLine = <FILE2>; chomp ($kUnitigPairLine);
# Now do the analysis for the resultsLine
&reportResults ($resultsLine);

close (FILE1);
close (FILE2);
close (INPUT_RESULTS_FILE);
close (OUTPUT_JOIN_FILE);
close (OUTPUT_FILE_WITH_KUNIS_TO_KEEP);

sub reportResults
{
    my ($resultsLines) = @_;
    my (@resultsLines, @kUnis, @oris, @numOvlsIn, @numOvlsOut);
    my ($resultsLine, @flds, $kUni, $ori, $numOvlsIn, $numOvlsOut);
    my ($firstBadLine, $lastBadLine, $i, $numPossibleLengths);

    chomp ($resultsLines);
    @resultsLines = split (/\n/, $resultsLines);
    return if ($#resultsLines < 3); # No connection possible
    @kUnis = @oris = @numOvlsIn = @numOvlsOut = ();
    $numPossibleLengths = 0;
    for (@resultsLines) {
	$resultsLine = $_;
	if ($resultsLine =~ /^\d+\s*$/) {
	    ++ $numPossibleLengths; }
	next unless ($resultsLine =~ /^uni/);
	@flds = split (" ", $resultsLine);
	($kUni) = ($flds[2] =~ /^(\d+)/);
	($ori) = ($flds[8] =~ /^(.)/);
	($numOvlsIn) = ($flds[17] =~ /^(\d+)/);
	($numOvlsOut) = ($flds[20] =~ /^(\d+)/);
	push (@kUnis, $kUni);
	push (@oris, $ori);
	push (@numOvlsIn, $numOvlsIn);
	push (@numOvlsOut, $numOvlsOut);
    }
    $firstBadLine = -1; $lastBadLine = -1;
    if ($numPossibleLengths > 1) {
	$firstBadLine = 0;
	$lastBadLine = $#numOvlsOut+1;
	goto outputABadPair; }
    for ($i=0; $i<=$#numOvlsOut; $i++) {
	if (($numOvlsOut[$i] > 1) || ($numOvlsIn[$i] > 1)) {
	    $firstBadLine = $i;
	    last; } }
    if ($firstBadLine < 0) { goto itsGood; }
    for ($i=$#numOvlsIn; $i>=0; $i--) {
	if (($numOvlsOut[$i] > 1) || ($numOvlsIn[$i] > 1)) {
	    $lastBadLine = $i;
	    last; } }
  outputABadPair:
    if ($lastBadLine <= $firstBadLine+1) { # Nothing acceptable
	print "NOTHING ACCEPTABLE (joins of too many possible lengths): $readPairLine ; $kUnitigPairLine\n";
	return; }
    print OUTPUT_FILE_WITH_KUNIS_TO_KEEP $readPairLine;
    print OUTPUT_FILE_WITH_KUNIS_TO_KEEP " $kUnis[0]";
    for ($i=$firstBadLine+1; $i<$lastBadLine; $i++) {
	print OUTPUT_FILE_WITH_KUNIS_TO_KEEP " $kUnis[$i]"; }
    print OUTPUT_FILE_WITH_KUNIS_TO_KEEP "\n";
    return;

    itsGood:
    print OUTPUT_JOIN_FILE $readPairLine;
    for ($i=0; $i<=$#kUnis; $i++) {
	if ($i > 0) { print OUTPUT_JOIN_FILE " $kmerLenMinus1"; }
	print OUTPUT_JOIN_FILE " $kUnis[$i] $oris[$i]"; }
    print OUTPUT_JOIN_FILE "\n";
}


