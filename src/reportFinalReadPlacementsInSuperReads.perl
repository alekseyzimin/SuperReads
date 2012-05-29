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
$kUnitigsLengthsFile = $ARGV[0];
$priorReadPlacementFile = $ARGV[1]; # e.g. allData5.5.v4/readPlacementsInSuperReads.preMateMerge.read.superRead.offset.ori.txt
$newSuperReadPerReadFile = $ARGV[2];

open (FILE, $kUnitigsLengthsFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $kUnitigLengths[$flds[0]] = $flds[1]; }
close (FILE);

open (FILE, $priorReadPlacementFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $readName = $flds[0];
    $origSuperRead{$readName} = $flds[1];
    $origOffset{$readName} = $flds[2];
    $origOri{$readName} = $flds[3]; }
close (FILE);

open (FILE, $newSuperReadPerReadFile);
$line = <FILE>;
close (FILE);
@flds = split (" ", $line);
if ($flds[0] eq "readNum") {
    $readNameIndex = 5; }
else {
    $readNameIndex = 2; }

open (FILE, $newSuperReadPerReadFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $readName = $flds[$readNameIndex];
    $superReadWithSpaces = <FILE>;
    chomp ($superReadWithSpaces);
    $priorSuperReadWithSpaces = $origSuperRead{$readName};
    $priorSuperReadWithSpaces =~ s/_/ /g;
    $priorSuperReadWithSpaces =~ s/F/ F/g;
    $priorSuperReadWithSpaces =~ s/R/ R/g;
    $searchStringF = " $priorSuperReadWithSpaces ";
    $biggerString = " $superReadWithSpaces ";
    $weWantTheOutputToTest = 0;
    $forwardIndex = index ($biggerString, $searchStringF);
    if ($forwardIndex < 0) { # It didn't match in the forward direction
	goto tryTheReverse; }
# The next section will be included once we are sure the rest is working
#    if ($forwardIndex == 0) {
#	$newOffset = $origOffset{$readName};
#	$newOri = $origOri{$readName};
#	goto outputTheRecord; }
    if ($forwardIndex >= 0) {
	$offsetAdjustment = &getOffsetAdjustment ($biggerString, $forwardIndex, "F");
#	if ($forwardIndex > 0) {
#	    print "readName = $readName ; searchStringF = \"$searchStringF\" ; biggerString = \"$biggerString\" ; offsetAdjustment = $offsetAdjustment\n"; }
	$newOffset = $offsetAdjustment + $origOffset{$readName};
	$newOri = $origOri{$readName};
	goto outputTheRecord;
    }
  tryTheReverse:
    $searchStringR = &returnReversedSuperReadString ($priorSuperReadWithSpaces);
    $searchStringR = " $searchStringR ";
    $reverseIndex = index ($biggerString, $searchStringR);
    if ($reverseIndex < 0) {
#	print STDERR "Read name = $readName\n";
#	print STDERR "Disaster: neither forward string \"$searchStringF\" nor reverse string \"$searchStringR\" found in new string \"$biggerString\". Quittinging.\n";
#	exit(1);
	# We now allow this since we are now forcing a super-read "merge" for mates with a common
	# k-unitig at an end where the orientation is correct
	$newOffset = -1000;
	$newOri = "F";
	goto outputTheRecord;
    }
    $offsetAdjustment = &getOffsetAdjustment ($biggerString, $reverseIndex + length ($searchStringR), "R");
#    if ($readName eq "SRR081522.593482/2") {
#    if (1) {
#	print "readName = $readName ; searchStringR = \"$searchStringR\" ; biggerString = \"$biggerString\" ; offsetAdjustment = $offsetAdjustment\n"; }
    $newOffset = $offsetAdjustment - $origOffset{$readName};
    if ($origOri{$readName} eq "F") {
	$newOri = "R"; }
    else {
	$newOri = "F"; }
    
  outputTheRecord:
    $superReadName = &makeSuperReadNameFromSuperReadWithSpaces ($superReadWithSpaces);
    if ($weWantTheOutputToTest) {
	print "testRec "; }
    print "$readName $superReadName $newOffset $newOri\n";
}

sub returnReversedSuperReadString
{
    my ($origSuperReadStr) = @_;
    my (@flds, $i, $modifiedSuperReadStr, $loopVar);

    @flds = split (" ", $origSuperReadStr);
    $i=$#flds-1;
    $modifiedSuperReadStr = $flds[$i] . " ";
    if ($flds[$i+1] eq "F") {
        $modifiedSuperReadStr .= "R"; }
    else {
        $modifiedSuperReadStr .= "F"; }
    for ($loopVar = $i-3; $loopVar >= 0; $loopVar -= 3) {
        $modifiedSuperReadStr .= (" " . $flds[$loopVar+2] . " " . $flds[$loopVar] . " ");
        if ($flds[$loopVar+1] eq "F")  {
            $modifiedSuperReadStr .= "R"; }
        else {
            $modifiedSuperReadStr .= "F"; } }
    return ($modifiedSuperReadStr);
}

sub getOffsetAdjustment
{
    my ($inputString, $numCharsUsed, $oriComparedToOrig) = @_;
    my ($tstr, @tflds, $startOffset, $i);

    if ($numCharsUsed <= 0) { # Only = 0 and only if ori is same as before
	return (0); }
#    if ($oriComparedToOrig ne "F") {
#	$weWantTheOutputToTest = 1; }
    $tstr = substr ($inputString, 0, $numCharsUsed);
    @tflds = split (" ", $tstr);
#    if ($readName eq "SRR081522.593482/2") {
#    if (1) {
#	print "readName = $readName ; start unitig = $tflds[0]\n"; }
    $startOffset = $kUnitigLengths[$tflds[0]];
    for ($i=2; $i<=$#tflds; $i+=3) {
#	if ($readName eq "SRR081522.593482/2") {
#	if (1) {
#	    print "readName = $readName ; startOffset = $startOffset ; ";
#	    if ($i+1 <= $#tflds) {
#		print "unitig# = $tflds[$i+1] ; "; } }
	$startOffset -= $tflds[$i];
	if ($i+1 <= $#tflds) {
	    $startOffset += $kUnitigLengths[$tflds[$i+1]]; }
#	if ($readName eq "SRR081522.593482/2") {
#	if (1) {
#	    print "endOffset = $startOffset\n"; }
	
    }
    return ($startOffset);
}

sub makeSuperReadNameFromSuperReadWithSpaces
{
    my ($inputStr) = @_;
    my (@flds, $outStr, $i);
    
    @flds = split (" ", $inputStr);
    $outStr = $flds[0] . $flds[1];
    for ($i=2; $i<$#flds; $i+=3) {
	$outStr .= ("_" . $flds[$i] . "_" . $flds[$i+1] . $flds[$i+2]);
    }

    return ($outStr);
}

