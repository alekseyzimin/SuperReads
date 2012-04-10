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
use File::Basename;

&processArgs;
$statsFile = &setStatsFilename;

open (FILE, $mateInfoFile);
$insertNum = 0;
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $strF = "$insertNum F";
    $strR = "$insertNum R";
    $readToInsertInfo{$flds[0]} = $strF;
    $readToInsertInfo{$flds[1]} = $strR;
    $insertInfoToRead{$strF} = $flds[0];
    $insertInfoToRead{$strR} = $flds[1];
    ++$insertNum; }
close (FILE);

if ($repetitiveKUnitigsFile) {
    open (FILE, $repetitiveKUnitigsFile);
    while ($line = <FILE>) {
	chomp ($line);
	@flds = split (" ", $line);
	next if ($flds[1] >= $minSepNonRptKUni);
	$isRepetitive{$flds[0]} = 1; }
    close (FILE); }

if ($chimericReadsFile) {
    open (FILE, $chimericReadsFile);
    while ($chimericRead = <FILE>) {
	chomp ($chimericRead);
	$isChimeric{$chimericRead} = 1; }
    close (FILE); }

open (FILE, $kUnitigLengthsFile);
while ($line = <FILE>) {
    chomp ($line);
    ($kUnitigNum, $kUnitigLength) = split (" ", $line);
    $kUnitigLengths[$kUnitigNum] = $kUnitigLength; }
close (FILE);

open (FILE, $superReadsFile);
while ($line = <FILE>) {
    $line2 = <FILE>;
    chomp ($line);
    ($readName, $origSuperReadStr) = ($line =~ /readName =\s(\S+)\s*:\s*(\S.+\S)\s*$/);
    if ($isChimeric{$readName}) {
	$wasOutput{$readName} = 1;
	next; }
    push (@readNames, $readName);
    $readInsertInfo = $readToInsertInfo{$readName};
    if (! $readInsertInfo) {
	print "$line\n", $line2;
	$wasOutput{$readName} = 1;
	next; }
    ($readInsertNum, $readInsertEnd) = split (" ", $readInsertInfo);
    ++$count{$readInsertNum};
    if ($count{$readInsertNum} == 1) {
	chomp ($line2);
	$lines{$readName} = "$line\n" . $line2;
	next; }
    # If we get here it's the second read for the insert (not checking if more)
    $curReadInsertInfo = $readToInsertInfo{$readName};
    chomp ($line2);
    $otherReadInsertInfo = $curReadInsertInfo;
    $otherReadInsertInfo =~ tr/FR/RF/;
    $otherRead = $insertInfoToRead{$otherReadInsertInfo};
    ($otherReadLine1, $otherReadLine2) = split (/\n/, $lines{$otherRead});
    if ($curReadInsertInfo =~ /F$/) {
	$forwardReadLine1 = $line;
	$forwardReadLine2 = $line2;
	$reverseReadLine1 = $otherReadLine1;
	$reverseReadLine2 = $otherReadLine2; }
    else {
	$reverseReadLine1 = $line;
	$reverseReadLine2 = $line2;
	$forwardReadLine1 = $otherReadLine1;
	$forwardReadLine2 = $otherReadLine2; }
    ($forwardReadName, $origSuperReadStr) = ($forwardReadLine1 =~ /readName =\s(\S+)\s*:\s*(\S.+\S)\s*$/);
    $forwardKUnitigString = $origSuperReadStr;
    ($reverseReadName, $origSuperReadStr) = ($reverseReadLine1 =~ /readName =\s(\S+)\s*:\s*(\S.+\S)\s*$/);
    $reverseKUnitigString = &returnReversedSuperReadString ($origSuperReadStr);

    $forwardKUnitigString = " $forwardKUnitigString ";
    $reverseKUnitigString = " $reverseKUnitigString ";
    if ($forwardKUnitigString eq $reverseKUnitigString) {
	print "$forwardReadLine1\n$forwardReadLine2\n";
	print "$reverseReadLine1\n$reverseReadLine2\n";
	$wasOutput{$forwardReadName} = $wasOutput{$reverseReadName} = 1;
	++$identicalMateSuperReads;
#	print "Both mates have the same super-read\n";
	next; }
    $lengthHold = length ($forwardKUnitigString);
    $length = length ($reverseKUnitigString);
    if ($lengthHold > $length) {
	$index = index ($forwardKUnitigString, $reverseKUnitigString);
	if ($index >= 0) {
#	    print "$forwardReadLine1\n$forwardReadLine2\n"; # DEBUGGING
#	    print "$reverseReadLine1\n$reverseReadLine2\n";  #DEBUGGING
	    @kUnitigsInCommon = &getKUnitigsFromSuperReadString ($reverseKUnitigString);
	    $good = &reportIfSomeKUnitigInArrayIsNonRepetitive (@kUnitigsInCommon);
#	    print "One super-read is a subset of the other (at 1). good = $good\n";
	    print "$forwardReadLine1\n$forwardReadLine2\n";
	    if (! $good) {
		print "$reverseReadLine1\n$reverseReadLine2\n";
	    }
	    else {
		++$numContainedSuperReadsAt1;
		($linePrefix) = ($reverseReadLine1 =~ /^(.*readName =\s\S+\s*:\s*)\S/);
		print $linePrefix, $forwardKUnitigString,"\n";
		print "$forwardReadLine2\n";
	    }
	    $wasOutput{$forwardReadName} = $wasOutput{$reverseReadName} = 1;
	    next; } }
    elsif ($lengthHold < $length) {
	$index = index ($reverseKUnitigString, $forwardKUnitigString);
	if ($index >= 0) {
#	    print "$forwardReadLine1\n$forwardReadLine2\n"; # DEBUGGING
#	    print "$reverseReadLine1\n$reverseReadLine2\n";  #DEBUGGING
	    @kUnitigsInCommon = &getKUnitigsFromSuperReadString ($forwardKUnitigString);
	    $good = &reportIfSomeKUnitigInArrayIsNonRepetitive (@kUnitigsInCommon);
#	    print "One super-read is a subset of the other (at 2). good = $good\n";
	    if (! $good) {
		print "$forwardReadLine1\n$forwardReadLine2\n";
	    }
	    else {
		++$numContainedSuperReadsAt2;
		($linePrefix) = ($forwardReadLine1 =~ /^(.*readName =\s\S+\s*:\s*)\S/);
		print $linePrefix, $reverseKUnitigString,"\n";
		print "$reverseReadLine2\n";
	    }
	    print "$reverseReadLine1\n$reverseReadLine2\n";
	    $wasOutput{$forwardReadName} = $wasOutput{$reverseReadName} = 1;
	    next; } }
    # The only time we should get here is if neither is a substring of the other
    # if $endID eq "a" then we saw the reverse read first;
    # otherwise we saw the forward read first
    # The following returns the string "FALSE" if it fails
    $returnSuperRead = &calculateInsertSuperRead;
#    print "$forwardReadLine1\n$forwardReadLine2\n"; DEBUGGING
#    print "$reverseReadLine1\n$reverseReadLine2\n"; DEBUGGING
#    print "We haven't covered this case yet!\n"; DEBUGGING
#    print "$forwardKUnitigString != $reverseKUnitigString\n"; DEBUGGING
#    print "We report a new super-read of \"$returnSuperRead\"\n"; DEBUGGING
    if ($returnSuperRead eq "FAIL") {
#	print "FAIL $forwardReadLine1\nFAIL $forwardReadLine2\n";
#	print "FAIL $reverseReadLine1\nFAIL $reverseReadLine2";
	@forwardFields = split (" ", $forwardKUnitigString);
	@reverseFields = split (" ", $reverseKUnitigString);
	$canForceMerge = 0;
	$tempForwardKUnitig = $forwardFields[0]; $tempForwardOri = $forwardFields[1];
	$tempReverseKUnitig = $reverseFields[0]; $tempReverseOri = $reverseFields[1];
	if (($tempForwardKUnitig == $tempReverseKUnitig) && ($tempForwardOri eq $tempReverseOri) && (! $isRepetitive{$tempForwardKUnitig})) {
#	    print "merge k-uniLen = $kUnitigLengths[$tempForwardKUnitig] ; ";
	    $canForceMerge = 1; }
	$tempForwardKUnitig = $forwardFields[-2]; $tempForwardOri = $forwardFields[-1];
	$tempReverseKUnitig = $reverseFields[-2]; $tempReverseOri = $reverseFields[-1];
	if (($tempForwardKUnitig == $tempReverseKUnitig) && ($tempForwardOri eq $tempReverseOri) && (! $isRepetitive{$tempForwardKUnitig})) {
#	    print " merge k-uniLen = $kUnitigLengths[$tempForwardKUnitig] ; ";
	    $canForceMerge = 1; }
	$forwardLen = &getSuperReadLengthFromSuperReadString ($forwardKUnitigString);
	$reverseLen = &getSuperReadLengthFromSuperReadString ($reverseKUnitigString);
#	print " ; lengths = ($forwardLen, $reverseLen)";
	$changeLine = 0;
	if ($canForceMerge) {
	    ++$forcedMerges;
	    if ($forwardLen >= $reverseLen) {
		$changeLine = 2; }
	    else {
		$changeLine = 1; }
	    ++$diffAndFailToMerge; }
	if ($changeLine != 1) {
	    print "$forwardReadLine1\n$forwardReadLine2\n"; }
	else {
	    ($linePrefix) = ($forwardReadLine1 =~ /^(.*readName =\s\S+\s*:\s*)\S/);
	    print $linePrefix, $reverseKUnitigString,"\n";
	    print "$reverseReadLine2\n";
	}
	if ($changeLine != 2) {
	    print "$reverseReadLine1\n$reverseReadLine2\n"; }
	else {
	    ($linePrefix) = ($reverseReadLine1 =~ /^(.*readName =\s\S+\s*:\s*)\S/);
	    print $linePrefix, $forwardKUnitigString,"\n";
	    print "$forwardReadLine2\n";
	}
#	if ($canForceMerge) {
#	    print " Can force merge"; }
#	print "\n";
    }
    else {
	$returnSuperRead = &getMinOfSuperReadAndReverseComplement ($returnSuperRead);
	++$newlyMerged;
	($linePrefix) = ($forwardReadLine1 =~ /^(.*readName =\s\S+\s*:\s*)\S/);
	print $linePrefix, $returnSuperRead,"\n";
	print $returnSuperRead,"\n";
	($linePrefix) = ($reverseReadLine1 =~ /^(.*readName =\s\S+\s*:\s*)\S/);
	print $linePrefix, $returnSuperRead,"\n";
	print $returnSuperRead,"\n"; }
    $wasOutput{$forwardReadName} = $wasOutput{$reverseReadName} = 1;
}

for (@readNames) {
    $readName = $_;
    next if ($wasOutput{$readName});
    print $lines{$readName}, "\n";
}

open (FILE, ">$statsFile");
print FILE "num identical super-reads: $identicalMateSuperReads\n";
print FILE "num contained at 1: $numContainedSuperReadsAt1\n";
print FILE "num contained at 2: $numContainedSuperReadsAt2\n";
print FILE "diff and merged: $newlyMerged\n";
print FILE "forced (dovetail) mate merges: $forcedMerges\n";
print FILE "diff but failed to merge: $diffAndFailToMerge\n";
close (FILE);

sub calculateInsertSuperRead
{
    my ($localForwardKUnitigString, $localReverseKUnitigString);
    my (@forwardFields, @reverseFields);
    my ($tstr, $tstr2, $offsetOfLastKUnitig, @tflds2);
    my ($numFieldsToAdd, $numFieldsTotal);
    my ($firstMatchingFieldInForward, $lastMatchingFieldInReverse);
    my ($i, $j, @newSuperReadFields, $newSuperReadString);
    my (@localListOfSharedKUnitigs, $good);

    $localForwardKUnitigString = " " . $forwardKUnitigString . " ";
    $localReverseKUnitigString = " " . $reverseKUnitigString . " ";
    @forwardFields = split (" ", $forwardKUnitigString);
    @reverseFields = split (" ", $reverseKUnitigString);
    $tstr = " " . $forwardFields[-2] . " " . $forwardFields[-1] . " ";
    $offsetOfLastKUnitig = index ($localReverseKUnitigString, $tstr);
    if ($offsetOfLastKUnitig < 0) {
	return ("FAIL"); }
    $tstr2 = substr ($localReverseKUnitigString, 0, $offsetOfLastKUnitig);
    @tflds2 = split (" ", $tstr2);
    $numFieldsToAdd = $#tflds2+1;
    $numFieldsTotal = $numFieldsToAdd+2;
    $firstMatchingFieldInForward = $#forwardFields - ($numFieldsTotal-1);
    $lastMatchingFieldInReverse = $numFieldsTotal-1;
    @localListOfSharedKUnitigs = ();
    for ($i=$firstMatchingFieldInForward, $j=0; $i<=$#forwardFields; $i++, $j++) {
	if ($j% 3 == 0) {
	    push (@localListOfSharedKUnitigs, $reverseFields[$j]); }
	if ($forwardFields[$i] ne $reverseFields[$j]) {
	    return ("FAIL"); } }
    # If we got here all the fields matched
    $good = &reportIfSomeKUnitigInArrayIsNonRepetitive (@localListOfSharedKUnitigs);
    if (! $good) {
#	print "Failed because the matching k-unitigs were all repetitive\n";
	return ("FAIL"); }
    @newSuperReadFields = ();
    for ($i=0; $i<$firstMatchingFieldInForward; $i++) {
	push (@newSuperReadFields, $forwardFields[$i]); }
    push (@newSuperReadFields, @reverseFields);
    $newSuperReadString = "@newSuperReadFields";
    return ($newSuperReadString);
}

sub processArgs
{
    my ($fileNum, $i);

    $minSepNonRptKUni = 1000000;
    $fileNum = 0;
    for ($i=0; $i<=$#ARGV; $i++) {
	if ($ARGV[$i] eq "-min-sep-non-rpt-k-uni") {
	    ++$i;
	    $minSepNonRptKUni = $ARGV[$i];
	    next; }
	if ($fileNum == 0) {
	    $superReadsFile = $ARGV[$i]; }
	elsif ($fileNum == 1) {
	    $mateInfoFile = $ARGV[$i]; }
	elsif ($fileNum == 2) {
	    $chimericReadsFile = $ARGV[$i]; }
	elsif ($fileNum == 3) {
	    $repetitiveKUnitigsFile = $ARGV[$i]; }
	elsif ($fileNum == 4) {
	    $kUnitigLengthsFile = $ARGV[$i]; }
	++$fileNum; }
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

sub getKUnitigsFromSuperReadString
{
    my ($superReadString) = @_;
    my (@flds, @flds2, $i);

    @flds = split (" ", $superReadString);
    @flds2 = ();
    for ($i=0; $i<=$#flds; $i++) {
	push (@flds2, $flds[$i]); }
    return (@flds2);
}

sub reportIfSomeKUnitigInArrayIsNonRepetitive
{
    my (@flds) = @_;
    my ($food, $fld);

    $good = 0;
    for (@flds) {
	$fld = $_;
	next if ($isRepetitive{$fld});
	$good = 1;
	last; }
    return ($good);
}

sub getMinOfSuperReadAndReverseComplement
{
    my ($superRead) = @_;
    my ($superReadReversed, @flds1, @flds2, $i);

    $superReadReversed = &returnReversedSuperReadString ($superRead);
    @flds1 = split (" ", $superRead);
    @flds2 = split (" ", $superReadReversed);
    for ($i=0; $i<=$#flds1; $i+=3) {
	if ($flds1[$i] < $flds2[$i]) {
	    return ($superRead); }
	if ($flds1[$i] > $flds2[$i]) {
	    return ($superReadReversed); }
	if ($flds1[$i+1] lt $flds2[$i+1]) {
	    return ($superRead); }
	if ($flds1[$i+1] gt $flds2[$i+1]) {
	    return ($superReadReversed); } }
    return ($superRead);
}

sub setStatsFilename
{
    my ($statsFileOutDir, $execPrefix, $statsFile);

    $statsFileOutDir = dirname ($superReadsFile);
    $execPrefix = basename ($0);
    ($execPrefix) = ($execPrefix =~ /^([^\.]+)\./);
    $statsFile = $statsFileOutDir . "/" . $execPrefix . ".stats.txt";
    return ($statsFile);
}

sub getSuperReadLengthFromSuperReadString
{
    my ($localSuperReadString) = @_;
    my ($localLength, @flds, $i);

    $localLength = 0;
    @flds = split (" ", $localSuperReadString);
    for ($i=0; $i<=$#flds; $i+=3) {
	$localLength += $kUnitigLengths[$flds[$i]]; }
    for ($i=2; $i<=$#flds; $i+=3) {
	$localLength -= $flds[$i]; }
    return ($localLength);
}

