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
$readMateFile = $ARGV[0];
$kUnitigVsReadNucmerFile = $ARGV[1];
$pass2file = "mateJoinerPrefixFile.txt";
$pass3file = "mateJoinerSuffixFile.txt";
$pass = 1;
if ($#ARGV >= 2) {
    $pass = $ARGV[2]; }


#readinfo hash contains:
#(readOri,insertSize,insertStdDev,isNeeded,kUniSeq,firstReadsKUniSeq)
%readInfo=();
#kUniInfo hash contains:
#(wasFound,kUni,lengthAdjustment,kUniOri)
%kUniInfo=();

open (FILE, $readMateFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    ($read1, $ori1) = ($flds[0] =~ /^(.+)(.)$/);
    if ($ori1 =~ /\d/) { $read1 .= $ori1; $ori1 = "F"; }
    ($read2, $ori2) = ($flds[1] =~ /^(.+)(.)$/);
    if ($ori2 =~ /\d/) { $read2 .= $ori2; $ori2 = "R"; }
    my @readinfo1_container;
    my @readinfo2_container;
    @readinfo1_container=($ori1,$flds[2],$flds[3],1);
    $readinfo2_container[0]=$ori2;
    $readinfo2_container[3]=2;
    
    push (@read1s, $read1);
    push (@read2s, $read2);
    if ($pass > 1) {
        $readinfo1_container[4]=" @flds[5..$#flds] ";
        $readinfo2_container[4]=" @flds[5..$#flds] ";
	if ($pass > 2) {
            $readinfo1_container[5]=$flds[4];}
    }
    $readInfo{$read1}=\@readinfo1_container;
    $readInfo{$read2}=\@readinfo2_container;
}
close (FILE);

if ($pass == 2) {
    open (SPCL_OUTFILE, ">$pass2file"); }
elsif ($pass == 3) {
    open (SPCL_OUTFILE, ">$pass3file"); }
# For the moment we will assume that all read1s are 'F' and read2s are 'R'
open (FILE, $kUnitigVsReadNucmerFile);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    next unless (${$readInfo{$flds[-1]}}[3]);
    $read = $flds[-1];
    next if (${$kUniInfo{$read}}[0] && (${$readInfo{$read}}[3] >= $pass));
    $kUni = $flds[-2];
    if ($read ne $readHold) {
	# Do the prior output for $readHold
	if ($pass > 1) {
	    &analyzeNucmerLines ($readHold); }
	# Now continue
	$readHold = $read;
	@heldNucmerLines = ();
    }
    my @kUniInfo_container;
    push (@heldNucmerLines, $line);
    if (($pass > 1) && ${$kUniInfo{$read}}[0]) {
	if ($pass - ${$readInfo{$read}}[3] == 1) {
	    next unless (index (${$readInfo{$read}}[4], " $kUni ") >= 0); }
	if ($pass - ${$readInfo{$read}}[3] == 2) {
	    next unless ($kUni eq ${$readInfo{$read}}[5]); }
    }
    $kUniInfo_container[1] = $kUni;
    if ($flds[6] < $flds[7]) {
	$kUniInfo_container[2] = $flds[4] - $flds[6]; }
    else {
	$kUniInfo_container[2] = ($flds[4]-1) - $flds[6]; }
    if ($flds[6] < $flds[7]) {
	$kUniInfo_container[3] = ${$readInfo{$read}}[0]; }
    else {
	if (${$readInfo{$read}}[3] == 1) { # The first read of a pair
	    $kUniInfo_container[3] = "R"; }
	else { # The second read of a pair
	    $kUniInfo_container[3] = "F"; }
    }
    $kUniInfo_container[0]=1;
    $kUniInfo{$read}=\@kUniInfo_container;
}
if ($pass > 1) {
    &analyzeNucmerLines ($readHold); }
close (FILE);
close (SPCL_OUTFILE);

for ($i=0; $i<=$#read1s; $i++) {
    $read1 = $read1s[$i];
    $read2 = $read2s[$i];
    $kUni1 = ${$kUniInfo{$read1}}[1];
    $kUni2 = ${$kUniInfo{$read2}}[1];
    $insertSize = ${$readInfo{$read1}}[1];
    $insertStdDev = ${$readInfo{$read1}}[2];
    $ori1 = ${$kUniInfo{$read1}}[3];
    $ori2 = ${$kUniInfo{$read2}}[3];
    $insertSize += (${$kUniInfo{$read1}}[2] + ${$kUniInfo{$read2}}[2]);
    print "${kUni1}$ori1 ${kUni2}$ori2 $insertSize $insertStdDev\n";
}

sub analyzeNucmerLines
{
    my ($localRead) = @_;
    my (@oris, @minOffsets, @maxOffsets, $line, @flds, $readOffset1, $readOffset2);
    my ($good, $i, $index, @overlaps, @kUnitigNums, $tOri);

    return unless ($localRead =~ /\S/);
    @kUnitigNums = @oris = @minOffsets = @maxOffsets = @overlaps = ();
    for (@heldNucmerLines) {
	$line = $_;
	@flds = split (" ", $line);
	$readOffset1 = $flds[6];
	$readOffset2 = $flds[7];
	$kUnitigNum = $flds[-2];
	push (@kUnitigNums, $kUnitigNum);
	if ($readOffset1 < $readOffset2) {
	    push (@oris, "F");
	    push (@minOffsets, $readOffset1);
	    push (@maxOffsets, $readOffset2); }
	else {
	    push (@oris, "R");
	    push (@minOffsets, $readOffset2);
	    push (@maxOffsets, $readOffset1); }
    }
    $good = 1;
    for ($i=1; $i<=$#minOffsets; $i++) {
	$overlaps[$i-1] = ($maxOffsets[$i-1]+1) - $minOffsets[$i];
	# The following test should never apply for what we're doing
	if ($overlaps[$i-1] < 0) {
	    $good = 0; }
    }
    print SPCL_OUTFILE $localRead;
    if ($pass == 2) {
	for ($i=0; 1; $i++) {
	    last if ($kUnitigNums[$i] eq ${$kUniInfo{$localRead}}[1]);
	    print SPCL_OUTFILE " $kUnitigNums[$i] $oris[$i] $overlaps[$i]"; }
    }
    elsif ($pass == 3) {
	$isStarted = 0;
	for ($i=$#kUnitigNums; $i>=0; $i--) {
	    if ($kUnitigNums[$i] eq ${$kUniInfo{$localRead}}[1]) {
		$isStarted = 1;
		next; }
	    if ($isStarted) {
		if ($oris[$i] eq "F") { $tOri = "R"; }
		else { $tOri = "F"; }
		print SPCL_OUTFILE " $overlaps[$i] $kUnitigNums[$i] $tOri";
	    }
	}
    }
    print SPCL_OUTFILE "\n";
}

