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
# Command line is
#  findEquivalentSuperReads.perl [flags] inputFiles
#
# Input files are of the type superReadGroupsForEachReadWhichHasAGroup*
#
# Output by default goes to STDOUT
#
# Flags are as follows
# -equiv-file $file : names a file to report the equivalent super reads;
#                        (otherwise it's not reported)
# -out-file $file : names a file equivalent to the input file, but with
#                   super-reads replaced by their shortest equivalent
# -equivalent-only : only output the equivalencies (-equiv-file must be set)
# -h : help
# 

# NOTE: This also requires kUnitigLengths.txt to exist; each line can be
# of the form 'kUni# length' or just 'length' for the unitigs where the 
# unitigs are numbered from 0.
$kUnitigLengthsFile = "kUnitigLengths.txt";
&processArgs;

# @superReadByReadFiles = ("superReadGroupsForEachReadWhichHasAGroup.final.txt", "superReadGroupsForEachReadWhichHasAGroup.overlapJoinedMates.txt");
&setKUnitigLengths ($kUnitigLengthsFile);

for (@superReadByReadFiles) {
    $superReadByReadFile = $_;
    open (FILE, $superReadByReadFile);
    while ($line = <FILE>) {
	$superRead = <FILE>;
	chomp ($superRead);
	next if ($exists{$superRead});
	$exists{$superRead} = 1;
	@flds = split (" ", $superRead);
	$firstKUni = $flds[0];
	$lastKUni = $flds[-2];
	$index = "$firstKUni $lastKUni";
	$superReadsForIndex{$index} .= "$superRead\n";
	++$count{$index};
    }
    close (FILE);
}

@keys = keys %count;
if ($equivFile) {
    open (EQUIV_FILE, ">$equivFile"); }
for (@keys) {
    $key = $_;
    next unless ($count{$key} > 1);
#    print "$superReadsForIndex{$key}\n";
    &determineIdenticalSuperReads ($superReadsForIndex{$key});
}
if ($equivFile) {
    close (EQUIV_FILE); }
exit (0) if ($equivalentOnly);

if ($outFile) {
    open (OUTFILE, ">$outFile"); }
for (@superReadByReadFiles) {
    $superReadByReadFile = $_;
    open (FILE, $superReadByReadFile);
    while ($line = <FILE>) {
	$superRead = <FILE>;
	chomp ($superRead);
	if (! $superReadEquiv{$superRead}) {
	    if ($outFile) {
		print OUTFILE $line, "$superRead\n"; }
	    else {
		print $line, "$superRead\n"; }
	    next; }
	$newSuperReadStr = $superReadEquiv{$superRead};
	chomp ($line);
	($first, $rest1) = ($line =~ /^(.+ : )(.*\S)\s*$/);
	if ($rest1 eq $superRead) {
	    $newRestOfLine1 = $newSuperReadStr; }
	else {
	    $newRestOfLine1 = &returnReversedSuperReadString ($newSuperReadStr); }
	if ($outFile) {
	    print OUTFILE $first, $newRestOfLine1, "\n", $superReadEquiv{$superRead}, "\n"; }
	else {
	    print $first, $newRestOfLine1, "\n", $superReadEquiv{$superRead}, "\n"; }
    }
    close (FILE);
}
if ($outFile) {
    close (OUTFILE); }

sub determineIdenticalSuperReads
{
    my ($superReadStringsForIndex) = @_;
    my (@superReadStringsForIndex, @indexArr, @outStr);
    local (@superReadLengthsForIndex);
    my ($i, $j, $k, @flds, $toffset, $tstr);
    my (@wasReduced, $indexElementI, $indexElementJ);
    my (@strsToCheck, @tarr, @tempArray1);

    chomp ($superReadStringsForIndex);
    @superReadStringsForIndex = split (/\n/, $superReadStringsForIndex);
    @superReadLengthsForIndex = ();
    for ($i=0; $i<=$#superReadStringsForIndex; $i++) {
	@flds = split (" ", $superReadStringsForIndex[$i]);
	$superReadLengthsForIndex[$i] = $#flds+1;
    }
    @indexArr = (0..$#superReadStringsForIndex);
    @indexArr = sort spcl @indexArr;
    for ($i=0; $i<=$#superReadStringsForIndex; $i++) {
	@flds = split (" ", $superReadStringsForIndex[$i]);
	$outStr[$i] = " $flds[0] 0 $flds[1]";
	$toffset = $kUniLen[$flds[0]];
	for ($j=2; $j<=$#flds; $j+=3) {
	    $toffset -= $flds[$j];
	    $outStr[$i] .= " $flds[$j+1] $toffset $flds[$j+2]";
	    $toffset += $kUniLen[$flds[$j+1]];
	}
	$outStr[$i] .= " ";
    }
    @wasReduced = ();
    for ($i=0; $i<=$#superReadStringsForIndex; $i++) {
	$indexElementI = $indexArr[$i];
	next if ($wasReduced[$indexElementI] =~ /\d/);
	@strsToCheck = ();
	@tempArray1 = split (" ", $outStr[$indexElementI]);
	for ($j=0; $j<=$superReadLengthsForIndex[$indexElementI]; $j+=3) {
	    @tarr = ($tempArray1[$j], $tempArray1[$j+1], $tempArray1[$j+2]);
	    $tstr = " @tarr ";
	    push (@strsToCheck, $tstr); }

	for ($j=$i+1; $j<=$#superReadStringsForIndex; $j++) {
	    $indexElementJ = $indexArr[$j];
	    next if ($wasReduced[$indexElementJ] =~ /\d/);
	    next if ($superReadLengthsForIndex[$indexElementI] == $superReadLengthsForIndex[$indexElementJ]);
	    for ($k=0; $k<=$#strsToCheck; $k++) {
#		print "Checking $strsToCheck[$k] in $outStr[$indexElementJ]\n";
		last unless (index ($outStr[$indexElementJ], $strsToCheck[$k]) >= 0); }
	    $wasReduced[$indexElementJ] = $indexElementI if ($k > $#strsToCheck);
	    
	}
    }
    
#    for ($i=0; $i<=$#superReadStringsForIndex; $i++) {
#	print "$outStr[$indexArr[$i]] ";
#	print "; $superReadLengthsForIndex[$indexArr[$i]] ";
#	print "; $wasReduced[$indexArr[$i]]\n"; }
#    print "\n";
#    print "Results:\n";
    for ($i=0; $i<=$#superReadStringsForIndex; $i++) {
	next unless ($wasReduced[$indexArr[$i]] =~ /\d/);
	
#	print "$outStr[$indexArr[$i]] -> $outStr[$wasReduced[$indexArr[$i]]]\n";
	$superReadEquiv{$superReadStringsForIndex[$indexArr[$i]]} = $superReadStringsForIndex[$wasReduced[$indexArr[$i]]];
	if ($equivFile) {
	    print EQUIV_FILE "$superReadStringsForIndex[$indexArr[$i]] -> $superReadStringsForIndex[$wasReduced[$indexArr[$i]]]\n"; }
#	print "len = $superReadLengthsForIndex[$indexArr[$i]]\n";
    }
#    print "\n";
}

sub spcl
{
    if ($superReadLengthsForIndex[$a] < $superReadLengthsForIndex[$b]) {
	return -1; }
    if ($superReadLengthsForIndex[$a] > $superReadLengthsForIndex[$b]) {
	return 1; }
    return ($a <=> $b);
}

sub setKUnitigLengths
{
    my ($lengthFile) = @_;
    my ($line, @flds, $count);

    open (FILE, $lengthFile);
    $line = <FILE>;
    chomp ($line); @flds = split (" ", $line);
    if ($#flds == 1) {
	$kUniLen[$flds[0]] = $flds[1];
	while ($line = <FILE>) {
	    chomp ($line);
	    @flds = split (" ", $line);
	    $kUniLen[$flds[0]] = $flds[1];
	}
    }
    else {
	$count = 0;
	$kUniLen[$count] = $flds[0];
	while ($line = <FILE>) {
	    chomp ($line);
	    ++$count;
	    $kUniLen[$count] = $flds[0];
	}
    }
    close (FILE);
}

sub processArgs
{
    if ($#ARGV < 0) {
	& reportUsageAndExit; }
    for ($i=0; $i<=$#ARGV; $i++) {
	if ($ARGV[$i] =~ /^[-]/) {
	    if ($ARGV[$i] =~ /^.h/) {
		&reportUsageAndExit; }
	    elsif ($ARGV[$i] eq "-equiv-file") {
		++$i;
		$equivFile = $ARGV[$i];
		next; }
	    elsif ($ARGV[$i] eq "-out-file") {
		++$i;
		$outFile = $ARGV[$i];
		next; }
	    elsif ($ARGV[$i] eq "-equivalent-only") {
		$equivalentOnly = 1;
		next; }
	    else {
		print STDERR "Unrecognized flag '$ARGV[$i]'\n";
		&reportUsageAndExit; }
	}
	else {
	    if (! -e $ARGV[$i]) {
		print STDERR "File '$ARGV[$i]' doesn't exist.\n";
		&reportUsageAndExit; }
	    push (@superReadByReadFiles, $ARGV[$i]);
	}
    }
    if ($equivalentOnly && (! $equivFile)) {
	print STDERR "The -equiv-file flag must be used to specify an output file when specifying -equivalent-only\n";
	&reportUsageAndExit; }
    if ($#superReadByReadFiles < 0) {
	print STDERR "You must give the name of an input file on the command line\n";
	&reportUsageAndExit; }
}
	      
sub reportUsageAndExit
{
    open (FILE, $0);
    print STDERR "\n";
    $line = <FILE>;
    while ($line = <FILE>) {
	last unless ($line =~ /^\#/);
	chomp ($line);
	($line) = ($line =~ /^..(.+)$/);
	print STDERR "$line\n";
    }
    close (FILE);

    exit (1);
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

