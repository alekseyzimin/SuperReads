#!/usr/bin/env perl
# We assume the following is one line of sequence per faux read
# $CeleraTerminatorDirectory is of the form */9-terminator
# Input files:
# 1) contig_end_pairs.fa
# 2) $workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt
# 3) $workingDirectory/superReadSequences.fasta
# 5) $CeleraTerminatorDirectory/genome.posmap.ctgscf
# 6) $CeleraTerminatorDirectory/genome.ctg.fasta
$contigEndPairsFile = "contig_end_pairs.fa";
&processArgs;
if (! $CeleraTerminatorDirectory) {
    print STDERR "No Celera terminator directory was specified. Bye!\n";
    exit (1);
}

# The following loop sets readLengths, contig, and contigOri for each faux read (contig ends)
open (FILE, $contigEndPairsFile);
while ($line = <FILE>) {
    chomp ($line);
    if ($line !~ /^>/) {
	$readLengths{$readName} = length ($line);
	next; }
    @flds = split (" ", $line);
    ($readName) = ($flds[0] =~ /^.(.+)$/);
    $contig{$readName} = $flds[1];
    $contigOriFromReadName{$readName} = $flds[2];
}
close (FILE);

# The following loads the super-read sequences corresponding to each super-read name
open (FILE, "$workingDirectory/superReadSequences.fasta");
while ($superRead = <FILE>) {
    chomp ($superRead);
    $superRead =~ s/^.//;
    $superReadSequence = <FILE>;
    chomp ($superReadSequence);
    $superReadSequence{$superRead} = $superReadSequence;
}
close (FILE);

# For the first of the two faux reads corresponding to each gap, this extracts the sequence
# to replace that gap, including the two faux reads around it (in %joiningSuperRead)
open (FILE, "$workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt");
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    if ($flds[0] =~ /[02468]$/) {
	$readName = $flds[0];
	($readNum) = ($readName =~ /^..(\d+)/);
	$readOffsetInSuperRead = $flds[2];
	$readOriInSuperRead = $flds[3]; }
    else {
	if ($readOriInSuperRead =~ /[Ff]/) {
	    $minOffset = $readOffsetInSuperRead;
	    $maxOffset = $flds[2];
	    $superReadOri = "F"; }
	else {
	    $minOffset = $flds[2];
	    $maxOffset = $readOffsetInSuperRead;
	    $superReadOri = "R"; }
	$joiningSuperRead{$readName} = substr ($superReadSequence{$flds[1]}, $minOffset, $maxOffset - $minOffset);
	if ($superReadOri =~ /[Rr]/) {
	    $joiningSuperRead{$readName} = revComp($joiningSuperRead{$readName}); 
	}
    }
}
close (FILE);

# Store the contig info in the order they appear in the scaffolds
# Columns of the input are (1) contig (2) scaff (3) min contig offset in scaff
#  (4) max contig offset in scaff (5) orientation of contig in scaff ('f' or 'r')
open (FILE, "$CeleraTerminatorDirectory/genome.posmap.ctgscf");
$orderedContigNum = 0;
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $scaff = $flds[1];
    # The following is so that the scaffs are output in order
    # It is only executed if this is the first time the scaff is
    # seen in the input file
    if (! $scfStr{$scaff}) {
	push (@scaffsToDo, $scaff); }
    ++$orderedContigNum;
    $orderedContigNum{$flds[0]} = $orderedContigNum; # Gives a number to the contig in order seen
    # The following adds the contig info to a string containing all the info for a scaffold
    # contigName, contigBegin, contigEnd, contigOri (all in relation to scaffold)
    $scfStr{$scaff} .= "$flds[0] @flds[2..4] "; }
close (FILE);

@connectedRead1s = keys %joiningSuperRead;
@connectedRead1s = sort spclByNum @connectedRead1s;
# This loop generates mergedContig which has triples of
# (frontFauxRead contig contigOri) for all the original contigs appearing in
# a final contig that has at least one merge, except that the
# first frontFauxRead isn't put in it; the domain is the set
# of contigs which are first in the final contigs
# mergedContigName is the new name of the contig, which has a string of
# contigNumbers and their scaffold orientations ('f' or 'r'),
# separated by underscores ('_'); the domain is the same as above
for (@connectedRead1s) {
    $readName = $_;
    $contig1 = $contig{$readName};
    $mateReadName = getMateReadName($readName);
    $contig2 = $contig{$mateReadName};
    # If the following is satisfied then contig1 has NOT been merged with the prior contig
    # $firstMergedContig{$contig1} is the first contig of the group containing $contig1
    if (! $firstMergedContig{$contig1}) {
	$tempContig = $firstMergedContig{$contig1} = $contig1;
	$mergedContig{$tempContig} = "$contig1 $contigOriFromReadName{$readName}";
	$mergedContigName{$tempContig} = $contig1 . $contigOriFromReadName{$readName}; }
    $tempContig = $firstMergedContig{$contig2} = $firstMergedContig{$contig1};
    $mergedContig{$tempContig} .= " $readName $contig2 $contigOriFromReadName{$mateReadName}";
    $mergedContigName{$tempContig} .= ("_" . $contig2 . $contigOriFromReadName{$mateReadName});
}

$contigSequenceInputFile = "$CeleraTerminatorDirectory/genome.ctg.fasta";
$basesPerLine = getBasesPerLine ($contigSequenceInputFile);
open (FILE, $contigSequenceInputFile);
open (OUTFILE, ">genome.ctg.fasta");
# This loops outputs all the contigs which are not merged with any other
while ($line = <FILE>) {
    chomp ($line);
    if ($line =~ /^>/) {
	($contig) = ($line =~ /^>\D*(\d+)/);
	if ($firstMergedContig{$contig}) {
	    $on = 0; }
	else {
	    $on = 1; }
    }
    if ($on) {
	print OUTFILE "$line\n"; }
    if ($line !~ /^>/) { # Add the sequence to that of the current contig (always
	# needed since we need to output scaffolds later
	$contigSeqLine{$contig} .= $line; }
}
close (FILE);

# Output the new (merged) contigs for those that actually had a merge
@mergedContigIndices = keys %mergedContig;
for (@mergedContigIndices) {
    $origContig = $_;
    $newContigInfo = $mergedContig{$origContig};
    $newContigName = $mergedContigName{$origContig};
    print OUTFILE ">jtg$newContigName\n";
    @flds = split (" ", $newContigInfo);
    $contig = $flds[0];
    $ori = $flds[1];
    # Get the appropriate sequence for the first contig based on the contig orientation
    if ($ori =~ /[fF]/) {
	$newContigSeq{$newContigName} = $contigSeqLine{$contig}; }
    else {
	$newContigSeq{$newContigName} = revComp ($contigSeqLine{$contig}); }
    for ($i=2; $i<$#flds; $i+=3) {
	$readName = $flds[$i];
	$contig = $flds[$i+1];
	$ori = $flds[$i+2];
	# The following kills off the sequence covered by the first faux read
	substr ($newContigSeq{$newContigName}, - $readLengths{$readName}) = "";
	$newContigSeq{$newContigName} .= $joiningSuperRead{$readName}; # Add the joining sequence
	if ($ori =~ /[fF]/) { # Taking off the second faux read from the (oriented) seq of the 2nd contig
	    $tstr = substr ($contigSeqLine{$contig}, $readLengths{getMateReadName($readName)}); }
	else {
	    $tstr = substr (revComp ($contigSeqLine{$contig}), $readLengths{getMateReadName($readName)}); }
	$newContigSeq{$newContigName} .= $tstr; # And now adding the sequence of the 2nd contig
    }
    outputFasta ($newContigSeq{$newContigName}, $basesPerLine);
}
close (OUTFILE);

# Output the scaffold sequence
open (OUTFILE, ">genome.scf.fasta");
open (CTG_SCF_FILE, ">genome.posmap.ctgscf");
for (@scaffsToDo) {
    $scaff = $_;
    $scaffHdr = ">jcf$scaff\n"; # Always changes the scaffold header to jcf
    $scaffSeq = "";
    $scaffStr = $scfStr{$scaff};
    @flds = split (" ", $scaffStr);
    $scaffOutputOffset = 0;
    $scaffCtgScfOutput = "";
    for ($i=0; $i<$#flds; $i+=4) {
	$contig = $flds[$i];
	if (($i-2 >= 0)) { # Add the appropriate number 'N's if the contig is not the first in the scaffold
	    # The first condition says that the contig was not merged with another, and the second
	    # says it was but was the first contig of the merged contig
	    # The combined condition only excludes contigs merged with others that are not the first
	    # contig in a merged contig
	    if ((! $firstMergedContig{$contig}) ||
		($firstMergedContig{$contig} eq $contig)) {
		$numNs = $flds[$i+1]-$flds[$i-2];
		$Ns = "N" x $numNs; # Add 'N's to the scaffold before the contig we're about to add
		$scaffOutputOffset += $numNs;
		$scaffSeq .= $Ns; } }
	if ($mergedContigName{$contig}) { # It's the first contig of a merged contig
	    $scaffSeq .= $newContigSeq{$mergedContigName{$contig}};
	    $scaffCtgScfOutput .= "$mergedContigName{$contig} $scaff $scaffOutputOffset ";
	    $scaffOutputOffset += length ($newContigSeq{$mergedContigName{$contig}});
	    $scaffCtgScfOutput .= "$scaffOutputOffset f\n";
	}
	elsif (! $firstMergedContig{$contig}) { # It's not merged with anything
	    if ($flds[$i+3] eq "f") {
		$scaffSeq .= $contigSeqLine{$contig}; }
	    else {
		$scaffSeq .= revComp($contigSeqLine{$contig}); }
	    $scaffCtgScfOutput .= "$contig $scaff $scaffOutputOffset ";
	    $scaffOutputOffset += length ($contigSeqLine{$contig});
	    $scaffCtgScfOutput .= ("$scaffOutputOffset " . $flds[$i+3] . "\n");
	}
    }
    print CTG_SCF_FILE $scaffCtgScfOutput;
    print OUTFILE $scaffHdr;
    outputFasta ($scaffSeq, $basesPerLine);
}
close (OUTFILE);
close (CTG_SCF_FILE);

sub spcl
{
    return ($readNum{$a} <=> $readNum{$b});
}

sub revComp
{
    my ($str) = @_;

    $str =~ tr/acgtACGT/tgcaTGCA/;
    $str = reverse ($str);
    return ($str);
}

sub getMateReadName
{
    my ($readName) = @_;
    my ($prefix, $last);

    ($prefix, $last) = ($readName =~ /^(.+)(.)$/);
    if ($last % 2 == 0) {
	++$last; }
    else {
	--$last; }
    $readName = $prefix . $last;
    return ($readName);
}

sub spclByNum
{
    return ($orderedContigNum{$contig{$a}} <=> $orderedContigNum{$contig{$b}});
}

# The following looks through a fasta file finding the standard line length used in the file
# To be used with Celera output; (must have some sequence that requires more than one line
# This goes through the file and looks for the first case where a fasta sequence is not followed
# by a header line and then uses the length of this first line to determine the number of bases
# per line used for the fasta file
sub getBasesPerLine
{
    my ($fn) = @_;

    open (FILE, $fn);
    while ($line = <FILE>) {
	chomp ($line);
	if ($line =~ /^>/) {
	    $state = 0; }
	else {
	    if ($state == 0) {
		$lineLen = length ($line);
		if ($lineLen > $maxLineLen) {
		    $maxLineLen = $lineLen; }
		++$state; }
	    else {
		last; } } }
    close (FILE);
    return ($maxLineLen);
}

sub outputFasta
{
    my ($seq, $lineLen) = @_;
    my ($seqLen, $i);
    $seqLen = length ($seq);
    for ($i=0; $i<$seqLen; $i+=$lineLen) { # Outputting using newlines
	print OUTFILE substr($seq, $i, $lineLen), "\n"; }
}

sub processArgs
{
    my ($arg);
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
	if (-d $arg) {
	    $CeleraTerminatorDirectory = $arg;
	    next; }
	if ($arg =~ /^contig_end_pairs/) {
	    $contigEndPairsFile = $arg;
	    next; }
	if ($arg eq "-contig-end-pairs-file") {
	    ++$i;
	    $contigEndPairsFile = $ARGV[$i];
	    next; }
	if ($arg eq "-working-directory") {
	    ++$i;
	    $workingDirectory = $ARGV[$i];
	    next; }
    }
}

sub byNum
{
    return ($b <=> $a);
}

