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

open (FILE, "$workingDirectory/superReadSequences.fasta");
while ($superRead = <FILE>) {
    chomp ($superRead);
    $superRead =~ s/^.//;
    $superReadSequence = <FILE>;
    chomp ($superReadSequence);
    $superReadSequence{$superRead} = $superReadSequence;
}
close (FILE);


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

# This is done early so we can put the contigs in scaffold order
open (FILE, "$CeleraTerminatorDirectory/genome.posmap.ctgscf");
$orderedContigNum = 0;
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $scaff = $flds[1];
    if (! $scfStr{$scaff}) {
	push (@scaffsToDo, $scaff); }
    ++$orderedContigNum;
    $orderedContigNum{$flds[0]} = $orderedContigNum;
    $scfStr{$scaff} .= "$flds[0] @flds[2..4] "; }
close (FILE);

@connectedRead1s = keys %joiningSuperRead;
@connectedRead1s = sort spclByNum @connectedRead1s;
for (@connectedRead1s) {
    $readName = $_;
    $contig1 = $contig{$readName};
    $mateReadName = getMateReadName($readName);
    $contig2 = $contig{$mateReadName};
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
    if ($line !~ /^>/) {
	$contigSeqLine{$contig} .= $line; }
}
close (FILE);

@mergedContigIndices = keys %mergedContig;
for (@mergedContigIndices) {
    $origContig = $_;
    $newContigInfo = $mergedContig{$origContig};
    $newContigName = $mergedContigName{$origContig};
    print OUTFILE ">jtg$newContigName\n";
    @flds = split (" ", $newContigInfo);
    $contig = $flds[0];
    $ori = $flds[1];
    if ($ori =~ /[fF]/) {
	$newContigSeq{$newContigName} = $contigSeqLine{$contig}; }
    else {
	$newContigSeq{$newContigName} = revComp ($contigSeqLine{$contig}); }
    for ($i=2; $i<$#flds; $i+=3) {
	$readName = $flds[$i];
	$contig = $flds[$i+1];
	$ori = $flds[$i+2];
	substr ($newContigSeq{$newContigName}, - $readLengths{$readName}) = "";
	$newContigSeq{$newContigName} .= $joiningSuperRead{$readName};
	if ($ori =~ /[fF]/) {
	    $tstr = substr ($contigSeqLine{$contig}, $readLengths{getMateReadName($readName)}); }
	else {
	    $tstr = substr (revComp ($contigSeqLine{$contig}), $readLengths{getMateReadName($readName)}); }
	$newContigSeq{$newContigName} .= $tstr;
    }
    outputFasta ($newContigSeq{$newContigName}, $basesPerLine);
}
close (OUTFILE);

open (OUTFILE, ">genome.scf.fasta");
open (CTG_SCF_FILE, ">genome.posmap.ctgscf");
for (@scaffsToDo) {
    $scaff = $_;
    $scaffHdr = ">jcf$scaff\n";
    $scaffSeq = "";
    $scaffStr = $scfStr{$scaff};
    @flds = split (" ", $scaffStr);
    $scaffOutputOffset = 0;
    $scaffCtgScfOutput = "";
    for ($i=0; $i<$#flds; $i+=4) {
	$contig = $flds[$i];
	if (($i-2 >= 0)) {
	    if ((! $firstMergedContig{$contig}) ||
		($firstMergedContig{$contig} eq $contig)) {
		$numNs = $flds[$i+1]-$flds[$i-2];
		$Ns = "N" x $numNs;
		$scaffOutputOffset += $numNs;
		$scaffSeq .= $Ns; } }
	if ($mergedContigName{$contig}) {
	    $scaffSeq .= $newContigSeq{$mergedContigName{$contig}};
	    $scaffCtgScfOutput .= "$mergedContigName{$contig} $scaff $scaffOutputOffset ";
	    $scaffOutputOffset += length ($newContigSeq{$mergedContigName{$contig}});
	    $scaffCtgScfOutput .= "$scaffOutputOffset f\n";
	}
	elsif (! $firstMergedContig{$contig}) {
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
    for ($i=0; $i<$seqLen; $i+=$lineLen) {
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

