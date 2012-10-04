#!/usr/bin/env perl
# We assume the following is one line of sequence per faux read
# $CeleraTerminatorDirectory is of the form */9-terminator
# Input files:
# 1) contig_end_pairs.fa
# 2) joined.${readsFile}_${fileNum}_2 (fileNum = 15..31), e.g.
#          joined.origReads.renamed.fasta_21_2
#     $dir = work_${readsFile}_${fileNum}_2, e.g.
#     $dir = work_origReads.renamed.fasta_21_2
# 3) $dir/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt
# 4) $dir/superReadSequences.fasta
# 5) $CeleraTerminatorDirectory/genome.posmap.ctgscf
# 6) $CeleraTerminatorDirectory/genome.ctg.fasta
$maxKMerLen = 31;
$minKMerLen = 19;
$kUnitigContinuationNumber = 2;
$readsFile = "origReads.renamed.fasta";
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

for ($i=$maxKMerLen; $i>=$minKMerLen; $i--) {
    $fn = createFilenameFromNum ($i);
    $dir = createDirFromNum ($i);
    @read1sToDoInPass = ();
    open (FILE, $fn);
    while ($line = <FILE>) {
	chomp ($line);
	if (! $alreadyDone{$line}) {
	    push (@read1sToDoInPass, $line);
	    $alreadyDone{$line} = 1; }
	$line = <FILE>; # Just skipping the second read of the mate pair
    }
    close (FILE);
    next unless ($#read1sToDoInPass >= 0);
    undef (%mustGetSequence);
    for (@read1sToDoInPass) {
	$mustGetSequence{$_} = 1; }

    open (FILE, "$dir/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt");
    undef (%superReadToRead);
    undef (%readOffsetInSuperRead);
    undef (%readNum);
    undef (%readOriInSuperRead);
    undef (%readToSuperRead);
    while ($line = <FILE>) {
	chomp ($line);
	@flds = split (" ", $line);
	$readName = $flds[0];
	$otherRead = getMateReadName ($readName);
	if (($mustGetSequence{$readName}) || ($mustGetSequence{$otherRead})) {
	    ($readNum) = ($readName =~ /^..(\d+)/);
	    $readOffsetInSuperRead{$readName} = $flds[2];
	    $readNum{$readName} = $readNum;
	    $readOriInSuperRead{$readName} = $flds[3];
	    $readToSuperRead{$readName} = $flds[1];
	    if ($mustGetSequence{$readName}) {
		$superReadToRead{$flds[1]} = $readName; } }
    }
    close (FILE);

    for (@read1sToDoInPass) {
	$read1 = $_;
	$read2 = getMateReadName ($read1);
	$isBad = 0;
	if ($readOriInSuperRead{$read1} eq $readOriInSuperRead{$read2}) {
	    $isBad = 1; }
	if ($readOriInSuperRead{$read1} =~ /[Ff]/) {
	    $minOffset = $readOffsetInSuperRead{$read1};
	    $maxOffset = $readOffsetInSuperRead{$read2};
	    $superReadOri = "F"; }
	else {
	    $minOffset = $readOffsetInSuperRead{$read2};
	    $maxOffset = $readOffsetInSuperRead{$read1};
	    $superReadOri = "R"; }
	if ($minOffset > $maxOffset) {
	    $isBad = 1; }
	if ($isBad) {
	    delete ($alreadyDone{$read1});
	}
	else {
	    $minOffset{$read1} = $minOffset;
	    $maxOffset{$read1} = $maxOffset;
	    $superReadOri{$read1} = $superReadOri; }
    }

    undef (%superReadSequence);

    open (FILE, "$dir/superReadSequences.fasta");
    while ($superRead = <FILE>) {
	chomp ($superRead);
	$superRead =~ s/^.//;
	$superReadSequence = <FILE>;
	chomp ($superReadSequence);
	if ($superReadToRead{$superRead}) {
	    $superReadSequence{$superRead} = $superReadSequence; }
    }
    close (FILE);

    for (@read1sToDoInPass) {
	$read1 = $_;
	next unless ($alreadyDone{$read1});
	$neededSequence = $superReadSequence{$readToSuperRead{$read1}};
	$joiningSuperRead{$read1} = substr ($neededSequence, $minOffset{$read1}, $maxOffset{$read1} - $minOffset{$read1});
	if ($superReadOri{$read1} =~ /[Rr]/) {
	    $joiningSuperRead{$read1} = revComp($joiningSuperRead{$read1}); }
    }
}

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

@connectedRead1s = keys %alreadyDone;
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

sub createFilenameFromNum
{
    my ($fileNum) = @_;
    my ($fn);

    $fn = "joined.${readsFile}_${fileNum}_${kUnitigContinuationNumber}";
    return ($fn);
}

sub createDirFromNum
{
    my ($fileNum) = @_;
    my ($fn);

    $fn = "work_${readsFile}_${fileNum}_${kUnitigContinuationNumber}";
    return ($fn);
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
    my ($arg, @kmerLens);
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
	if (-d $arg) {
	    $CeleraTerminatorDirectory = $arg;
	    next; }
	if ($arg =~ /^contig_end_pairs/) {
	    $contigEndPairsFile = $arg;
	    next; }
	if (-f $arg) {
	    $readsFile = $arg;
	    next; }
	if ($arg eq "-reads-file") {
	    ++$i;
	    $readsFile = $ARGV[$i];
	    next; }
	if ($arg eq "-contig-end-pairs-file") {
	    ++$i;
	    $contigEndPairsFile = $ARGV[$i];
	    next; }
	if ($arg =~ /^\d+$/) {
	    push (@kmerLens, $arg); } }
    @kmerLens = sort byNum @kmerLens;
    if ($#kmerLens >= 0) {
        $maxKMerLen = $kmerLens[0]; }
    if ($#kmerLens >= 1) {
        $minKMerLen = $kmerLens[1]; }
    if ($#kmerLens >= 2) {
        $kUnitigContinuationNumber = $kmerLens[2]; }
}

sub byNum
{
    return ($b <=> $a);
}

