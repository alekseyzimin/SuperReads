#! /usr/bin/env perl
# Report expected joining sequence length for each Celera gap where, for the
# correct orientation, the gap has a unique length implied by the mate pairs
# as reported by the Celera assembler
# Mandatory argument
# The 9-terminator directory of the Celera run
# Mandatory option:
# --contig-end-seq-file filename ; the file that contains the contig end
#    sequences to be joined by the gap closer
# Optional flag:
# --reduced-column-output : Output just the gapNum mean and stdev for each gap
#   

# $defaultMean = 600;
# $defaultStdev = 200;
&processArgs;
# $contigEndSeqFile = "/genome3/raid/tri/localJoinTests/mouse/testMouse121130.1000fishingCutoff_10kNodes/contig_end_pairs.100.fa";
$infile = $contigEndSeqFile;
# The following assumes that the contig sequence file has alternating
# lines of header with sequence
$cmd = "grep \"^>\" $infile |";
# open (FILE, $cmd);
open (FILE, $infile);
$gapNum = 0;
while ($line = <FILE>) {
    chomp ($line);
    @flds = split(" ", $line);
    $ctg1 = $flds[1]; $ctg1ori = $flds[2];
    $seq1 = <FILE>;
    $contigLen1[$gapNum] = length ($seq1)-1;
    $line = <FILE>;
    @flds = split (" ", $line);
    $ctg2 = $flds[1]; $ctg2ori = $flds[2];
    $seq2 = <FILE>;
    $contigLen2[$gapNum] = length ($seq2)-1;
    $str1 = "$ctg1 $ctg2";
    $str2 = "$ctg2 $ctg1";
    $gapNum{$str1} = $gapNum{$str2} = $gapNum;
    $contig1[$gapNum] = $ctg1; $contig2[$gapNum] = $ctg2;
    $contig1ori[$gapNum] = $ctg1ori; $contig2ori[$gapNum] = $ctg2ori;
    ++$gapNum;
}

# open (FILE, "finishedSeqLens.txt"); while ($line = <FILE>) { chomp ($line); ($gapNum) = ($line =~ /^\s*(\d+)\s/); $line[$gapNum] = $line; } close (FILE);
open (FILE, "$CeleraTerminatorDirectory/genome.posmap.ctglkg");
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $str = "@flds[0..1]";
    next unless ($gapNum{$str} =~ /\d/);
#    next unless ($flds[8] eq "U"); # Logic now refined and later in program
    next unless ($flds[4] eq ".");
    $gapNum = $gapNum{$str};
    $fail = 0;
    if ($flds[0] == $contig1[$gapNum]) {
	if ($contig1ori[$gapNum] eq "f") {
	    if ($contig2ori[$gapNum] eq "f") { $expectedOri = "N"; $fail = 1 unless ($flds[8] eq "U"); }
	    else { $expectedOri = "I"; } }
	else {
	    if ($contig2ori[$gapNum] eq "f") { $expectedOri = "O"; }
	    else { $expectedOri = "A"; } }
    }
    else { # Contigs in the join file are ctg after the gap first, then before
	if ($contig1ori[$gapNum] eq "f") {
	    if ($contig2ori[$gapNum] eq "f") { $expectedOri = "A"; }
	    else { $expectedOri = "I"; } }
	else {
	    if ($contig2ori[$gapNum] eq "f") { $expectedOri = "O"; }
	    else { $expectedOri = "N"; } }
    }
    next unless ($flds[2] eq $expectedOri);
    next if ($fail);
    $flds[5] += $contigLen1[$gapNum] + $contigLen2[$gapNum];
    if (! $reducedColumnOutput) {
	$tline = "$gapNum @flds[2..3] @flds[5..8]\n"; }
    else {
	$tline = "$gapNum @flds[5..6]\n"; }
    if ($flds[8] eq "U") { $Ulines[$gapNum] .= $tline; ++$numULines[$gapNum]; }
    else { $ABlines[$gapNum] .= $tline; ++$numABLines[$gapNum]; }
    if ($gapNum > $maxGapNum) {
	$maxGapNum = $gapNum; }
#    print "$gapNum $str $contig1[$gapNum] $contig1ori[$gapNum] $contig2[$gapNum] $contig2ori[$gapNum] @flds[2..3] @flds[5..8]\n";
    
}
close (FILE);
for ($i=0; $i<=$maxGapNum; $i++) {
    if ($Ulines[$i]) { print $Ulines[$i] if ($numULines[$i] == 1); }
    else { print $ABlines[$i] if ($numABLines[$i] == 1); }
}
# | perl -ane 'print if (($F[4] <= 200) && ($F[4] >= -136) && ($F[5] <= 50));' | ~/ccat.perl | sort -n -k1,1 > contigFinishedSeqJoinInfo.gapNum.finLen.ctglkgLine.txt

sub processArgs
{
    my ($i);

    for ($i=0; $i<=$#ARGV; $i++) {
	if ($ARGV[$i] eq "--contig-end-seq-file") {
	    ++$i;
	    $contigEndSeqFile = $ARGV[$i];
	    next; }
#	elsif ($ARGV[$i] eq "--default-mean") {
#	    ++$i;
#	    $defaultMeanGiven = 1;
#	    $defaultMean = $ARGV[$i];
#	    next; }
#	elsif ($ARGV[$i] eq "--default-stdev") {
#	    ++$i;
#	    $defaultStdevGiven = 1;
#	    $defaultStdev = $ARGV[$i];
#	    next; }
#	elsif ($ARGV[$i] eq "--keep-if-no-direct-mate-joins") {
#	    $keepIfNoDirectMateJoins = 1;
#	    next; }
	elsif ($ARGV[$i] eq "--reduced-column-output") {
	    $reducedColumnOutput = 1;
	    next; }
	# Add for a help statement here
	elsif ($ARGV[$i] =~ /^\-h/i) { &reportUsage; }
	elsif ($ARGV[$i] =~ /^\-\-h/i) { &reportUsage; }
	elsif ($ARGV[$i] =~ /^\-/) { print STDERR "Unrecognized flag ",$ARGV[$i], ".\n"; &reportUsage; }
	else { 
	    $CeleraTerminatorDirectory = $ARGV[$i]; }
    }
#    if ($defaultMeanGiven && ! $defaultStdevGiven) {
#	$defaultStdev = $defaultMeanGiven / 3; }
    if (! $contigEndSeqFile) {
	print STDERR "You must use the --contig-end-seq-file to report the name of the fasta file of contig ends used when joining.\n";
	$error = 1; }
    elsif (! -e $contigEndSeqFile) {
	print STDERR "The contig end fasta file '$contigEndSeqFile' doesn't exist.\n";
	$error = 1; }
    if (! $CeleraTerminatorDirectory) {
	print STDERR "You must specify the Celera terminator directory as an argument on the command line.\n";
	$error = 1; }
    elsif (! -e $CeleraTerminatorDirectory) {
	print STDERR "The Celera terminator directory '$CeleraTerminatorDirectory' doesn't exist.\n";
	$error = 1; }
    elsif (! -d $CeleraTerminatorDirectory) {
	print STDERR "The Celera terminator directory '$CeleraTerminatorDirectory' isn't a directory.\n";
	$error = 1; }
    if ($error) {
	&reportUsage; }
}

sub reportUsage
{
    open (FILE, $0);
    $line = <FILE>;
    while ($line = <FILE>) {
	chomp ($line);
	last unless ($line =~ /^\#/);
	($line) = ($line =~ /^..(.*)$/);
	print "$line\n"; }
    close (FILE);
    exit (0);
}

