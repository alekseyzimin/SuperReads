#! /usr/bin/env perl
# Report expected joining sequence length for each Celera gap where, for the
# correct orientation, the gap has a unique length implied by the mate pairs
# as reported by the Celera assembler
# Mandatory argument
# The 9-terminator directory of the Celera run
# Mandatory option:
# --contig-end-seq-file filename ; the file that contains the contig end
#    sequences to be joined by the gap closer

&processArgs;
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
    $ctg1 = $flds[1]; 
    $seq1 = <FILE>;
    $contigLen1[$gapNum] = length ($seq1)-1;
    $line = <FILE>;
    @flds = split (" ", $line);
    $ctg2 = $flds[1]; 
    $seq2 = <FILE>;
    $contigLen2[$gapNum] = length ($seq2)-1;
    $str1 = "$ctg1 $ctg2";
    $str2 = "$ctg2 $ctg1";
    $gapNum{$str1} = $gapNum{$str2} = $gapNum;
    ++$gapNum;
}

$cmd = "grep --text -A 10 -E \"^\\\{\(SCF|CTP\)\$\" $CeleraTerminatorDirectory/genome.asm |";
open (FILE, $cmd);
while ($line = <FILE>) {
    next unless ($line =~ /^ct1:/);
    chomp ($line);
    ($ct1) = ($line =~ /^....(\S+)\s*$/);
    while ($line = <FILE>) {
	last if ($line =~ /^ct2:/); }
    chomp ($line);
    ($ct2) = ($line =~ /^....(\S+)\s*$/);
    while ($line = <FILE>) {
	last if ($line =~ /^mea:/); }
    chomp ($line);
    ($mean) = ($line =~ /^....(\S+)\s*$/);
    while ($line = <FILE>) {
	last if ($line =~ /^std:/); }
    chomp ($line);
    ($stdev) = ($line =~ /^....(\S+)\s*$/);
    while ($line = <FILE>) {
	last if ($line =~ /^ori:/); }
    chomp ($line);
#    ($ori) = ($line =~ /^....(\S+)\s*$/);
    next if ($ct1 eq $ct2);
    $str = "$ct1 $ct2";
    next unless ($gapNum{$str} =~ /\d/);
    $gapNum = $gapNum{$str};
    $mean[$gapNum] = $mean;
    $stdev[$gapNum] = $stdev;

    if ($gapNum > $maxGapNum) {
	$maxGapNum = $gapNum; }
    
}
close (FILE);
for ($i=0; $i<=$maxGapNum; $i++) {
    if ($mean[$i] =~ /\d/) {
	$mean[$i] += ($contigLen1[$i] + $contigLen2[$i]);
	}else{
	$mean[$i] = ($contigLen1[$i] + $contigLen2[$i])+500;#assume 500 bp gap
        $stdev[$i]=200;#assume 200bp stdev
	}
    print "$i $mean[$i] $stdev[$i]\n";
}

sub processArgs
{
    my ($i);

    for ($i=0; $i<=$#ARGV; $i++) {
	if ($ARGV[$i] eq "--contig-end-seq-file") {
	    ++$i;
	    $contigEndSeqFile = $ARGV[$i];
	    next; }
	# Add for a help statement here
	elsif ($ARGV[$i] =~ /^\-h/i) { &reportUsage; }
	elsif ($ARGV[$i] =~ /^\-\-h/i) { &reportUsage; }
	elsif ($ARGV[$i] =~ /^\-/) { print STDERR "Unrecognized flag ",$ARGV[$i], ".\n"; &reportUsage; }
	else { 
	    $CeleraTerminatorDirectory = $ARGV[$i]; }
    }
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

