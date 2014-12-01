#!/usr/bin/env perl
# This program is used to close gaps in a scaffold sequence file.
# Example invocation:
# closeGapsInScaffFastaFile.perl --scaffold-fasta-file 9-terminator/genome.scf.fa --min-kmer-len 17 --max-ker-len 31 --work-directory . --num-threads 16 --output-directory outputDir --reads-file pe.cor.fa --reads-file sj.cor.fa
#
# There are no args, only flags (mostly) with arguments. They are as follows:
#
# Required flags:
# --scaffold-fasta-file filename : file containing the scaffold sequences
# --reads-file filename : specify a read file to use (multiple files allowed,
#                            so long as the flag is repeated)
#
# Flags for required values which have defaults (i.e. flag not necessary)
# --work-directory dir : directory where the work will occur (default: '.')
# --min-kmer-len # : specify the min kmer len used (default: 19)
# --max-kmer-len # : specify the max kmer len used (default: 85)
# --output-directory dir : specify the output directory (default: 10-gapclose)
# --num-threads # : specify the number of threads (default: 1)
# -t # : same as --num-threads #
# -s # : the jellyfish hash size (default: 200000000)
# --contig-length-for-joining # : The length of sequence at the ends of the contigs
#                     which create the faux mate pairs which are joined (default: 100)
# --contig-length-for-fishing # : The length of sequence at the ends of the contigs
#                     to be used to find reads which might fit in the gaps (default: 100)
# --maxnodes # : The maximum number of nodes allowed when trying to join the
#                     faux reads (default: 200000)
# --reduce-read-set-kmer-size # : The k-mer size for fishing reads into buckets.
#                     (default: 21)
#
# The flags may be in any order.

use File::Basename;
$exeDir = dirname ($0);
&processArgs;
$cmd = "$exeDir/splitFileAtNs $scaffoldFastaFile > genome.ctg.fasta";
print "$cmd\n"; system ($cmd);

open (FILE, "splitScaffoldPlacementFile.txt");
open (OUTFILE, ">genome.posmap.ctgscf");
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $flds[0] =~ s/^...//;
    $flds[1] =~ s/^...//;
    print OUTFILE "@flds\n";
}
close (FILE); close (OUTFILE);

open (FILE, "genome.posmap.ctgscf");
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    push (@contigs, $flds[0]);
    push (@begins, $flds[2]);
    push (@ends, $flds[3]);
}
close (FILE);

open (OUTFILE, ">genome.asm");
for ($i=0; $i<$#contigs; ++$i) {
    $ctg1 = $contigs[$i];
    $ctg2 = $contigs[$i+1];
    $mean = $begins[$i+1] - $ends[$i];
    $std = $mean * .1;
    if ($std < 100) {
	$std = 100; }
    $std = int ($std);
    print OUTFILE "{SCF\nct1:$ctg1\nct2:$ctg2\nmea:$mean\nstd:$std\n";
}

for (@readsFiles) {
    $readFile = $_;
    $readFileStr .= "--reads-file $readFile "; }
$cmd = "$exeDir/closeGapsLocally.perl -s $jellyfishHashSize --Celera-terminator-directory $workDirectory $readFileStr --output-directory $outputDirectory --min-kmer-len $minKMerLen --max-kmer-len $maxKMerLen --num-threads $numThreads --contig-length-for-joining $contigLengthForJoining --contig-length-for-fishing $contigLengthForFishing --reduce-read-set-kmer-size $reduceReadSetKMerSize";
print "$cmd\n"; system ($cmd);

sub processArgs
{
    my ($arg, $tfile, $cmd, @kmerLens, $i);
    $reduceReadSetKMerSize = 21;
    $kUnitigContinuationNumber = 2;
    $maxKMerLen = 85;
    $jellyfishHashSize = 200000000;
    $minKMerLen = 19;
    $numThreads = 1;
    $maxFishingKMerCount = 1000;
    $maxReadsInMemory = 100000000;
    $contigLengthForJoining = $contigLengthForFishing = 100;
    $maxNodes = 10000;
    $fauxInsertMean = 600;
    $fauxInsertStdev = 200;
    $numStdevsAllowed = 5;
    $workDirectory = ".";
    $outputDirectory = "10-gapclose";
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
	if ($arg eq "--scaffold-fasta-file") {
	    ++$i;
	    $scaffoldFastaFile = $ARGV[$i];
	    next; }
	if ($arg eq "--min-kmer-len") {
	    ++$i;
	    $minKMerLen = $ARGV[$i];
	    next; }
	if ($arg eq "--max-kmer-len") {
	    ++$i;
	    $maxKMerLen = $ARGV[$i];
	    next; }
	if ($arg eq "--work-directory") {
	    ++$i;
	    $workDirectory = $ARGV[$i];
	    next; }
	if ($arg eq "--reads-file") {
	    ++$i;
	    push (@readsFiles, $ARGV[$i]);
	    next; }
	if ($arg eq "--output-directory") {
	    ++$i;
	    $outputDirectory = $ARGV[$i];
	    next; }
	if (($arg eq "--num-threads") || ($arg eq "-t")) {
	    ++$i;
	    $numThreads = $ARGV[$i];
	    next; }
	if ($arg eq "--reduce-read-set-kmer-size") {
	    ++$i;
	    $reduceReadSetKMerSize = $ARGV[$i];
	    next; }
	if ($arg eq "--contig-length-for-fishing") {
	    ++$i;
	    $contigLengthForFishing = $ARGV[$i];
	    next; }
	if ($arg eq "--contig-length-for-joining") {
	    ++$i;
	    $contigLengthForJoining = $ARGV[$i];
	    next; }
	if ($arg =~ /-h/i) {
	    &reportUsage; }
#	if ($arg eq "--maxnodes") {
#	    ++$i;
#	    $maxNodes = $ARGV[$i];
#	    next; }
#	if ($arg eq "--max-reads-in-memory") {
#	    ++$i;
#	    $maxReadsInMemory = $ARGV[$i];
#	    next; }
	# The following is the max number of (Celera assembler) standard
	# deviations the joining sequence length can be from the length
	# predicted by the Celera assembler
#	if ($arg eq "--num-stdevs-allowed") {
#	    ++$i;
#	    $numStdevsAllowed = $ARGV[$i];
#	    next; }
#       if ($arg eq "--faux-insert-mean") {
#            ++$i;
#            $fauxInsertMean = $ARGV[$i];
#            next; }
#        if ($arg eq "--faux-insert-stdev") {
#            ++$i;
#            $fauxInsertStdev = $ARGV[$i];
#            next; }
        if ($arg eq "-s") {
            ++$i;
            $jellyfishHashSize = $ARGV[$i];
            next; }
#	if ($arg eq "--keep-directories") {
#	    $keepDirectoriesFlag = $arg;
#	    next; }
	if (-f $arg) {
	    push (@readsFiles, $arg);
	    next; }
#	if ($arg =~ /^\d+/) {
#	    push (@kmerLens, $arg);
#	    next; }
	if (-d $arg) {
	    $tfile = "$arg/genome.posmap.scflen";
	    if (-e $tfile) {
		$workDirectory = $arg;
		next; }
	    $outputDirectory = $arg;
	    $cmd = "\\rm -r $outputDirectory";
	    print "$cmd\n"; system ($cmd);
	    next; }
	$outputDirectory = $arg;
    }
    if ($scaffoldFastaFile !~ /\S/) {
	print STDERR "You must enter a scaffold sequence file. Bye!\n";
	&reportUsage; }
    if (! -e $scaffoldFastaFile) {
	print STDERR "Scaffold sequence file \"$scaffoldFastaFile\" doesn't exist. Bye!\n";
	&reportUsage; }
    if (-s $scaffoldFastaFile == 0) {
	print STDERR "Scaffold sequence file \"$scaffoldFastaFile\" has 0 size. Bye!\n";
	&reportUsage; }
    if ($#readsFiles < 0) {
	print STDERR "You must enter at least one read file with the flag --reads-file. Bye\n";
	&reportUsage; }
    $isGood = 1;
    for (@readsFiles) {
	$treadFile = $_;
	if (! -e $treadFile) {
	    print "Read file \"$treadFile\" doesn't exist.\n";
	    $isGood = 0;
	    next; }
	if (-s $treadFile == 0) {
	    print STDERR "Read file\"$treadFile\" has 0 size.\n";
	    $isGood = 0; }
    }
    if (! $isGood) {
	print STDERR "Bye!\n";
	&reportUsage; }
}

sub reportUsage
{
    open (FILE, $0);
    $line = <FILE>;
    while ($line = <FILE>) {
	last unless ($line =~ /^\#/);
	chomp ($line);
	($line2) = ($line =~ /^..(.*)$/);
	print STDERR "$line2\n";
    }
    close (FILE);
    exit (1);
}

