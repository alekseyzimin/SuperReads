#!/usr/bin/env perl
use File::Basename;
$exeDir = dirname ($0);
&processArgs;
$scaffoldFastaFile = $ARGV[0];
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
$cmd = "$exeDir/closeGapsLocally.perl -s $jellyfishHashSize --Celera-terminator-directory $CeleraTerminatorDirectory $readFileStr --output-directory $outputDirectory --min-kmer-len $minKMerLen --max-kmer-len $maxKMerLen --num-threads $numThreads --contig-length-for-joining $contigLengthForJoining --contig-length-for-fishing $contigLengthForFishing --reduce-read-set-kmer-size $reduceReadSetKMerSize";
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
    $CeleraTerminatorDirectory = ".";
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
	if ($arg eq "--Celera-terminator-directory") {
	    ++$i;
	    $CeleraTerminatorDirectory = $ARGV[$i];
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
		$CeleraTerminatorDirectory = $arg;
		next; }
	    $outputDirectory = $arg;
	    $cmd = "\\rm -r $outputDirectory";
	    print "$cmd\n"; system ($cmd);
	    next; }
	$outputDirectory = $arg;
    }
    @kmerLens = sort byNum @kmerLens;
    if ($#kmerLens >= 0) {
	$maxKMerLen = $kmerLens[0]; }
    if ($#kmerLens >= 1) {
	$minKMerLen = $kmerLens[1]; }
#   if (! $jellyfishHashSize) {
#	print STDERR "You must enter a jellyfish hash size. Bye!\n";
#	&reportUsage; }
#    if (! $CeleraTerminatorDirectory) {
#	print STDERR "You must enter a 9-terminator directory from a Celera run. Bye!\n";
#	&reportUsage; }
}
