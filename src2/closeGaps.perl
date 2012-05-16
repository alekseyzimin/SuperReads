#!/usr/bin/env perl
# Example invocation:
# closeGaps.perl 15 31 origReads.renamed.fasta .../CA/9-terminator outputDir
# where 15 is the min k-mer size to use (default is 15)
#       31 is the max k-mer size to use (default is 31)
#       origReads.renamed.fasta is the name of the (input) read file
#       .../CA/9-terminator is the directory with the Celera output
#       outputDir is the directory where the output will be sent
# The args may be in any order.
use Cwd;
use File::Basename;
$exeDir = dirname ($0);
$meanForGap = 500;
$stdevForGap = 200;
# $tempExeDir = "/home/tri/superReadPipeline/SuperReads/build/install_root/bin";
$tempExeDir = $exeDir;
#Create absolute paths where necessary
# Must allow one to specify the dir with the hash and input k-unitigs file
# Must allow specification of the min k-unitig continuation values (default 2)
$cwd = cwd;
&processArgs;
print "";
$CeleraTerminatorDirectory = returnAbsolutePath ($CeleraTerminatorDirectory);
$readsFile = returnAbsolutePath ($readsFile);
if ($dirForKUnitigs) {
    $dirForKUnitigs = returnAbsolutePath ($dirForKUnitigs); }
$localReadsFile = basename ($readsFile);
if (! -e $outputDirectory) {
    $cmd = "mkdir $outputDirectory"; print "$cmd\n"; system ($cmd); }
chdir ($outputDirectory);
$cmd = "$exeDir/part01.perl $CeleraTerminatorDirectory";
print "$cmd\n"; system ($cmd);
$cmd = "$exeDir/create_end_pairs.perl $CeleraTerminatorDirectory > contig_end_pairs.fa";
print "$cmd\n"; system ($cmd);
$cmd = "ln -s $readsFile .";
print "$cmd\n"; system ($cmd);
$cmd = "echo \"cc $meanForGap $stdevForGap\" > meanAndStdevByPrefix.cc.txt";
print "$cmd\n"; system ($cmd);

&runMainLoop;

$cmd = "$exeDir/getSequenceForClosedGaps.perl $CeleraTerminatorDirectory $localReadsFile $minKMerLen $maxKMerLen";
print "$cmd\n"; system ($cmd);

&createMergedJoinFile;

&cleanUp;

sub cleanUp
{
    my ($k, $suffix, $cmd);
    for ($k=$maxKMerLen; $k>=$minKMerLen; $k--) {
	$suffix = $localReadsFile . "_" . $k . "_" . $kUnitigContinuationNumber;
	$cmd = "\\rm -rf work_$suffix"; system ($cmd);
	$cmd = "\\rm k_u_hash_${suffix}_0"; system ($cmd);
	$cmd = "\\rm  out.$suffix"; system ($cmd);
	$cmd = "\\rm joined.$suffix"; system ($cmd);
    }
}

sub createMergedJoinFile
{
    my ($outfile, $k, $suffix, $fn, $line, %wasOutput);
    $outfile = "joined.${localReadsFile}_all.txt";
    open (OUTFILE, ">$outfile");
    for ($k=$maxKMerLen; $k>=$minKMerLen; $k--) {
	$suffix = $localReadsFile . "_" . $k . "_" . $kUnitigContinuationNumber;
	$fn = "joined.$suffix";
	open (FILE, $fn);
	while ($line = <FILE>) {
	    next if ($wasOutput{$line});
	    print OUTFILE $line;
	    $wasOutput{$line} = 1; }
	close (FILE); }
    close (OUTFILE);
}

sub runMainLoop
{
    my ($k, $kMinus1, $kPlus1, $suffix, $cmd);
    for ($k=$maxKMerLen; $k>=$minKMerLen; $k--) {
	$kMinus1 = $k-1;
	$kPlus1 = $k+1;
	$suffix = $localReadsFile . "_" . $k . "_" . $kUnitigContinuationNumber;
	if ($dirForKUnitigs) {
	    $cmd = "ln -s $dirForKUnitigs/k_u_hash_${suffix}_0 .";
	    print "$cmd\n"; system ($cmd);
	    $cmd = "ln -s $dirForKUnitigs/k_unitigs_${suffix}.fa .";
	    print "$cmd\n"; system ($cmd); }
	else {
	    $cmd = "jellyfish count -m $k -t 16 -C -r -s 200000000 -o k_u_hash_${suffix} $readsFile";
	    print "$cmd\n"; system ("time $cmd");
	    $cmd = "$tempExeDir/create_k_unitigs --cont-on-low  --low-stretch=$kMinus1 -C -t 16 -l $kPlus1 -o k_unitigs_${suffix} -m $kUnitigContinuationNumber -M $kUnitigContinuationNumber k_u_hash_${suffix}_0";
	    print "$cmd\n"; system ("time $cmd"); }
	$cmd = "\\rm -rf out.$suffix"; system ($cmd);
	$cmd = "\\rm -rf work_$suffix"; system ($cmd);
	
	$cmd = "$tempExeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -minreadsinsuperread 1 -kunitigsfile k_unitigs_${suffix}.fa -low-memory -l $k -s 200000000 -t 1 -mkudisr 0 work_${suffix} contig_end_pairs.fa 1>>out.$suffix 2>>out.$suffix";
	print "$cmd\n"; system ("time $cmd");

	$cmd = "grep '^Num ' out.$suffix"; system ($cmd);

	$cmd = "cat work_$suffix/readPositionsInSuperReads | $exeDir/outputJoinedPairs.perl > joined.$suffix";
	print "$cmd\n"; system ("time $cmd");
    }
}


sub returnAbsolutePath
{
    my ($file) = @_;
    if ($file !~ /^\//) {
	$file = "$cwd/$file"; }
    return ($file);
}

sub processArgs
{
    my ($arg, $tfile, $cmd, @kmerLens, $i);
    $kUnitigContinuationNumber = 2;
    $maxKMerLen = 31;
    $minKMerLen = 15;
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
	if ($arg eq "-dir-for-kunitigs") {
	    ++$i;
	    $dirForKUnitigs = $ARGV[$i];
	    next; }
	if (-f $arg) {
	    $readsFile = $arg;
	    next; }
	if ($arg =~ /^\d+/) {
	    push (@kmerLens, $arg);
	    next; }
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
    if ($#kmerLens >= 2) {
	$kUnitigContinuationNumber = $kmerLens[2]; }
    if (! $CeleraTerminatorDirectory) {
	print STDERR "You must enter a 9-terminator directory from a Celera run. Bye!\n";
	&reportUsage; }
    if (! $readsFile) {
	print STDERR "You must enter the name of an (existing) input reads file. Bye!\n";
	&reportUsage; }
}

sub byNum
{
    return ($b <=> $a);
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

