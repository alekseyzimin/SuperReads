#!/usr/bin/env perl
#
# This program is used to find reverse-pointing jumping reads
#
# Example invocation:
# findReversePointingJumpingReads.perl --min-kmer-len 17 --max-ker-len 31 --Celera-terminator-directory .../CA/9-terminator --num-threads 16 --output-directory outputDir --reads-file pe.cor.fa --reads-file sj.cor.fa
#
# There are no args, only flags (mostly) with arguments. They are as follows:
#
# Required flags:
# --jumping-library-read-file file : name of file with jumping-library reads 
# --Celera-terminator-directory dir : specify the Celera terminator directory
#                           where the assembly whose gaps must be closed exists
# --reads-file filename : specify a read file to use (multiple files allowed,
#                            so long as the flag is repeated)
# --output-directory dir : specify the output directory
# -s # : the jellyfish hash size
#
# Flags for required values which have defaults (i.e. flag not necessary)
# --min-kmer-len # : specify the min kmer len used (default: 17)
# --max-kmer-len # : specify the max kmer len used (default: 65)
# --num-threads # : specify the number of threads (default: 1)
# -t # : same as --num-threads #
# --maxnodes # : The maximum number of nodes allowed when trying to join the
#                     faux reads (default: 200000)
# --reduce-read-set-kmer-size # : The k-mer size for fishing reads into buckets.
#                     (default: 21)
# --keep-directories : Keep the local directories (default: false)
# --max-reads-in-memory : The maximum number of reads whose sequence can be kept
#                         in memory at one time (default: 100000000)
# --faux-insert-mean : The mean of the insert size used for the faux reads around
#                      a gap (default: 500)
# --faux-insert-stdev : The stdev of the insert size used for the faux reads around
#                      a gap (default: 200)
# --num-stdevs-allowed # : The maximum number of standard deviations that the
#                      length of the join can deviate from the estimate output
#                      by the Celera assembler. The standard deviation generated
#                      by the Celera assembler is used. (default: 3)
#
# --num-joins-per-directory # : How many joins to attempt in each directory (default: 100)
#
# The flags may be in any order.
use Cwd;
use File::Basename;
$exeDir = dirname ($0);
# Create absolute paths where necessary
# Must allow one to specify the dir with the hash and input k-unitigs file
# Must allow specification of the min k-unitig continuation values (default 2)
$cwd = cwd;
$exeDir = returnAbsolutePath ($exeDir);
&processArgs;
$shm = "/dev/shm";
if ((! -d $shm) || $keepDirectoriesFlag) {
    $subdir2 = "subdir2"; }
else {
    $tempTime = time;
    $subdir2 = "$shm/$tempTime"; }
print "";
$CeleraTerminatorDirectory = returnAbsolutePath ($CeleraTerminatorDirectory);
$jumpingLibraryReadFile = returnAbsolutePath ($jumpingLibraryReadFile);
for ($i=0; $i<=$#readsFiles; $i++) {
    $readsFile = $readsFiles[$i];
    $readsFile = returnAbsolutePath ($readsFile);
    $readsFiles[$i] = $readsFile; }

$localReadsFile = "localReadsFile";
print "outputDirectory = $outputDirectory\n";
if (! -e $outputDirectory) {
    $cmd = "mkdir $outputDirectory"; runCommandAndExitIfBad($cmd); }
$reverseComplementedJumpingReadFile = "$outputDirectory/reverseComplemented.jumpingReadFile.fa";
$reverseComplementedJumpingReadFile = returnAbsolutePath ($reverseComplementedJumpingReadFile);
$joiningEndPairs = $reverseComplementedJumpingReadFile;
$fishingEndPairs = $joiningEndPairs;
chdir ($outputDirectory);
$cmd = "$exeDir/outputMatedReadsAsReverseComplement.perl $jumpingLibraryReadFile > $joiningEndPairs";
runCommandAndExitIfBad ($cmd);

# $meanAndStdevJoinSeqLenByGapFile = "gap.insertMeanAndStdev.txt";
# $cmd = "$exeDir/getMeanAndStdevForGapsByGapNumUsingCeleraAsmFile.perl $CeleraTerminatorDirectory --contig-end-seq-file $joiningEndPairs > $meanAndStdevJoinSeqLenByGapFile";
# runCommandAndExitIfBad ($cmd);

$cmd = "echo \"cc $fauxInsertMean $fauxInsertStdev\" > meanAndStdevByPrefix.cc.txt";
runCommandAndExitIfBad ($cmd);

# Doing the section to avoid fishing using k-mers that occur (too) many times in the reads
$cmd = "jellyfish-2.0 count -s $jellyfishHashSize -C -t $numThreads -m $reduceReadSetKMerSize -L 100 -o restrictKmers.jf @readsFiles";
runCommandAndExitIfBad ($cmd);

$cmd = "jellyfish-2.0 dump -L $maxFishingKMerCount restrictKmers.jf -c > highCountKmers.txt";
runCommandAndExitIfBad ($cmd);

$localJellyfishHashSize = -s $fishingEndPairs;
$localJellyfishHashSize = int ($localJellyfishHashSize / .75);
$cmd = "jellyfish-2.0 count -s $localJellyfishHashSize -C -t $numThreads -m $reduceReadSetKMerSize -o fishingAll.jf $fishingEndPairs";
runCommandAndExitIfBad ($cmd);

open (FILE, "highCountKmers.txt");
while ($line = <FILE>) {
    ($kmer) = ($line =~ /^(\S+)\s/);
    $isHighCountKmer{$kmer} = 1; }
close (FILE);

# The following outputs a file containing all the k-mers occurring in the contig ends used for fishing that
# do not occur overly many times in the read database (as specified by our parameters)
$suffix = $localReadsFile . "_" . $reduceReadSetKMerSize . "_" . $kUnitigContinuationNumber;
$cmd = "jellyfish-2.0 dump fishingAll.jf -c |";
open (FILE, $cmd);
$outcmd = "| jellyfish-2.0 count -s $localJellyfishHashSize -C -t $numThreads -m $reduceReadSetKMerSize -o k_u_hash_${suffix}_faux_reads.jf /dev/fd/0";
open (OUTFILE, $outcmd);
print OUTFILE ">strangeFishingOutfile\n";
while ($line = <FILE>) {
    ($kmer) = ($line =~ /^(\S+)\s/);
    next if ($isHighCountKmer{$kmer});
    print OUTFILE $kmer,"N\n"; }
# print OUTFILE "\n";
close (FILE); close (OUTFILE);

$kUnitigFilename = "k_unitigs_${suffix}_faux_reads.fa";

# End section that eliminates k-mers that occur too often
$cmd = "jellyfish-2.0 dump -c k_u_hash_${suffix}_faux_reads.jf | awk 'BEGIN{n=0}{n++;print \">\"n\"length:$reduceReadSetKMerSize\\n\"\$1}' > $kUnitigFilename";
runCommandAndExitIfBad ($cmd);

$kUnitigFilesize = -s $kUnitigFilename;
$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -minreadsinsuperread 1 -kunitigsfile $kUnitigFilename -s $kUnitigFilesize -low-memory -l $reduceReadSetKMerSize --stopAfter findReadKUnitigMatches -t $numThreads -mkudisr 0 workFauxVsFaux $fishingEndPairs 1>>out.${suffix}_workFauxVsFaux 2>>out.${suffix}_workFauxVsFaux";
runCommandAndExitIfBad ($cmd);
$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -minreadsinsuperread 1 -kunitigsfile $kUnitigFilename -s $kUnitigFilesize -low-memory -l $reduceReadSetKMerSize --stopAfter findReadKUnitigMatches -t $numThreads -mkudisr 0 workReadsVsFaux @readsFiles 1>>out.${suffix}_workReadsVsFaux 2>>out.${suffix}_workReadsVsFaux";
runCommandAndExitIfBad ($cmd);
# We're still in the output directory
$cmd = "$exeDir/collectReadSequencesForLocalGapClosing --faux-reads-file $fishingEndPairs --faux-read-matches-to-kunis-file workFauxVsFaux/newTestOutput.nucmerLinesOnly --read-matches-to-kunis-file workReadsVsFaux/newTestOutput.nucmerLinesOnly --num-joins-per-directory $numJoinsPerDirectory";
for (@readsFiles) {
    $readFile = $_;
    $cmd .= " --reads-file $readFile"; }
$cmd .= " --max-reads-in-memory $maxReadsInMemory --dir-for-gaps .";
runCommandAndExitIfBad ($cmd);

# Now run the directories
# $cmd = "strace -o strace_out -e trace=process,open -f -ttt ";
$cmd = "$exeDir/runByDirectory -t $numThreads --jumping-read-joining-run $keepDirectoriesFlag --Celera-terminator-directory $CeleraTerminatorDirectory --max-nodes $maxNodes --min-kmer-len $minKMerLen --max-kmer-len $maxKMerLen --join-aggressive 1 --mean-for-faux-inserts $fauxInsertMean --stdev-for-faux-inserts $fauxInsertStdev --num-stdevs-allowed $numStdevsAllowed --output-dir $subdir2 --contig-end-sequence-file $joiningEndPairs --dir-for-read-sequences . -o output.txt --num-joins-per-directory $numJoinsPerDirectory";
runCommandAndExitIfBad ($cmd);

$cmd = "$exeDir/reportReadsToExclude.perl output.txt reverseComplemented.jumpingReadFile.fa > readsToExclude.txt";
runCommandAndExitIfBad ($cmd);

if (! $keepDirectoriesFlag) {
    $cmd = "\\rm -rf $subdir2"; print "$cmd\n"; system ($cmd); }

sub getJellyfishHashSizeNeeded
{
    my ($jellyfishHash) = @_;
    my ($cmd, $numRecs, $line, @flds, $jellyfishSizeNeeded);

    $cmd = "jellyfish-2.0 histo -t $numThreads -h 1 $jellyfishHash |";
    open (FILE, $cmd);
    $numRecs = 0;
    while ($line = <FILE>) {
	@flds = split (" ", $line);
	$numRecs += $flds[1]; }
    $jellyfishSizeNeeded = int ($numRecs/.7)+1;
    return ($jellyfishSizeNeeded);
}

sub getReadFileSize
{
    my (@files) = @_;
    my ($totSize, $file);
    $totSize = 0;
    for (@files) {
	$file = $_;
	$totSize += (-s $file); }
    return ($totSize);
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
    $reduceReadSetKMerSize = 21;
    $kUnitigContinuationNumber = 2;
    $maxKMerLen = 65;
    $minKMerLen = 17;
    $numThreads = 1;
    $maxFishingKMerCount = 300;
    $maxReadsInMemory = 100000000;
    $maxNodes = 1000;
    $fauxInsertMean = 350;
    $fauxInsertStdev = 50;
    $numStdevsAllowed = 3;
    $numJoinsPerDirectory = 100;
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
        if ($arg eq "--max-fishing-mer-count") {
            ++$i;
            $maxFishingKMerCount = $ARGV[$i];
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
	if ($arg eq "--jumping-library-read-file") {
	    ++$i;
	    $jumpingLibraryReadFile = $ARGV[$i];
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
	if ($arg eq "--maxnodes") {
	    ++$i;
	    $maxNodes = $ARGV[$i];
	    next; }
	if ($arg eq "--max-reads-in-memory") {
	    ++$i;
	    $maxReadsInMemory = $ARGV[$i];
	    next; }
	# The following is the max number of (Celera assembler) standard
	# deviations the joining sequence length can be from the length
	# predicted by the Celera assembler
	if ($arg eq "--num-stdevs-allowed") {
	    ++$i;
	    $numStdevsAllowed = $ARGV[$i];
	    next; }
        if ($arg eq "--faux-insert-mean") {
            ++$i;
            $fauxInsertMean = $ARGV[$i];
            next; }
        if ($arg eq "--faux-insert-stdev") {
            ++$i;
            $fauxInsertStdev = $ARGV[$i];
            next; }
        if ($arg eq "-s") {
            ++$i;
            $jellyfishHashSize = $ARGV[$i];
            next; }
	if ($arg eq "--keep-directories") {
	    $keepDirectoriesFlag = $arg;
	    next; }
	if ($arg eq "--num-joins-per-directory") {
	    ++$i;
	    $numJoinsPerDirectory = $ARGV[$i];
	    next; }
	if (-f $arg) {
	    push (@readsFiles, $arg);
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
    if (! $jellyfishHashSize) {
	print STDERR "You must enter a jellyfish hash size. Bye!\n";
	&reportUsage; }
    if (! $CeleraTerminatorDirectory) {
	print STDERR "You must enter a 9-terminator directory from a Celera run. Bye!\n";
	&reportUsage; }
    if (! $jumpingLibraryReadFile) {
	print STDERR "You must enter a jumping library read file to find bad mates. Bye!\n";
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

sub getReadMateName
{
    my ($rdName) = @_;
    my ($prefix, $num);
    ($prefix, $num) = ($rdName =~ /^(.+)(.)$/);
    if ($num % 2 == 0) {
        ++$num; }
    else {
        --$num; }
    $rdName = $prefix . $num;
    return ($rdName);
}

# This runs localCmd and captures the return code and makes
# sure it ran properly. Otherwise it exits.
sub runCommandAndExitIfBad
{
    my ($localCmd) = @_;
    my ($retCode, $exitValue);
    
    if ($localCmd =~ /\S/) {
        print STDERR "$localCmd\n";
	system ($localCmd);
        $retCode = $?;
        if ($retCode == -1) {
            print STDERR "failed to execute: $!\n";
	    goto failingRun; }
        elsif ($retCode & 127) {
            printf STDERR "child died with signal %d, %s coredump\n",
            ($retCode & 127), ($retCode & 128) ? 'with' : 'without';
	    goto failingRun; }
        else {
            $exitValue = $retCode >> 8;
            if ($exitValue == 255) { $exitValue = -1; }
            if ($exitValue != 0) {
                printf STDERR "child exited with value %d\n", $exitValue;
                print STDERR "Command \"$localCmd\" failed. Bye!\n";
		$retCode = $exitValue;
		goto failingRun; }
        }
    }
  successfulRun:
    return;
    
  failingRun:
    exit ($retCode);
}

