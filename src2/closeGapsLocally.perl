#!/usr/bin/env perl
#
# This program is used to close gaps in our assembly.
#
# Example invocation:
# closeGapsLocally.perl --min-kmer-len 17 --max-ker-len 31 --Celera-terminator-directory .../CA/9-terminator --num-threads 16 --output-directory outputDir --reads-file pe.cor.fa --reads-file sj.cor.fa
#
# There are no args, only flags (mostly) with arguments. They are as follows:
#
# Required flags:
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
# --contig-length-for-joining # : The length of sequence at the ends of the contigs
#                     which create the faux mate pairs which are joined (default: 100)
# --contig-length-for-fishing # : The length of sequence at the ends of the contigs
#                     to be used to find reads which might fit in the gaps (default: 100)
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
#
# The flags may be in any order.
use Cwd;
use File::Basename;
$exeDir = dirname ($0);
# Create absolute paths where necessary
# Must allow one to specify the dir with the hash and input k-unitigs file
# Must allow specification of the min k-unitig continuation values (default 2)
$cwd = cwd;
&processArgs;
$shm = "/dev/shm";
if ((! -d $shm) || $keepDirectoriesFlag) {
    $subdir2 = "subdir2"; }
else {
    $tempTime = time;
    $subdir2 = "$shm/$tempTime"; }
print "";
$CeleraTerminatorDirectory = returnAbsolutePath ($CeleraTerminatorDirectory);
for ($i=0; $i<=$#readsFiles; $i++) {
    $readsFile = $readsFiles[$i];
    $readsFile = returnAbsolutePath ($readsFile);
    $readsFiles[$i] = $readsFile; }

$localReadsFile = "localReadsFile";
print "outputDirectory = $outputDirectory\n";
if (! -e $outputDirectory) {
    $cmd = "mkdir $outputDirectory"; runCommandAndExitIfBad($cmd); }
$fishingEndPairs = "contig_end_pairs.${contigLengthForFishing}.fa";
$joiningEndPairs = "$outputDirectory/contig_end_pairs.${contigLengthForJoining}.fa";
$joiningEndPairs = returnAbsolutePath ($joiningEndPairs);
chdir ($outputDirectory);
$cmd = "$exeDir/getEndSequencesOfContigs.perl $CeleraTerminatorDirectory $contigLengthForJoining $contigLengthForFishing";
runCommandAndExitIfBad ($cmd);

$cmd = "$exeDir/create_end_pairs.perl $CeleraTerminatorDirectory $contigLengthForJoining > $joiningEndPairs";
runCommandAndExitIfBad ($cmd);

if ($contigLengthForJoining != $contigLengthForFishing) {
    $cmd = "$exeDir/create_end_pairs.perl $CeleraTerminatorDirectory $contigLengthForFishing > $fishingEndPairs";
    runCommandAndExitIfBad ($cmd);
}
$cmd = "echo \"cc $fauxInsertMean $fauxInsertStdev\" > meanAndStdevByPrefix.cc.txt";
runCommandAndExitIfBad ($cmd);

# Doing the section to avoid fishing using k-mers that occur (too) many times in the reads
$cmd = "jellyfish count -s $jellyfishHashSize -C -t $numThreads -m $reduceReadSetKMerSize -L 100 -o restrictKmers @readsFiles";
runCommandAndExitIfBad ($cmd);

# Must put in the merge command here
$localOutfile = "restrictKmers";
$testFile = "${localOutfile}_1";
if (-e $testFile) {
    $cmd = "jellyfish merge -o $localOutfile ${localOutfile}_*"; }
else {
    $cmd = "ln -s ${localOutfile}_0 $localOutfile"; }
runCommandAndExitIfBad ($cmd);

$cmd = "jellyfish dump -L $maxFishingKMerCount restrictKmers_0 -c > highCountKmers.txt";
runCommandAndExitIfBad ($cmd);

$localJellyfishHashSize = -s $fishingEndPairs;
$localJellyfishHashSize = int ($localJellyfishHashSize / .79) + 1;
$cmd = "jellyfish count -s $localJellyfishHashSize -C -t $numThreads -m $reduceReadSetKMerSize -o fishingAll $fishingEndPairs";
runCommandAndExitIfBad ($cmd);

open (FILE, "highCountKmers.txt");
while ($line = <FILE>) {
    ($kmer) = ($line =~ /^(\S+)\s/);
    $isHighCountKmer{$kmer} = 1; }
close (FILE);

# The following outputs a file containing all the k-mers occurring in the contig ends used for fishing that
# do not occur overly many times in the read database (as specified by our parameters)
$localJellyfishHashSize = 10000000;
$suffix = $localReadsFile . "_" . $reduceReadSetKMerSize . "_" . $kUnitigContinuationNumber;
$cmd = "jellyfish dump fishingAll_0 -c |";
open (FILE, $cmd);
$outcmd = "| jellyfish count -s $localJellyfishHashSize -C -r -t $numThreads -m $reduceReadSetKMerSize -o k_u_hash_${suffix}_faux_reads /dev/fd/0";
open (OUTFILE, $outcmd);
print OUTFILE ">strangeFishingOutfile\n";
while ($line = <FILE>) {
    ($kmer) = ($line =~ /^(\S+)\s/);
    next if ($isHighCountKmer{$kmer});
    print OUTFILE $kmer,"N\n"; }
# print OUTFILE "\n";
close (FILE); close (OUTFILE);
# End section that eliminates k-mers that occur too often

if (0) {
$suffix = $localReadsFile . "_" . $reduceReadSetKMerSize . "_" . $kUnitigContinuationNumber;
$localJellyfishHashSize = -s $fishingEndPairs;
$localJellyfishHashSize = int ($localJellyfishHashSize / .79) + 1;
$cmd = "jellyfish count -m $reduceReadSetKMerSize -t $numThreads -C -r -s $localJellyfishHashSize -o k_u_hash_${suffix}_faux_reads $fishingEndPairs";
runCommandAndExitIfBad ($cmd);
$tfile = "k_u_hash_${suffix}_faux_reads_1";
if (-e $tfile) {
    print STDERR "The jellyfish hash size must be made larger. Bye!\n";
    exit (1); }
} # End of 'if (0)' section

$cmd = "$exeDir/create_k_unitigs -C -t $numThreads -l $reduceReadSetKMerSize -m 1 -M 1 -o k_unitigs_${suffix}_faux_reads k_u_hash_${suffix}_faux_reads_0";
runCommandAndExitIfBad ($cmd);

$kUnitigFilename = "k_unitigs_${suffix}_faux_reads.fa";
$kUnitigFilesize = -s $kUnitigFilename;
$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -minreadsinsuperread 1 -kunitigsfile $kUnitigFilename -s $kUnitigFilesize -low-memory -l $reduceReadSetKMerSize --stopAfter findReadKUnitigMatches -t $numThreads -mkudisr 0 workFauxVsFaux $fishingEndPairs 1>>out.${suffix}_workFauxVsFaux 2>>out.${suffix}_workFauxVsFaux";
runCommandAndExitIfBad ($cmd);
$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -minreadsinsuperread 1 -kunitigsfile $kUnitigFilename -s $kUnitigFilesize -low-memory -l $reduceReadSetKMerSize --stopAfter findReadKUnitigMatches -t $numThreads -mkudisr 0 workReadsVsFaux @readsFiles 1>>out.${suffix}_workReadsVsFaux 2>>out.${suffix}_workReadsVsFaux";
runCommandAndExitIfBad ($cmd);
# We're still in the output directory
$cmd = "$exeDir/collectReadSequencesForLocalGapClosing --faux-reads-file $fishingEndPairs --faux-read-matches-to-kunis-file workFauxVsFaux/newTestOutput.nucmerLinesOnly --read-matches-to-kunis-file workReadsVsFaux/newTestOutput.nucmerLinesOnly  -t $numThreads $keepDirectoriesFlag --Celera-terminator-directory $CeleraTerminatorDirectory --max-nodes $maxNodes --min-kmer-len $minKMerLen --max-kmer-len $maxKMerLen --mean-for-faux-inserts $fauxInsertMean --stdev-for-faux-inserts $fauxInsertStdev --output-dir $subdir2 --contig-end-sequence-file $joiningEndPairs ";
for (@readsFiles) {
    $readFile = $_;
    $cmd .= " --reads-file $readFile"; }
# $cmd .= " --max-reads-in-memory $maxReadsInMemory --dir-for-gaps .";
$cmd .= " --dir-for-gaps . 2> out.err ;sort -gk1,1 out.err |perl -e '{open(OUTFILE1,\">superReadSequences.fasta\");open(OUTFILE2,\">readPlacementsInSuperReads.final.read.superRead.offset.ori.txt\");while(\$line=<STDIN>){chomp(\$line);\@f=split(/\\s+/,\$line);print OUTFILE1 \">\$f[0]\\n\$f[1]\\n\";print OUTFILE2 \"\$f[2] \$f[3] \$f[4] \$f[5]\\n\$f[6] \$f[7] \$f[8] \$f[9]\\n\";}}'";
runCommandAndExitIfBad ($cmd);

# Now run the directories
# $cmd = "$exeDir/runByDirectory -t $numThreads $keepDirectoriesFlag --Celera-terminator-directory $CeleraTerminatorDirectory --max-nodes $maxNodes --min-kmer-len $minKMerLen --max-kmer-len $maxKMerLen --mean-for-faux-inserts $fauxInsertMean --stdev-for-faux-inserts $fauxInsertStdev --output-dir $subdir2 --contig-end-sequence-file $joiningEndPairs --dir-for-read-sequences .";
# runCommandAndExitIfBad ($cmd);

$cmd = "$exeDir/getSequenceForLocallyClosedGaps.perl $CeleraTerminatorDirectory -contig-end-pairs-file $joiningEndPairs -working-directory .";
runCommandAndExitIfBad ($cmd);

&createMergedJoinFile;

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

sub getJellyfishHashSizeNeeded
{
    my ($jellyfishHash) = @_;
    my ($cmd, $numRecs, $line, @flds, $jellyfishSizeNeeded);

    $cmd = "jellyfish histo -t $numThreads -h 1 $jellyfishHash |";
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
    $maxFishingKMerCount = 1000;
#    $maxReadsInMemory = 100000000;
    $contigLengthForJoining = $contigLengthForFishing = 100;
    $maxNodes = 200000;
    $fauxInsertMean = 500;
    $fauxInsertStdev = 200;
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
	if ($arg eq "--maxnodes") {
	    ++$i;
	    $maxNodes = $ARGV[$i];
	    next; }
	if ($arg eq "--max-reads-in-memory") {
	    ++$i;
#	    $maxReadsInMemory = $ARGV[$i];
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

