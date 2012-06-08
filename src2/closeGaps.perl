#!/usr/bin/env perl
#
# This program is used to close gaps in our assembly.
#
# Example invocation:
# closeGaps.perl --min-kmer-len 15 --max-ker-len 31 --Celera-terminator-directory .../CA/9-terminator --jellyfish-hash-size 2000000000 --num-threads 16 --output-directory outputDir --reads-file pe.cor.fa --reads-file sj.cor.fa
#
# There are no args, only flags (mostly) with arguments. They are as follows:
#
# Required flags:
# --Celera-terminator-directory dir : specify the Celera terminator directory
#                           where the assembly whose gaps must be closed exists
# --reads-file filename : specify a read file to use (multiple files allowed,
#                            so long as the flag is repeated)
# --output-directory dir : specify the output directory
# --jellyfish-hash-size # : specify the jellyfish hash size
#
# Flags for required values which have defaults (i.e. flag not necessary)
# --min-kmer-len # : specify the min kmer len used (default: 15)
# --max-kmer-len # : specify the max kmer len used (default: 31)
# --num-threads # : specify the number of threads (default: 1)
# -t # : same as --num-threads #
# --contig-length-for-joining # : The length of sequence at the ends of the contigs
#                     which create the faux mate pairs which are joined (default: 100)
# --use-all-kunitigs : Use k-unitigs which are the k-mer length as well as all those longer than
#                     the k-mer length. (The default is not to use k-unis of the k-mer length)
# --maxnodes # : The maximum number of nodes allowed when trying to join the
#                     faux reads (default: 2000)
# --kunitig-continuation-number # : specify the number to continue when running
#                     create_k_unitigs (the -m and -M options to that program)
#                     (currently "invalidated") (default: 2)
#
# Flags for optional aspects
# --dir-for-kunitigs dir : specifies the directory to get k-unitigs from
#                           if we have them
# --reduce-read-set # : Start by reducing the read set to only those that
#                     match a k-unitig from the genomic sequences surrounding
#                     a gap. The number specifies the k-mer size used to
#                     find these matches. (Don't make it too small.)
# --contig-length-for-fishing # : The length of sequence at the ends of the contigs
#                     to be used to find reads which might fit in the gaps (default: 100)
# --noclean : Don't clean up after the run
#
# The flags may be in any order.
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
if (! $reduceReadSetKMerSize) {
    $contigLengthForFishing = $contigLengthForJoining; }
print "";
$CeleraTerminatorDirectory = returnAbsolutePath ($CeleraTerminatorDirectory);
for ($i=0; $i<=$#readsFiles; $i++) {
    $readsFile = $readsFiles[$i];
    $readsFile = returnAbsolutePath ($readsFile);
    $readsFiles[$i] = $readsFile; }

if ($dirForKUnitigs) {
    $dirForKUnitigs = returnAbsolutePath ($dirForKUnitigs); }
$localReadsFile = "localReadsFile";
if (! -e $outputDirectory) {
    $cmd = "mkdir $outputDirectory"; runCommandAndExitIfBad($cmd); }
chdir ($outputDirectory);
$fishingEndPairs = "contig_end_pairs.${contigLengthForFishing}.fa";
$joiningEndPairs = "contig_end_pairs.${contigLengthForJoining}.fa";
$cmd = "$exeDir/getEndSequencesOfContigs.perl $CeleraTerminatorDirectory $contigLengthForJoining $contigLengthForFishing";
runCommandAndExitIfBad ($cmd);

$cmd = "$exeDir/create_end_pairs.perl $CeleraTerminatorDirectory $contigLengthForJoining > $joiningEndPairs";
runCommandAndExitIfBad ($cmd);

if ($contigLengthForJoining != $contigLengthForFishing) {
    $cmd = "$exeDir/create_end_pairs.perl $CeleraTerminatorDirectory $contigLengthForFishing > $fishingEndPairs";
    runCommandAndExitIfBad ($cmd);
}
$cmd = "echo \"cc $meanForGap $stdevForGap\" > meanAndStdevByPrefix.cc.txt";
runCommandAndExitIfBad ($cmd);

if ($reduceReadSetKMerSize) {
    $suffix = $localReadsFile . "_" . $reduceReadSetKMerSize . "_" . $kUnitigContinuationNumber;
    $localJellyfishHashSize = -s $fishingEndPairs;
    $localJellyfishHashSize = int ($localJellyfishHashSize / .79) + 1;
    $cmd = "jellyfish count -m $reduceReadSetKMerSize -t $numThreads -C -r -s $localJellyfishHashSize -o k_u_hash_${suffix}_all_faux_reads $fishingEndPairs";
    runCommandAndExitIfBad ($cmd);
    $cmd = "jellyfish dump -U $maxFishingKMerCount -c k_u_hash_${suffix}_all_faux_reads_0 |perl -ane 'BEGIN{print \">0\\n\"}{for(\$i=0;\$i<\$F[1];\$i++){print \$F[0],\"N\"}}' | jellyfish count -t $numThreads -p 126 -C -r -o k_u_hash_${suffix}_faux_reads -s $localJellyfishHashSize -m $reduceReadSetKMerSize /dev/fd/0";
    runCommandAndExitIfBad ($cmd);
    $tfile = "k_u_hash_${suffix}_faux_reads_1";
    if (-e $tfile) {
	print STDERR "The jellyfish hash size must be made larger. Bye!\n";
	exit (1); }
    $cmd = "$tempExeDir/create_k_unitigs -C -t $numThreads -l $reduceReadSetKMerSize -o k_unitigs_${suffix}_faux_reads -m 1 -M 1 k_u_hash_${suffix}_faux_reads_0";
    runCommandAndExitIfBad ($cmd);

    $cmd = "$tempExeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -minreadsinsuperread 1 -kunitigsfile k_unitigs_${suffix}_faux_reads.fa -low-memory -l $reduceReadSetKMerSize --stopAfter joinKUnitigs -t $numThreads -mkudisr 0 work_${suffix}_faux_reads @readsFiles 1>>out.${suffix}_faux_reads 2>>out.${suffix}_faux_reads";
    runCommandAndExitIfBad ($cmd);

    $tfile = "work_${suffix}_faux_reads/newTestOutput.nucmerLinesOnly";
    open (FILE, $tfile);
    while ($line = <FILE>) {
	chomp ($line);
	@flds = split (" ", $line);
	$readName = $flds[0];
	$isNeeded{$readName} = 1;
	$readName = getReadMateName($readName);
	$isNeeded{$readName} = 1;
    }
    close (FILE);

    &createReducedReadSet;
}

&runMainLoop;

$cmd = "$exeDir/getSequenceForClosedGaps.perl $CeleraTerminatorDirectory $joiningEndPairs -reads-file $localReadsFile $minKMerLen $maxKMerLen";
runCommandAndExitIfBad ($cmd);

&createMergedJoinFile;

if (! $noClean) { &cleanUp; }

sub cleanUp
{
    my ($k, $suffix, $cmd);
    for ($k=$maxKMerLen; $k>=$minKMerLen; $k--) {
	$suffix = $localReadsFile . "_" . $k . "_" . $kUnitigContinuationNumber;
	$cmd = "\\rm -rf work_$suffix"; system ($cmd);
	$cmd = "\\rm k_u_hash_${suffix}_0"; system ($cmd);
	$cmd = "\\rm k_u_hash_${suffix}_all_faux_reads"; system ($cmd);
	$cmd = "\\rm  out.$suffix"; system ($cmd);
	$cmd = "\\rm joined.$suffix"; system ($cmd);
    }
}

sub createReducedReadSet
{
    my ($outputReadsFile, $line, $type, $readName, $isNeeded, $seq, $state, $qual, $cmd);

    $cmd = "cat @readsFiles |";
    open (FILE, $cmd);
    $outputReadsFile = "reducedReadsFile.fa";
    $outputReadsFile = returnAbsolutePath ($outputReadsFile);
    @readsFiles = ($outputReadsFile);
    open (OUTFILE, ">$outputReadsFile");
    while ($line=<FILE>) {
	if ($line =~ /^>/) { # then fasta
	    $type = "fasta";
	    ($readName) = ($line =~ /^.(\S+)/);
	    $isNeeded = $isNeeded{$readName};
	    print OUTFILE ">$readName\n" if($isNeeded);
	    next;
	}
	if (($line !~ /^\@/) && ($type eq "fasta")) {
	    print OUTFILE $line if ($isNeeded);
	    next; }
	if ($line =~ /^\@/){ # then fastq
	    $type = "fastq";
	    ($readName) = ($line =~ /^.(\S+)/);
	    $isNeeded = $isNeeded{$readName};
	    print OUTFILE ">$readName\n" if($isNeeded);
	    $seq = $qual = "";
	    $state = "seq";
	    next;
	}
	if (($state eq "seq") && ($line =~ /^\+/)) {
	    $state = "qual";
	    next; }
	if ($state eq "seq") {
	    chomp ($line);
	    $seq .= $line;
	    next; }
	if ($state eq "qual") {
	    chomp ($line);
	    $qual .= $line;
	    if (length($seq) == length ($qual)) {
		print OUTFILE "$seq\n" if ($isNeeded); }
	    next;
	}
	print STDERR "We should never get here. Error on input read '${readName}'!\n";
	exit (1);
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
    my ($k, $kMinus1, $suffix, $cmd, $tfile, $minKUniLengthForPass);
    
    for ($k=$maxKMerLen; $k>=$minKMerLen; $k--) {
	$kMinus1 = $k-1;
	$suffix = $localReadsFile . "_" . $k . "_" . $kUnitigContinuationNumber;
	if ($dirForKUnitigs) {
	    $cmd = "ln -s $dirForKUnitigs/k_u_hash_${suffix}_0 .";
	    runCommandAndExitIfBad ($cmd);
	    
	    $cmd = "ln -s $dirForKUnitigs/k_unitigs_${suffix}.fa .";
	    runCommandAndExitIfBad ($cmd); }
	else {
	    if ($k == $maxKMerLen) {
		$totInputSize = getReadFileSize (@readsFiles);
		$maxJellyfishHashSize = int ($totInputSize / .8)+1;
		if ($maxJellyfishHashSize < $jellyfishHashSize) {
		    $jellyfishHashSize = $maxJellyfishHashSize; }
	    }
	    $cmd = "jellyfish count -m $k -t $numThreads -C -r -s $jellyfishHashSize -o k_u_hash_${suffix} @readsFiles";
	    runCommandAndExitIfBad ($cmd);
	    $tfile = "k_u_hash_${suffix}_1";
	    if (-e $tfile) {
		print STDERR "The jellyfish hash size must be made larger. Bye!\n";
		exit (1); }
	    if ($k == $maxKMerLen) {
		$jellyfishHashSize = getJellyfishHashSizeNeeded ("k_u_hash_${suffix}_0"); }
	    if ($useAllKUnitigs == 1) {
		$minKUniLengthForPass = $k; }
	    else {
		$minKUniLengthForPass = $k+1; }
	    $cmd = "$tempExeDir/create_k_unitigs --cont-on-low  --low-stretch=$kMinus1 -C -t $numThreads -l $minKUniLengthForPass -o k_unitigs_${suffix} -m $kUnitigContinuationNumber -M $kUnitigContinuationNumber k_u_hash_${suffix}_0";
	    runCommandAndExitIfBad ($cmd);
	}
	$cmd = "\\rm -rf out.$suffix"; system ($cmd);
	$cmd = "\\rm -rf work_$suffix"; system ($cmd);
	
	$cmd = "$tempExeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -minreadsinsuperread 1 -kunitigsfile k_unitigs_${suffix}.fa -low-memory -l $k -t $numThreads -maxnodes $maxNodes -mkudisr 0 work_${suffix} $joiningEndPairs 1>>out.$suffix 2>>out.$suffix";
	runCommandAndExitIfBad ($cmd);

	$cmd = "grep '^Num ' out.$suffix"; system ($cmd);

	$cmd = "cat work_$suffix/readPositionsInSuperReads | $exeDir/outputJoinedPairs.perl > joined.$suffix";
	runCommandAndExitIfBad ($cmd);
    }
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
    $kUnitigContinuationNumber = 2;
    $maxKMerLen = 31;
    $minKMerLen = 15;
    $numThreads = 1;
    $maxFishingKMerCount = 5;
    $useAllKUnitigs = 0;
    $noClean = 0;
    $contigLengthForJoining = $contigLengthForFishing = 100;
    $maxNodes = 2000;
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
	if ($arg eq "--dir-for-kunitigs") {
	    ++$i;
	    $dirForKUnitigs = $ARGV[$i];
	    next; }
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
	if ($arg eq "--kunitig-continuation-number") {
	    ++$i;
	    $kUnitigContinuationNumber = $ARGV[$i];
	    next; }
	if ($arg eq "--output-directory") {
	    ++$i;
	    $outputDirectory = $ARGV[$i];
	    next; }
	if ($arg eq "--jellyfish-hash-size") {
	    ++$i;
	    $jellyfishHashSize = $ARGV[$i];
	    next; }
	if (($arg eq "--num-threads") || ($arg eq "-t")) {
	    ++$i;
	    $numThreads = $ARGV[$i];
	    next; }
	if ($arg eq "--reduce-read-set") {
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
	if ($arg eq "--use-all-kunitigs") {
	    $useAllKUnitigs = 1;
	    next; }
	if ($arg eq "--noclean") {
	    $noClean = 1;
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
    if ($#kmerLens >= 2) {
	$kUnitigContinuationNumber = $kmerLens[2]; }
    if (! $CeleraTerminatorDirectory) {
	print STDERR "You must enter a 9-terminator directory from a Celera run. Bye!\n";
	&reportUsage; }
    if (! $jellyfishHashSize) {
	print STDERR "You must enter a jellyfish hash size. Bye!\n";
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
        print "$localCmd\n";
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

