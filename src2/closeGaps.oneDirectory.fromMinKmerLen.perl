#!/usr/bin/env perl
use File::Basename;
$exeDir = dirname ($0);
$localReadsFile = "localReadsFile";
$joiningEndPairs = "fauxReads.fasta";
&processArgs;

if ($dirToChangeTo) {
    chdir ($dirToChangeTo); }

&runMainLoop;

sub runMainLoop
{
    my ($k, $suffix, $cmd, $directoryPassed);
    my ($fn1, $fn2, $isFirstLoop);
    my (@localLines, $localLine, $tempLine);
    my ($outFn, $sz, $sz2, $totInputSize, $minContinuation);

    $directoryPassed=0; 
    $isFirstLoop = 1;
    # Here's the main loop (k-mer values going up)
    for ($k = $minKMerLen; $k<=$maxKMerLen; ++$k) {
	$suffix = $localReadsFile . "_" . $k . "_" . $kUnitigContinuationNumber;
	if ($isFirstLoop) {
	    $isFirstLoop = 0;
	    $totInputSize = getReadFileSize (@readsFiles); }
#	$minContinuation = int ($k/2);
	$minContinuation = $k-1;
	$cmd = "$exeDir/create_k_unitigs_large_k -c $minContinuation -t $numThreads -m $k -n $totInputSize -l $k -f 0.000001 @readsFiles  |  grep -v '^>' | perl -ane '{\$seq=\$F[0]; \$F[0]=~tr/ACTGacgt/TGACtgac/;\$revseq=reverse(\$F[0]); \$h{(\$seq ge \$revseq)?\$seq:\$revseq}=1;}END{\$n=0;foreach \$k(keys \%h){print \">\",\$n++,\" length:\",length(\$k),\"\\n\$k\\n\"}}' >> k_unitigs_${suffix}.fa";
	if (runCommandAndReturnIfBad ($cmd)) {
	    last; }
	$cmd = "\\rm -rf out.$suffix"; system ($cmd);
	$cmd = "\\rm -rf work_$suffix"; system ($cmd);

	$tempFilename = "meanAndStdevByPrefix.cc.txt";
	if (! -e $tempFilename) {
	    $cmd = "echo cc $meanForFauxInserts $stdevForFauxInserts > $tempFilename"; 
	    runCommandAndExitIfBad ($cmd);
	}
	$tempKUnitigsFile = "k_unitigs_${suffix}.fa";
#	$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -num-stdevs-allowed $numStdevsAllowed -minreadsinsuperread 1 -kunitigsfile $tempKUnitigsFile -low-memory -l $k -t $numThreads -maxnodes $maxNodes -mkudisr 0 work_${suffix} $joiningEndPairs 1>>out.$suffix 2>>out.$suffix";
	$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -num-stdevs-allowed $numStdevsAllowed -closegaps -minreadsinsuperread 1 -kunitigsfile $tempKUnitigsFile -low-memory -l $k -t $numThreads -join-aggressive $joinAggressive -maxnodes $maxNodes -mkudisr 0 work_${suffix} $joiningEndPairs 1>>out.$suffix 2>>out.$suffix";
	if (runCommandAndReturnIfBad ($cmd)) {
	    last; }
	# Now looking for the report line from joinKUnitigs indicating what happened with the
	# faux mate pair (joined, missing nodes, too many nodes, etc.)
	
	open (FILE,"out.$suffix");
	$localLine = "";
        while($tempLine=<FILE>){
	    next if (not($tempLine =~ /^Num/));
            print $tempLine;
	    next unless ($tempLine =~ / 1\s*$/);
	    $localLine = $tempLine;
	    last; }
        close(FILE);

	# Initial check to see if the faux mates were joinable
        $lineForMatchHold="";
        open(FILE,"work_$suffix/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt");
        while ($line = <FILE>) {
            chomp ($line);
            @flds = split (" ", $line);
            ($val) = ($flds[0] =~ /^..(\d+)/);
            $val = int ($val/2);
            $lineForMatch = "$val $flds[1]";
            if ($lineForMatch eq $lineForMatchHold) {
                $directoryPassed = 1;
	    }
            $lineForMatchHold = $lineForMatch;
        }
        close(FILE);
	if ($localLine =~ /same unitig/){ $directoryPassed = 0; next; } 
	elsif ($directoryPassed == 1) { last; }
	elsif ($localLine =~ /not uniquely joinable/ || $localLine =~ /too many nodes/) { next; }
	elsif ($localLine =~ /missing sequence/) { last; }
	else { last; }
    }
    
    $passingKMerFile = "passingKMer.txt";
    open (FILE, ">$passingKMerFile");
    if ($directoryPassed) { 
        print FILE "$k\n"; }
    else {
	print FILE "11\n"; }
    close(FILE);
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
    return (0);
    
  failingRun:
    exit ($retCode);
}

sub runCommandAndReturnIfBad
{
    my ($localCmd) = @_;
    my ($retCode, $exitValue);
    
    if ($localCmd =~ /\S/) {
        print "$localCmd\n";
	$localCmd = "time $localCmd"; # For debugging 11/23/12
	system ($localCmd);
        $retCode = $?;
        if ($retCode == -1) {
            print STDERR "failed to execute: $!\n";
	    goto failingRun2; }
        elsif ($retCode & 127) {
            printf STDERR "child died with signal %d, %s coredump\n",
            ($retCode & 127), ($retCode & 128) ? 'with' : 'without';
	    goto failingRun2; }
        else {
            $exitValue = $retCode >> 8;
            if ($exitValue == 255) { $exitValue = -1; }
            if ($exitValue != 0) {
                printf STDERR "child exited with value %d\n", $exitValue;
                print STDERR "Command \"$localCmd\" failed. Bye!\n";
		$retCode = $exitValue;
		goto failingRun2; }
        }
    }
  successfulRun2:
    return (0);
    
  failingRun2:
    return ($retCode);
}

sub processArgs
{
    my ($arg, $cmd, @kmerLens, $i);
    $kUnitigContinuationNumber = 2;
    $maxKMerLen = 65;
    $minKMerLen = 17;
    $numThreads = 1;
    $useAllKUnitigs = 1;
    $maxNodes = 2000;
    $meanForFauxInserts = 500;
    $stdevForFauxInserts = 200;
    $numStdevsAllowed = 5;
    $joinAggressive = 0;
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
	if ($arg eq "--dir-for-kunitigs") {
	    ++$i;
	    $dirForKUnitigs = $ARGV[$i];
	    next; }
	if ($arg eq "--mean-for-faux-inserts") {
	    ++$i;
	    $meanForFauxInserts = $ARGV[$i];
	    next; }
	if ($arg eq "--stdev-for-faux-inserts") {
	    ++$i;
	    $stdevForFauxInserts = $ARGV[$i];
	    next; }
	if ($arg eq "--num-stdevs-allowed") {
	    ++$i;
	    $numStdevsAllowed = $ARGV[$i];
	    next; }
        if ($ARGV[$i] eq "--join-aggressive") {
            ++$i;
            $joinAggressive = $ARGV[$i];
            next; }
	if ($arg eq "--min-kmer-len") {
	    ++$i;
	    $minKMerLen = $ARGV[$i];
	    next; }
	if ($arg eq "--max-kmer-len") {
	    ++$i;
	    $maxKMerLen = $ARGV[$i];
	    next; }
	if ($arg eq "--reads-file") {
	    ++$i;
	    push (@readsFiles, $ARGV[$i]);
	    next; }
	if ($arg eq "--kunitig-continuation-number") {
	    ++$i;
	    $kUnitigContinuationNumber = $ARGV[$i];
	    next; }
	if ($arg eq "--maxnodes") {
	    ++$i;
	    $maxNodes = $ARGV[$i];
	    next; }
	if ($arg eq "--dir-to-change-to") {
	    ++$i;
	    $dirToChangeTo = $ARGV[$i];
	    next; }
	if (-f $arg) {
	    push (@readsFiles, $arg);
	    next; }
	if ($arg =~ /^\d+/) {
	    push (@kmerLens, $arg);
	    next; }
    }
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

