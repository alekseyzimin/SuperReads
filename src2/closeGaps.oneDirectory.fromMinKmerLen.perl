#!/usr/bin/env perl
use File::Basename;
use File::Copy;
$exeDir = dirname ($0);
$localReadsFile = "localReadsFile";
$joiningEndPairs = "fauxReads.fasta";
&processArgs;

if ($dirToChangeTo) {
    chdir ($dirToChangeTo); }

$multipleJoinRun = &determineIfMultipleJoinsAreRun;
$joiningEndPairsOrig = $joiningEndPairs;

&runMainLoop;

sub runMainLoop
{
    my ($k, $suffix, $cmd, $directoryPassed);
    my ($fn1, $fn2);
    my (@localLines, $localLine, $tempLine);
    my ($outFn, $sz, $sz2, $totInputSize, $minContinuation);

    $directoryPassed=0; 
    # Here's the main loop (k-mer values going up)
    for ($k = $minKMerLen; $k<=$maxKMerLen; ++$k) {
	$suffix = $localReadsFile . "_" . $k . "_" . $kUnitigContinuationNumber;
	if (($k == $minKMerLen) || (($k == $minKMerLen+1) && $multipleJoinRun)) {
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
	if ($multipleJoinRun) {
	    $outputJoinFlag = "-output-join-result-for-each-join"; }
#	$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -num-stdevs-allowed $numStdevsAllowed -minreadsinsuperread 1 -kunitigsfile $tempKUnitigsFile -low-memory -l $k -t $numThreads -maxnodes $maxNodes -mkudisr 0 work_${suffix} $joiningEndPairs 1>>out.$suffix 2>>out.$suffix";
	$cmd = "$exeDir/createSuperReadsForDirectory.perl -mikedebug -noreduce -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.cc.txt -num-stdevs-allowed $numStdevsAllowed -closegaps -minreadsinsuperread 1 -kunitigsfile $tempKUnitigsFile -low-memory -l $k -t $numThreads -join-aggressive $joinAggressive $outputJoinFlag -maxnodes $maxNodes -mkudisr 0 work_${suffix} $joiningEndPairs 1>>out.$suffix 2>>out.$suffix";
	if (runCommandAndReturnIfBad ($cmd)) {
	    last; }
	# Now looking for the report line from joinKUnitigs indicating what happened with the
	# faux mate pair (joined, missing nodes, too many nodes, etc.)
	
	undef %type;
	@fwdReadsForNextPass = ();
	@goodFwdMates = ();
	open (FILE,"out.$suffix");
	if (! $multipleJoinRun) {
	    $localLine = "";
	    while ($tempLine = <FILE>) {
		next if ($tempLine !~ /^Num/);
		print $tempLine;
		next unless ($tempLine =~ / 1\s*$/);
		$localLine = $tempLine;
		last; }
	}
	else { # We must check what to print and if undef works with a hash
	    while ($tempLine = <FILE>) {
		chomp ($tempLine);
		next if ($tempLine =~ /^Num/);
		@flds = split (" ", $tempLine);
		if ($tempLine =~ /unresolved at end/) {
		    push (@goodFwdMates, $flds[0]);
		    $type{$flds[0]} = "unresolved"; }
		elsif ($tempLine =~ /unjoinable over max nodes/) {
		    push (@fwdReadsForNextPass, $flds[0]);
		    $type{$flds[0]} = "over max nodes"; }
		elsif ($tempLine =~ /unjoinable missing sequence/) {
		    if ($type{$flds[0]} !~ /\S/) {
			$type{$flds[0]} = "missing sequence"; } }
		elsif ($tempLine =~ /same unitig/) {
		    $type{$flds[0]} = "same unitig"; }
		elsif ($tempLine =~ /uniquely joinable/) {
		    push (@goodFwdMates, $flds[0]);
		    $type{$flds[0]} = "joinable"; }
		elsif ($tempLine =~ /disambiguated/) {
		    push (@goodFwdMates, $flds[0]);
		    $type{$flds[0]} = "disambiguated"; }
	    }
	}

	close(FILE);

	# Initial check to see if the faux mates were joinable
        $lineForMatchHold="";
        open(FILE,"work_$suffix/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt");
	if (! $multipleJoinRun) {
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
	}
	else {
	    undef %indexExists;
	    @goodJoins = ();
	    while ($line = <FILE>) {
		chomp ($line);
		@flds = split (" ", $line);
		($prefix, $val) = ($flds[0] =~ /^(..)(\d+)/);
		$val = int ($val/2);
		$lineForMatch = "$val $flds[1]";
		if ($indexExists{$lineForMatch}) {
		    $val2 = 2 * $val;
		    push (@goodJoins, "${prefix}$val2");
		    ++$val2;
		    push (@goodJoins, "${prefix}$val2"); }
		++$indexExists{$lineForMatch};
	    }
	}
        close(FILE);

	if (! $multipleJoinRun) {
	    if ($localLine =~ /same unitig/){ $directoryPassed = 0; next; } 
	    elsif ($directoryPassed == 1) { last; }
	    elsif ($localLine =~ /not uniquely joinable/ || $localLine =~ /too many nodes/) { next; }
	    elsif ($localLine =~ /missing sequence/) { last; }
	    else { last; }
	}
	else {
	    $passingReadsFile = "passingReadsFile.txt";
	    open (FILE, ">>$passingReadsFile");
	    for (@goodFwdMates) {
		$read = $_;
		print FILE "$read\n";
		$mateRead = &getMateFromReadName ($read);
		print FILE "$mateRead\n"; }
	    close (FILE);
	    $readsForNextPassFile = "readsForNextPass.txt";
	    open (FILE, ">$readsForNextPassFile");
	    undef %readIsNeededForNextPass;
	    for (@fwdReadsForNextPass) {
		$read = $_;
		print FILE "$read\n";
		$readIsNeededForNextPass{$read} = 1;
		$mateRead = &getMateFromReadName ($read);
		print FILE "$mateRead\n";
		$readIsNeededForNextPass{$mateRead} = 1; }
	    close (FILE);
	    $sz = -s $readsForNextPassFile;
	    if ($sz == 0) {
		unlink $readsForNextPassFile;
		last; }
	    else {
		open (FILE, $joiningEndPairs);
		$newJoiningEndPairs = $joiningEndPairsOrig . ".$k";
		open (OUTFILE, ">$newJoiningEndPairs");
		@readGroupsNeededForNextPass = ();
		while ($line = <FILE>) {
		    chomp ($line);
		    @flds = split (" ", $line);
		    ($readName) = ($flds[0] =~ /^.(.+)$/);
		    if ($readIsNeededForNextPass{$readName}) {
			$on = 1;
			if ($readName =~ /[02468]$/) {
			    ($tempFld) = ($readName =~ /^..(.+)$/);
			    $tempFld /= 2;
			    push (@readGroupsNeededForNextPass, $tempFld); }
		    }
		    else {
			$on = 0; }
		    print OUTFILE "$line\n" if ($on);
		    $line = <FILE>;
		    print OUTFILE $line if ($on); }
		close (FILE); close (OUTFILE);
		$joiningEndPairs = $newJoiningEndPairs;
	    }
	    if ($k == $minKMerLen) {
		$origReadFileHold = $readsFiles[0] . "Hold";
		move ($readsFiles[0], $origReadFileHold);
		$spclCmd = "grep -A 1 -P \" \(" . $readGroupsNeededForNextPass[0];
		for ($l=1; $l<=$readGroupsNeededForNextPass; ++$l) {
		    $spclCmd .= ("|" . $readGroupsNeededForNextPass[$l]); }
		$spclCmd .= "\) \" $origReadFileHold > $readsFiles[0]";
		system ($spclCmd); }
	} # Ends section if it is a multiple-join run
    }
    
    $passingKMerFile = "passingKMer.txt";
    open (FILE, ">$passingKMerFile");
    if ($directoryPassed) { 
        print FILE "$k\n"; }
    else {
	print FILE "11\n"; }
    close(FILE);
    $passingReadsFile = "passingReadsFile.txt";
    $passingReadsFileHold = "passingReadsFile.orig.txt";
    $passingReadsFileSize = -s $passingReadsFile;
    if ($passingReadsFileSize > 0) {
	undef %isPassing;
	move ($passingReadsFile, $passingReadsFileHold);
	open (FILE, $passingReadsFileHold);
	while ($line = <FILE>) {
	    chomp ($line);
	    $isPassing{$line} = 1; }
	close (FILE);
	open (OUTFILE, ">$passingReadsFile");
	open (FILE, $joiningEndPairs);
	while ($line = <FILE>) {
	    chomp ($line);
	    @flds = split (" ", $line);
	    ($name) = ($flds[0] =~ /^.(.+)$/);
	    if ($isPassing{$name}) {
		print OUTFILE "$flds[1]\n"; }
	}
	close (FILE);
    }
}

sub getMateFromReadName
{
    my ($readName) = @_;
    my ($prefix, $last, $mateName);
    ($prefix, $last) = ($readName =~ /^(.+)(.)$/);
    if ($last % 2 == 0) {
	++$last; }
    else {
	--$last; }
    $mateName = $prefix . $last;
    return ($mateName);
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
    $maxNodes = 2000;
    $meanForFauxInserts = 500;
    $stdevForFauxInserts = 200;
    $numStdevsAllowed = 5;
    $joinAggressive = 0;
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
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

sub determineIfMultipleJoinsAreRun
{
    my ($cmd, $tline, @flds);
    $cmd = "wc -l $joiningEndPairs";
    open (FILE, "$cmd |"); $tline = <FILE>; close (FILE);
    @flds = split (" ", $tline);
    if ($flds[0] > 4) {
	return (1); }
    else {
	return (0); }
}

