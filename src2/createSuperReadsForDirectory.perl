#!/usr/bin/perl
#
# This exec takes a (set of) input read files (in fasta format) and a
# file of input k-unitigs (specified with the switch -kunitigsfile) and outputs
# the set of super-reads for these reads (in fasta format).
#
# The args are the input fasta files as well as (optionally) the directory
# where you want the work to occur. If the directory is not specified, the
# work is done in the current directory.
# If the work directory doesn't exist then it is created.
# 
# The flags are as follows:
# -l merLen : the length of the k-mer to use for the calculations (31)
# -s tableSize : the size of the table when running jellyfish (2,000,000,000)
# -t numProcessors : the number of processors to run jellyfish and create_k_unitigs (16)
# -kunitigsfile filename : a user-given k-unitigs file; otherwise we calculate
# -mean-and-stdev-by-prefix-file filename : a file giving mate info about each
#                      library. Each line is the 2-letter prefix for the reads
#                      in the library followed by its mean and stdev. This
#                      file is mandatory unless -jumplibraryreads is specified
# -mkudisr numBaseDiffs : max base diffs between overlapping k-unitigs in super-reads (0)
# -minreadsinsuperread minReads : super-reads containing fewer than numReads
#                                reads will be eliminated (2)
# -merged-unitig-data-prefix prefix : the prefix for the filenames relating to
#                      merged unitig data. We assume that the k-unitig sequence
#                      is in  'prefix'.fasta, and the map file from orig to
#                      merged k-unitigs is in 'prefix'.map.
# --stopAfter target : Stop the run after one of the following "target" names:
#               createLengthStatisticsFiles
#               createKUnitigHashTable
#               addMissingMates
#               findReadKUnitigMatches
#               createLengthStatisticsForMergedKUnitigsFiles
#               createKUnitigMaxOverlaps
#               joinKUnitigs
#               getSuperReadInsertCounts
#               createFastaSuperReadSequences
#               reduceSuperReads
#               createFinalReadPlacementFile
#               createFinalSuperReadFastaSequences
# -noclean : don't clean up the files afterwards
# -mikedebug : don't kill off intermediate results
# -jumplibraryreads : we are generating for jump-library reads; a k-unitigs
#                                 file must be specified
# -h : help 
use File::Basename;
use Cwd;
$exeDir = dirname ($0);
$pwd = cwd;
if ($exeDir !~ /^\//) {
    $exeDir = "$pwd/$exeDir"; }

$maxNodes=2000;
$noReduce=0;
&processArgs;
if ($jumpLibraryReads) {
    $maxNodes = 0; }
$merLenMinus1 = $merLen - 1;
$maxHashFillFactor = .8;

$successFile = "$workingDirectory/superReads.success";
unlink ($successFile) if (-e $successFile);
# The following is set to 1 when the first success file for a step is missing
$mustRun = 0;
if (! -d $workingDirectory) {
    $cmd = "mkdir $workingDirectory";
    print "$cmd\n"; system ($cmd); }
# We now require that a k-unitigs file was passed on the command line
if ($kUnitigsFile !~ /^\//) {
    $kUnitigsFile = "$pwd/$kUnitigsFile"; }
$jellyfishKUnitigDataPrefix = "$workingDirectory/organismMerCountsForKUnitigs";
$jellyfishKUnitigHashFile = $jellyfishKUnitigDataPrefix . "_0";
$kUnitigLengthsFile = "$workingDirectory/kUnitigLengths.txt";
# The following stores the actual number of k-unitigs
$numKUnitigsFile = "$workingDirectory/numKUnitigs.txt";
# The following stores the largest k-unitig number (+1)
$maxKUnitigNumberFile = "$workingDirectory/maxKUnitigNumber.txt";
$totBasesInKUnitigsFile = "$workingDirectory/totBasesInKUnitigs.txt";
if ($mergedUnitigDataPrefix) {
    $mergedUnitigInputKUnitigsFile = $mergedUnitigDataPrefix . ".fasta";
    $mergedUnitigInputKUnitigMappingFile = $mergedUnitigDataPrefix . ".map";
    &runCommandAndExitIfBad ("", $mergedUnitigInputKUnitigsFile, 1, "mergedKUnitigFastaFileExists");
    &runCommandAndExitIfBad ("", $mergedUnitigInputKUnitigMappingFile, 1, "mergedKUnitigMapFileExists");
    $mergedKUnitigLengthsFile = "$workingDirectory/mergedKUnitigs.kUnitigLengths.txt";
    $mergedNumKUnitigsFile = "$workingDirectory/mergedKUnitigs.numKUnitigs.txt";
    $mergedMaxKUnitigNumberFile = "$workingDirectory/mergedKUnitigs.maxKUnitigNumber.txt";
    $mergedTotBasesInKUnitigsFile = "$workingDirectory/mergedKUnitigs.totBasesInKUnitigs.txt"; }
else {
    $mergedUnitigInputKUnitigsFile = $kUnitigsFile;
    $mergedKUnitigLengthsFile = $kUnitigLengthsFile;
    $mergedNumKUnitigsFile = $numKUnitigsFile;
    $mergedMaxKUnitigNumberFile = $maxKUnitigNumberFile;
    $mergedTotBasesInKUnitigsFile = $totBasesInKUnitigsFile; }
    
$totReadFile = "$workingDirectory/inputReads.fasta";
if ($#fastaFiles == 0) {
    $totReadFile = $fastaFiles[0]; }
$readsAfterAddingMissingMates = "$workingDirectory/inputreads.fa";
$numReadsFile = "$workingDirectory/numReads.txt";
$prefixForOverlapsBetweenKUnitigs = "$workingDirectory/overlap";
$kUnitigOverlapsFile = "${prefixForOverlapsBetweenKUnitigs}.overlaps";
$superReadCountsFile = "$workingDirectory/superReadCounts.all";

$joinerOutput = "$workingDirectory/readPositionsInSuperReads";
$readKUnitigMatchOutput = "$workingDirectory/newTestOutput.nucmerLinesOnly";
$sequenceCreationErrorFile = "$workingDirectory/createFastaSuperReadSequences.errors.txt";
# $myProgOutput2 = "$workingDirectory/readPlacementsInSuperReads.postMateMerge.read.superRead.offset.ori.txt";
# $myProgOutput3 = "$workingDirectory/superReadCounts.count.superRead.txt";
$finalSuperReadSequenceFile = "$workingDirectory/superReadSequences.fasta";
$finalReadPlacementFile = "$workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt";
# $reducedReadPlacementFile = "$workingDirectory/readPlacementsInSuperReads.reduced.read.superRead.offset.ori.txt";

if ($jumpLibraryReads) {
    $finalSuperReadSequenceFile = "$workingDirectory/superReadSequences.jumpLibrary.fasta";
}

# goto startHere;


&cleanUpFailureDirectories;

# In addition to obvious output file, this also generates the files
# numKUnitigs.txt, maxKUnitigNumber.txt, and totBasesInKUnitigs.txt in
# $workingDirectory
$cmd = "cat $kUnitigsFile | $exeDir/getLengthStatisticsForKUnitigsFile.perl $workingDirectory > $kUnitigLengthsFile";
&runCommandAndExitIfBad ($cmd, $kUnitigLengthsFile, 1, "createLengthStatisticsFiles", $totBasesInKUnitigsFile, $numKUnitigsFile, $maxKUnitigNumberFile, $kUnitigLengthsFile);

$minSizeNeededForTable = &reportMinJellyfishTableSizeForKUnitigs;
 redoKUnitigsJellyfish:
    $cmd = "jellyfish count -m $merLen -r -o $jellyfishKUnitigDataPrefix -c 6 -p 126 --both-strands -s $minSizeNeededForTable -t $numProcessors $kUnitigsFile";
&runCommandAndExitIfBad ($cmd, $jellyfishKUnitigHashFile, 1, "createKUnitigHashTable", "$workingDirectory/organismMerCountsForKUnitigs_0");

$tableResizeFactor = &returnTableResizeAmount ($jellyfishKUnitigDataPrefix, $jellyfishKUnitigHashFile);
if ($tableResizeFactor > 1) {
    $tableSize *= 2;
    print "Resizing the table to $tableSize for the k-unitig jellyfish run\n";
    goto redoKUnitigsJellyfish; }

# In addition to obvious output file, this also generates the files
# mergedKUnitigs.numKUnitigs.txt, mergedKUnitigs.maxKUnitigNumber.txt, and mergedKUnitigs.totBasesInKUnitigs.txt in
# $workingDirectory
if ($mergedUnitigDataPrefix) {
    $cmd = "cat $mergedUnitigInputKUnitigsFile | $exeDir/getLengthStatisticsForKUnitigsFile.perl -output-prefix mergedKUnitigs $workingDirectory > $mergedKUnitigLengthsFile";
    &runCommandAndExitIfBad ($cmd, $mergedKUnitigLengthsFile, 1, "createLengthStatisticsForMergedKUnitigsFiles", $mergedTotBasesInKUnitigsFile, $mergedNumKUnitigsFile, $mergedMaxKUnitigNumberFile, $mergedKUnitigLengthsFile); }
else { # The following is so we stop here if we are using --stopAfter createLengthStatisticsForMergedKUnitigsFiles but don't use the merged k-unitigs
    &runCommandAndExitIfBad ("", "", 0, "createLengthStatisticsForMergedKUnitigsFiles"); }

open (FILE, $mergedMaxKUnitigNumberFile); $maxKUnitigNumber = <FILE>; chomp ($maxKUnitigNumber); close (FILE);
$cmd = "$exeDir/createKUnitigMaxOverlaps $mergedUnitigInputKUnitigsFile -kmervalue $merLen -largest-kunitig-number ".(int($maxKUnitigNumber)+1)." $prefixForOverlapsBetweenKUnitigs";
&runCommandAndExitIfBad($cmd, $kUnitigOverlapsFile, 1, "createKUnitigMaxOverlaps", $kUnitigOverlapsFile, "$workingDirectory/overlap.coords");

$cmd = "ln -s $totReadFile $readsAfterAddingMissingMates";
&runCommandAndExitIfBad ($cmd, $readsAfterAddingMissingMates, 1, "addMissingMates", $readsAfterAddingMissingMates);

$cmd = "$exeDir/findMatchesBetweenKUnitigsAndReads $jellyfishKUnitigHashFile -t $numProcessors -o $readKUnitigMatchOutput $kUnitigsFile $maxKUnitigNumberFile  $readsAfterAddingMissingMates";
&runCommandAndExitIfBad ($cmd, $readKUnitigMatchOutput, 1, "findReadKUnitigMatches", $readKUnitigMatchOutput);
if (! $mikedebug) { &killFiles ($jellyfishKUnitigHashFile, $readsAfterAddingMissingMates); }

# Do the shooting method here
if ($mergedUnitigDataPrefix) {
    $mergedUnitigDataFileStr = "--kunitigs-translation-file $mergedUnitigInputKUnitigMappingFile"; }
$cmd = "$exeDir/joinKUnitigs_v3 --max-nodes-allowed $maxNodes --mean-and-stdev-by-prefix-file $meanAndStdevByPrefixFile --unitig-lengths-file $mergedKUnitigLengthsFile --num-kunitigs-file $mergedMaxKUnitigNumberFile --overlaps-file $kUnitigOverlapsFile --min-overlap-length $merLenMinus1 -o $joinerOutput $mergedUnitigDataFileStr -t $numProcessors $readKUnitigMatchOutput";
&runCommandAndExitIfBad ($cmd, $joinerOutput, 1, "joinKUnitigs", $joinerOutput);

if ($jumpLibraryReads) {
    goto jumpLibraryCalculations; }

$cmd= "$exeDir/getSuperReadInsertCountsFromReadPlacementFileTwoPasses -n `cat $numKUnitigsFile | awk '{print \$1*100}'` -o $superReadCountsFile $joinerOutput";
&runCommandAndExitIfBad ($cmd, $superReadCountsFile, 1, "getSuperReadInsertCounts", $superReadCountsFile);

if ($mergedUnitigDataPrefix) {
    $mergedUnitigDataFileStr = "-maxunitignumberfile $mergedMaxKUnitigNumberFile"; }

$goodSuperReadsNamesFile = "$workingDirectory/superReadNames.txt";
$fastaSuperReadErrorsFile = "$workingDirectory/createFastaSuperReadSequences.errors.txt";

if($noReduce==0) {
    $localGoodSequenceOutputFile = "${finalSuperReadSequenceFile}.all";
    $superReadNameAndLengthsFile = "$workingDirectory/sr_sizes.tmp";
    $reduceFile = "$workingDirectory/reduce.tmp";
    $cmd = "cat $superReadCountsFile | $exeDir/createFastaSuperReadSequences $workingDirectory /dev/fd/0 -seqdiffmax $seqDiffMax -min-ovl-len $merLenMinus1 -minreadsinsuperread $minReadsInSuperRead $mergedUnitigDataFileStr -good-sr-filename $goodSuperReadsNamesFile -kunitigsfile $mergedUnitigInputKUnitigsFile -good-sequence-output-file $localGoodSequenceOutputFile -super-read-name-and-lengths-file $superReadNameAndLengthsFile 2> $sequenceCreationErrorFile";
    &runCommandAndExitIfBad ($cmd, $superReadNameAndLengthsFile, 1, "createFastaSuperReadSequences", $localGoodSequenceOutputFile, $goodSuperReadsNamesFile, $superReadNameAndLengthsFile);

    $cmd = "$exeDir/reduce_sr $maxKUnitigNumber $mergedKUnitigLengthsFile $merLen $superReadNameAndLengthsFile -o $reduceFile";
    &runCommandAndExitIfBad ($cmd, $reduceFile, 1, "reduceSuperReads", $reduceFile, $fastaSuperReadErrorsFile);

    $cmd = "$exeDir/eliminateBadSuperReadsUsingList $joinerOutput $goodSuperReadsNamesFile -super-read-reduction-file $reduceFile > $finalReadPlacementFile";
    &runCommandAndExitIfBad ($cmd, $finalReadPlacementFile, 1, "createFinalReadPlacementFile", $finalReadPlacementFile);

    $cmd = "$exeDir/outputRecordsNotOnList $reduceFile $localGoodSequenceOutputFile 0 --fld-num 0 > $finalSuperReadSequenceFile";
    &runCommandAndExitIfBad ($cmd, $finalSuperReadSequenceFile, 1, "createFinalSuperReadFastaSequences", $finalSuperReadSequenceFile); }
else {
    $cmd = "cat $superReadCountsFile | $exeDir/createFastaSuperReadSequences $workingDirectory /dev/fd/0 -seqdiffmax $seqDiffMax -min-ovl-len $merLenMinus1 -minreadsinsuperread $minReadsInSuperRead $mergedUnitigDataFileStr -good-sr-filename $goodSuperReadsNamesFile -kunitigsfile $mergedUnitigInputKUnitigsFile 2> $sequenceCreationErrorFile > $finalSuperReadSequenceFile";
    &runCommandAndExitIfBad ($cmd, $finalSuperReadSequenceFile, 1, "createFastaSuperReadSequences", $finalSuperReadSequenceFile, $goodSuperReadsNamesFile);

    $cmd = "$exeDir/eliminateBadSuperReadsUsingList $joinerOutput $goodSuperReadsNamesFile > $finalReadPlacementFile";
    &runCommandAndExitIfBad ($cmd, $finalReadPlacementFile, 1, "createFinalReadPlacementFile", $finalReadPlacementFile, $fastaSuperReadErrorsFile);
}

$cmd = "touch $successFile";
system ($cmd);

exit (0);

 jumpLibraryCalculations:
$minReadLength = 64;
$maxReadLength = 2047;
$cmd = "$exeDir/createFastaSuperReadSequences -jump-library $workingDirectory $joinerOutput -seqdiffmax $seqDiffMax -min-ovl-len $merLenMinus1 -minreadsinsuperread 1 $mergedUnitigDataFileStr -kunitigsfile $mergedUnitigInputKUnitigsFile -maxunitignumberfile $mergedMaxKUnitigNumberFile -min-read-length $minReadLength -max-read-length $maxReadLength -good-sequence-output-file $finalSuperReadSequenceFile  2> $sequenceCreationErrorFile";
&runCommandAndExitIfBad ($cmd, $finalSuperReadSequenceFile, 1, "createFastaSequencesForJumpingLibrary", $finalSuperReadSequenceFile, $sequenceCreationErrorFile);
    
$cmd = "touch $successFile";
system ($cmd);

sub processArgs
{
    $tableSize = 2000000000;
    $merLen = 31;
    $numProcessors = 16;
    $minReadsInSuperRead = 2;
    $seqDiffMax = 0;
    $help = 0;
    if ($#ARGV < 0) {
	$help = 1; }
    for ($i=0; $i<=$#ARGV; $i++) {
	if ($ARGV[$i] eq "-l") {
	    ++$i;
	    $merLen = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-s") {
	    ++$i;
	    $tableSize = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-t") {
	    ++$i;
	    $numProcessors = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-kunitigsfile") {
	    ++$i;
	    $kUnitigsFile = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-mkudisr") {
	    ++$i;
	    $seqDiffMax = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-minreadsinsuperread") {
	    ++$i;
	    $minReadsInSuperRead = $ARGV[$i];
	    next; }
        elsif ($ARGV[$i] eq "-maxnodes") {
            ++$i;
            $maxNodes = $ARGV[$i];
            next; }
	elsif ($ARGV[$i] eq "-merged-unitig-data-prefix") {
	    ++$i;
	    $mergedUnitigDataPrefix = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-mean-and-stdev-by-prefix-file") {
	    ++$i;
	    $meanAndStdevByPrefixFile = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "--stopAfter") {
	    ++$i;
	    $stopAfter = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-jumplibraryreads") {
	    $jumpLibraryReads = 1;
	    next; }
	elsif ($ARGV[$i] eq "-noreduce") {
            $noReduce = 1;
            next; }
	elsif ($ARGV[$i] eq "-h") {
	    $help = 1;
	    next; }
	elsif ($ARGV[$i] eq "-mikedebug") {
	    $mikedebug = 1;
	    next; }
	elsif (-f $ARGV[$i]) {
	    push (@fastaFiles, $ARGV[$i]);
	    next; }
	else {
	    if (! $workingDirectory) {
		$workingDirectory = $ARGV[$i]; }
	    else {
		print STDERR "Working directory was specified as \"$workingDirectory\", now specified as ",$ARGV[$i],".\nThis message could also occur if an input read file doesn't exist.\n";
		$help = 1; } } }
    if ($#fastaFiles < 0) {
	$help = 1; }
    if (! $kUnitigsFile) {
	print STDERR "A k-unitigs file must be supplied\n";
	$help = 1; }
    if ((! $meanAndStdevByPrefixFile) && (! $jumpLibraryReads)) {
	print STDERR "A file specifying mean and stdev of each library must be given (using the \n  '-mean-and-stdev-by-prefix-file' switch) if -jumplibraryreads is not set.\n";
	$help = 1; }
    if ($help) {
	&giveUsageAndExit; }
}

sub giveUsageAndExit
{
    open (FILE, $0);
    $line = <FILE>;
    while ($line = <FILE>) {
	chomp ($line);
	last unless ($line =~ /^\#/);
	substr ($line, 0, 2) = "";
	print "$line\n"; }
    exit (0);
}

sub returnTableResizeAmount
{
    my ($dataPrefix, $hashFile) = @_;
    
    for ($i=1; 1; $i++) {
	$fn = $dataPrefix . "_$i";
	last unless (-e $fn);
	unlink ($fn); }
    unlink ($hashFile) if ($i > 1);
    return ($i);
}

sub getNumLines
{
    my ($fn) = @_;
    my ($cmd, $line, @flds);
    
    $cmd = "wc -l $fn |";
    open (CRAZYFILE, $cmd);
    $line = <CRAZYFILE>;
    close (CRAZYFILE);
    @flds = split (" ", $line);
    return ($flds[0]);
}

sub reportMinJellyfishTableSizeForKUnitigs
{
    my ($numKMersInKUnitigs, $numKUnitigs, $minSizeNeededForTable);
    
    open (FILE, $totBasesInKUnitigsFile);
    $numKMersInKUnitigs = <FILE>;
    close (FILE);
    chomp ($numKMersInKUnitigs);
    
    open (FILE, $numKUnitigsFile);
    $numKUnitigs = <FILE>;
    chomp ($numKUnitigs);
    close (FILE);
    $numKMersInKUnitigs -= ($numKUnitigs * ($merLenMinus1));
    $minSizeNeededForTable = int ($numKMersInKUnitigs/$maxHashFillFactor + .5);
    return ($minSizeNeededForTable);
}

sub killFiles
{
    my (@filesToKill) = @_;
    my ($file);
    
    for (@filesToKill) {
	$file = $_;
	if (-e $file) {
	    unlink ($file); } }
}

# If localCmd is set, it captures the return code and makes
# sure it ran properly. Otherwise it exits.
# If one sets the fileName then we assume it must exist
# With minSize one can set the minimum output file size
sub runCommandAndExitIfBad
{
    my ($localCmd, $fileName, $minSize, $stepName, @filesCreated) = @_;
    my ($retCode, $exitValue, $sz);
    my ($totSize, $tempFilename, $cmd);
    my ($successFile, $failDir);
    
#    sleep (5); # For testing only
    $failDir = $workingDirectory . "/" . $stepName . ".Failed";
    if (-e $failDir) {
	$cmd = "rm -rf $failDir"; print "$cmd\n"; system ($cmd); }
    $successFile = $workingDirectory . "/" . $stepName . ".success";
    if (-e $successFile) {
	if ($mustRun) {
	    unlink ($successFile); }
	else {
	    print STDERR "Step $stepName already completed. Continuing.\n";
	    goto successfulRun; } }
    else {
	$mustRun = 1; }
    
    if ($localCmd =~ /\S/) {
        print "$localCmd\n";
        system ("time $localCmd");
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
    goto successfulRun unless ($fileName =~ /\S/);
    if ($fileName =~ /\*/) {
	goto multipleFiles; }
    if (! -e $fileName) {
        print STDERR "Output file \"$fileName\" doesn't exist. Bye!\n";
	$retCode = 1;
	goto failingRun; }
    $sz = -s $fileName;
    if ($sz < $minSize) {
        print STDERR "Output file \"$fileName\" is of size $sz, must be at least of size $minSize. Bye!\n";
	$retCode = 1;
	goto failingRun; }
    goto successfulRun;
  multipleFiles:
    $cmd = "ls $fileName |";
    $totSize = 0;
    open (CMD, $cmd);
    while ($tempFilename = <CMD>) {
	chomp ($tempFilename);
	$sz = -s $tempFilename;
	$totSize += $sz; }
    close (CMD);
    if ($totSize < $minSize) {
	print STDERR "The combined output files from \"$fileName\" have a total size of $totSize, must be at least of size $minSize. Bye!\n";
	$retCode = 1;
	goto failingRun; }
  successfulRun:
    if (! -e $successFile) {
	$cmd = "touch $successFile"; print "$cmd\n"; system ($cmd); }
    if ($stopAfter eq $stepName) {
	print STDERR "Stopping after step ${stepName}. Bye!\n";
	exit (0); }
    return;
    
  failingRun:
    $outdir = "$workingDirectory/${stepName}.Failed";
    mkdir ($outdir);
    for (@filesCreated) {
	$cmd = "mv $_ $outdir";
	print "$cmd\n"; system ($cmd); }
    exit ($retCode);
}

sub cleanUpFailureDirectories
{
    my ($tfile, $totFile, $cmd);

    if (! -e $workingDirectory) {
	return; }
    opendir (DIR, $workingDirectory);
    while ($tfile = readdir (DIR)) {
	next unless ($tfile =~ /\.Failed$/);
	$totFile = "$workingDirectory/$tfile";
	next unless (-d $totFile);
	$cmd = "rm -rf $totFile";
	print "$cmd\n"; system ($cmd); }
    closedir (DIR);
}


    
