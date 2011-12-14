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
# -M minCoverageToContinueAUnitig : same as the parameter into create_k_unitigs (3)
# -m minCoverageToCreateAFork : same as the parameter into create_k_unitigs (2)
# -kunitigsfile filename : a user-given k-unitigs file; otherwise we calculate
# -mean-and-stdev-by-prefix-file filename : a file giving mate info about each
#                      library. Each line is the 2-letter prefix for the reads
#                      in the library followed by its mean and stdev. This
#                      file is mandatory unless -jumplibraryreads is specified
# -mkudisr numBaseDiffs : max base diffs between overlapping k-unitigs in super-reads (0)
# -minreadsinsuperread minReads : super-reads containing fewer than numReads
#                                reads will be eliminated (2)
# -noclean : don't clean up the files afterwards
# -mikedebug : don't kill off intermediate results
# -jumplibraryreads : we are generating for jump-library reads; a k-unitigs
#                                 file must be specified
# -elim-dupls : eliminate duplicate reads ahead of time (not implemented yet).
# -h : help 
use File::Basename;
use Cwd;
$exeDir = dirname ($0);
# $exeDir = "/home/tri/superReadPipeline/SuperReads/build/src2/";
$pwd = cwd;
if ($exeDir !~ /^\//) {
    $exeDir = "$pwd/$exeDir"; }

&processArgs;
$merLenMinus1 = $merLen - 1;

$maxHashFillFactor = .8;

$successFile = "$workingDirectory/superReads.success";
unlink ($successFile) if (-e $successFile);

if (! -d $workingDirectory) {
    $cmd = "mkdir $workingDirectory";
    print "$cmd\n"; system ($cmd); }
$jellyfishKUnitigDataPrefix = "$workingDirectory/organismMerCountsForKUnitigs";
$jellyfishKUnitigHashFile = $jellyfishKUnitigDataPrefix . "_0";
$kUnitigLengthsFile = "$workingDirectory/kUnitigLengths.txt";
# The following stores the actual number of k-unitigs
$numKUnitigsFile = "$workingDirectory/numKUnitigs.txt";
# The following stores the largest k-unitig number (+1)
$maxKUnitigNumberFile = "$workingDirectory/maxKUnitigNumber.txt";
$totBasesInKUnitigsFile = "$workingDirectory/totBasesInKUnitigs.txt";
$totReadFile = "$workingDirectory/inputReads.fasta";
if ($#fastaFiles == 0) {
    $totReadFile = $fastaFiles[0]; }
$numReadsFile = "$workingDirectory/numReads.txt";
$prefixForOverlapsBetweenKUnitigs = "$workingDirectory/overlap";
$kUnitigOverlapsFile = "${prefixForOverlapsBetweenKUnitigs}.overlaps";

$myProgOutput1prefix = "$workingDirectory/newTestOutput";
$joinerOutputPrefix = "$workingDirectory/readPositionsInSuperReads";
$myProgOutput1_1 = "$workingDirectory/testOutput.nucmerLinesOnly.txt";
$myProgOutput1_1prefix = "$workingDirectory/newTestOutput.nucmerLinesOnly";
$sequenceCreationErrorFile = "$workingDirectory/createFastaSuperReadSequences.errors.txt";
$myProgOutput2 = "$workingDirectory/readPlacementsInSuperReads.postMateMerge.read.superRead.offset.ori.txt";
$myProgOutput3 = "$workingDirectory/superReadCounts.count.superRead.txt";
$finalSuperReadSequenceFile = "$workingDirectory/superReadSequences.fasta";
$finalReadPlacementFileUsingReadNumbers = "$workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.usingReadNumbers.txt";
$finalReadPlacementFile = "$workingDirectory/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt";
$reducedReadPlacementFile = "$workingDirectory/readPlacementsInSuperReads.reduced.read.superRead.offset.ori.txt";

if ($jumpLibraryReads) {
    $myProgOutput1_1 = "$workingDirectory/testOutput.nucmerLinesOnly.jumpLibrary.txt";
    $myProgOutput1_1prefix = "$workingDirectory/newTestOutput.nucmerLinesOnly.jumpLibrary";
    $myProgOutput8 = "$workingDirectory/superReadSequences.jumpLibrary.bothSidesExtended.fasta";
    $myProgOutput12 = "$workingDirectory/readPlacementsInSuperReads.forJumpLibraryWBothSidesExtended.read.superRead.offset.ori.txt";
    $finalSuperReadSequenceFile = "$workingDirectory/superReadSequences.jumpLibrary.fasta";
}

# goto startHere;

# We now require that a k-unitigs file was passed on the command line
if ($kUnitigsFile !~ /^\//) {
    $kUnitigsFile = "$pwd/$kUnitigsFile"; }

# In addition to obvious output file, this also generates the files
# numKUnitigs.txt, maxKUnitigNumber.txt, and totBasesInKUnitigs.txt in
# $workingDirectory
$cmd = "cat $kUnitigsFile | $exeDir/getLengthStatisticsForKUnitigsFile.perl $workingDirectory > $kUnitigLengthsFile";
&runCommandAndExitIfBad ($cmd, $kUnitigLengthsFile, 1);

$minSizeNeededForTable = &reportMinJellyfishTableSizeForKUnitigs;
redoKUnitigsJellyfish:
    $cmd = "jellyfish count -m $merLen -r -o $jellyfishKUnitigDataPrefix -c 6 -p 126 --both-strands -s $minSizeNeededForTable -t $numProcessors $kUnitigsFile";
&runCommandAndExitIfBad ($cmd, $jellyfishKUnitigHashFile, 1);

$tableResizeFactor = &returnTableResizeAmount ($jellyfishKUnitigDataPrefix, $jellyfishKUnitigHashFile);
if ($tableResizeFactor > 1) {
    $tableSize *= 2;
    print "Resizing the table to $tableSize for the k-unitig jellyfish run\n";
    goto redoKUnitigsJellyfish; }

$cmd = "$exeDir/findMatchesBetweenKUnitigsAndReads $jellyfishKUnitigHashFile -t $numProcessors -p $myProgOutput1_1prefix $kUnitigsFile $maxKUnitigNumberFile $totReadFile";
&runCommandAndExitIfBad ($cmd, "", 0);
if (! $mikedebug) { &killFiles ($jellyfishKUnitigHashFile, $totReadFile); }

if ($jumpLibraryReads) {
    goto jumpLibraryCalculations; }

open (FILE, $maxKUnitigNumberFile); $maxKUnitigNumber = <FILE>; chomp ($maxKUnitigNumber); close (FILE);
$cmd = "$exeDir/createKUnitigMaxOverlaps $kUnitigsFile -kmervalue $merLen -largest-kunitig-number ".(int($maxKUnitigNumber)+1)." $prefixForOverlapsBetweenKUnitigs";
&runCommandAndExitIfBad($cmd, $kUnitigOverlapsFile, 1);

# Do the shooting method here
$cmd = "$exeDir/joinKUnitigs_v3 -mean-and-stdev-by-prefix-file $meanAndStdevByPrefixFile -unitig-lengths-file $kUnitigLengthsFile -num-kunitigs-file $maxKUnitigNumberFile -overlaps-file $kUnitigOverlapsFile -min-overlap-length $merLenMinus1 -report-paths -prefix $joinerOutputPrefix -num-file-names $numProcessors $myProgOutput1_1prefix";
print "$cmd\n";
system("time $cmd");

#AZ -- put super read reduction before the counts
my $numConcurrentJobs=2;
open(FILE,">$workingDirectory/commands.sh");
print FILE "#!/bin/bash\n";
for($k=0;$k<$numProcessors;$k+=$numConcurrentJobs){
    for($j=0;$j<$numConcurrentJobs;$j++){
	print FILE "if [ -e ${joinerOutputPrefix}_".($k+$j)." ];then cat ${joinerOutputPrefix}_".($k+$j)."  | $exeDir/getSuperReadInsertCountsFromReadPlacementFile >  $workingDirectory/superReadCounts_".($k+$j).";fi & \n";
	print FILE "PID$j=\$!\n"
    }
    print FILE "wait ";
    for($j=0;$j<$numConcurrentJobs;$j++){
	print FILE "\$PID$j ";
    }
    print FILE "\n";
}
close(FILE);
system("cat $workingDirectory/commands.sh");
$cmd = "chmod 0755 $workingDirectory/commands.sh;  $workingDirectory/commands.sh";
print "$cmd\n";
system($cmd);

open(FILE,">$workingDirectory/commands.sh");
print FILE "#!/bin/bash\n";
$numConcurrentJobs=16;
print FILE "sort -m -k2,2 ";
for($k=0;$k<$numProcessors;$k+=$numConcurrentJobs){
    print FILE "<(sort -m -k2,2 $workingDirectory/superReadCounts_{".$k;
    for($j=1;$j<$numConcurrentJobs;$j++){
	print FILE ",".($k+$j) if( -e "$workingDirectory/superReadCounts_".($k+$j));
    }
    print FILE "}| awk 'BEGIN{l=\"-1\";c=0}{if(l==\$2){c+=\$1}else{if(l!=\"-1\"){print c\" \"l;}c=\$1;l=\$2}}END{print c\" \"l}') ";
}
print FILE " | awk 'BEGIN{l=\"-1\";c=0}{if(l==\$2){c+=\$1}else{if(l!=\"-1\" && c >= $minReadsInSuperRead ){print c\" \"l;}c=\$1;l=\$2}}END{print c\" \"l}' >  $workingDirectory/superReadCounts.all\n";
close(FILE);

system("cat $workingDirectory/commands.sh");
$cmd = "chmod 0755 $workingDirectory/commands.sh;  $workingDirectory/commands.sh";
print "$cmd\n";
system($cmd);

$cmd = "cat $workingDirectory/superReadCounts.all | $exeDir/createFastaSuperReadSequences $workingDirectory /dev/fd/0 -seqdiffmax $seqDiffMax -min-ovl-len $merLenMinus1 -minreadsinsuperread $minReadsInSuperRead -good-sr-filename $workingDirectory/superReadNames.txt -kunitigsfile $kUnitigsFile 2> $sequenceCreationErrorFile | tee $finalSuperReadSequenceFile.all | perl -ane 'BEGIN{my \$seq_length=0}{if(\$F[0] =~ /^>/){if(\$seq_length>0){print \$seq_length,\"\\n\";} print substr(\$F[0],1),\" \";\$seq_length=0;}else{\$seq_length+=length(\$F[0]);}}END{if(\$seq_length>0){print \$seq_length,\"\\n\";}}' | sort -nrk2,2 -S 40% -T ./ > $workingDirectory/sr_sizes.tmp";
&runCommandAndExitIfBad ($cmd,"$workingDirectory/sr_sizes.tmp", 1);

$cmd = "cat $workingDirectory/sr_sizes.tmp| $exeDir/reduce_sr $maxKUnitigNumber  > $workingDirectory/reduce.tmp";
&runCommandAndExitIfBad ($cmd,"$workingDirectory/reduce.tmp", 1);

$cmd = "cat ${joinerOutputPrefix}_* | $exeDir/eliminateBadSuperReadsUsingList /dev/fd/0 $workingDirectory/superReadNames.txt | perl -e '{open(FILE,\$ARGV[0]);while(\$line=<FILE>){chomp(\$line);\@F=split(\" \",\$line);\$sr{\$F[0]}=\$F[1]} while(\$line=<STDIN>){chomp(\$line);\@l=split(\" \",\$line);if(defined(\$sr{\$l[1]})){print \"\$l[0] \",\$sr{\$l[1]},\" \$l[2] \$l[3]\\n\";}else{print \"\$line\\n\";}}}' $workingDirectory/reduce.tmp >  $finalReadPlacementFile"; 
&runCommandAndExitIfBad ($cmd, $finalReadPlacementFile, 1);

$cmd = "awk '{print \$1}' $workingDirectory/reduce.tmp > $workingDirectory/sr_to_eliminate.tmp; $exeDir/extractreads_not.pl $workingDirectory/sr_to_eliminate.tmp $finalSuperReadSequenceFile.all 1 > $finalSuperReadSequenceFile";
&runCommandAndExitIfBad ($cmd,$finalSuperReadSequenceFile, 1);

$cmd = "touch $successFile";
system ($cmd);

exit (0);

jumpLibraryCalculations:
# Must run createFastaSuperReadSequences here

    $cmd = "$exeDir/outputSuperReadSeqForJumpLibrary.perl $myProgOutput8 $myProgOutput12 > $finalSuperReadSequenceFile";
&runCommandAndExitIfBad ($cmd, $finalSuperReadSequenceFile, 1);

$cmd = "touch $successFile";
system ($cmd);

sub processArgs
{
    $tableSize = 2000000000;
    $merLen = 31;
    $numProcessors = 16;
    $minCoverageToCreateAFork = 2;
    $minCoverageToContinueAUnitig = 3;
    $minReadsInSuperRead = 2;
    $seqDiffMax = 0;
    $elimDupls = 0;
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
	elsif ($ARGV[$i] eq "-M") {
	    ++$i;
	    $minCoverageToContinueAUnitig = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-kunitigsfile") {
	    ++$i;
	    $kUnitigsFile = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-m") {
	    ++$i;
	    $minCoverageToCreateAFork = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-mkudisr") {
	    ++$i;
	    $seqDiffMax = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-minreadsinsuperread") {
	    ++$i;
	    $minReadsInSuperRead = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-mean-and-stdev-by-prefix-file") {
	    ++$i;
	    $meanAndStdevByPrefixFile = $ARGV[$i];
	    next; }
	elsif ($ARGV[$i] eq "-elim-dupls") {
	    $elimDupls = 1;
	    next; }
	elsif ($ARGV[$i] eq "-jumplibraryreads") {
	    $jumpLibraryReads = 1;
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
    my ($localCmd, $fileName, $minSize) = @_;
    my ($retCode, $exitValue, $sz);

    if ($localCmd =~ /\S/) {
        print "$localCmd\n";
        system ("time $localCmd");
        $retCode = $?;
        if ($retCode == -1) {
            print STDERR "failed to execute: $!\n";
            exit ($retCode); }
        elsif ($retCode & 127) {
            printf STDERR "child died with signal %d, %s coredump\n",
            ($retCode & 127), ($retCode & 128) ? 'with' : 'without';
            exit ($retCode); }
        else {
            $exitValue = $retCode >> 8;
            if ($exitValue == 255) { $exitValue = -1; }
            if ($exitValue != 0) {
                printf STDERR "child exited with value %d\n", $exitValue;
                print STDERR "Command \"$localCmd\" failed. Bye!\n";
                exit ($exitValue); }
        }
    }
    return unless ($fileName =~ /\S/);
    if (! -e $fileName) {
        print STDERR "Output file \"$fileName\" doesn't exist. Bye!\n";
        exit (1); }
    $sz = -s $fileName;
    if ($sz < $minSize) {
        print STDERR "Output file \"$fileName\" is of size $sz, must be at least of size $minSize. Bye!\n";
        exit (1); }
}

