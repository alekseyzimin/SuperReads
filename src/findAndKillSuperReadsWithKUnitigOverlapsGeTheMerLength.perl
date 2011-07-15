#!/usr/bin/perl
# This takes the input from the first file named on the 
# command line, (a file like
# superReadGroupsForEachReadWhichHasAGroup.txt) and kills those that
# have an overlap between k-unitigs of the superRead that is at least the
# length of the k-mer (entered as the second arg). The list of super-reads
# that pass the tests are output to the third arg, and the input about the
# super-reads that were eliminated (including how many reads "belonged"
# to this super-read) are the fourth arg
$inputSuperReadGroupFile = $ARGV[0];
$kmerLen = $ARGV[1];
$outputSuperReadGroupFile = $ARGV[2];
$killedSuperReadFile = $ARGV[3];

open (FILE, $inputSuperReadGroupFile);
open (OUTFILE, ">$outputSuperReadGroupFile");
while ($lineHold = <FILE>) {
    $line = <FILE>;
    chomp ($line);
    @flds = split (" ", $line);
    $bad = 0;
    for ($i=2; $i<=$#flds; $i+=3) {
	if ($flds[$i] >= $kmerLen) {
	    $bad = 1; } }
    if ($bad) {
	if (! $badCount{$line}) {
	    push (@badSuperReads, $line); }
	++$badCount{$line}; }
    else {
	print OUTFILE $lineHold, "$line\n"; }
}

@indices = (0..$#badSuperReads);
@indices = sort spcl @indices;

open (OUTFILE2, ">$killedSuperReadFile");
for (@indices) {
    $index = $_;
    $badSuperRead = $badSuperReads[$index];
    @flds = split (" ", $badSuperRead);
    $badSuperReadName = $flds[0] . $flds[1];
    for ($i=2; $i<$#flds; $i+=3) {
	$badSuperReadName .= ("_" . $flds[$i] . "_" . $flds[$i+1] . $flds[$i+2]); }
    print OUTFILE2 $badCount{$badSuperRead}," ",$badSuperReadName,"\n";
}

sub spcl
{
    if ($badCount{$badSuperReads[$a]} != $badCount{$badSuperReads[$b]}) {
	return ($badCount{$badSuperReads[$b]} <=> $badCount{$badSuperReads[$a]}); }
    return ($a <=> $b);
}

