#!/usr/bin/env perl
# This program is called with arguments (not flags)
# The first argument is the name of the Celera terminator directory (ends in 9-terminator)
# The remaining arguments are the arguments which report all the values of the lengths
# of the ends of the contigs you want to output. These go into output files,
# whose names end in the number of bases kept for each of the contigs.
# If a length is specified multiple times it is only output once.
$CeleraTerminatorDirectory = "/genome2/raid/tri/assembly/rhodobacter/runDir01/CA/9-terminator";
if ($#ARGV >= 0) {
    $CeleraTerminatorDirectory = $ARGV[0]; }
else {
    print STDERR "In ${0}: no Celera terminator directory specified. Bye!\n";
    exit (1); }
if ($#ARGV == 1) {
    @lengthsToOutput = (100); }
else {
    @lengthsToOutput = (@ARGV[1..$#ARGV]); }
$maxLength = findMax(@lengthsToOutput);
$file = "$CeleraTerminatorDirectory/genome.ctg.fasta";
# open (OUTFILE, ">genome.ctg.fwd.pairs");
open (FILE, $file);
$line = <FILE>;
chomp ($line);
($ctg) = ($line =~ /^....(\S+)/);
#print OUTFILE $ctg;
push (@ctgs, $ctg);
while ($line = <FILE>) {
    chomp ($line);
    if ($line =~ /^>/) {
	&doOutputContigSeqs;
	$contigSeq = "";
	($ctg) = ($line =~ /^....(\S+)/);
#	print OUTFILE $ctg;
	push (@ctgs, $ctg);
	next; }
    $contigSeq .= $line;
}
&doOutputContigSeqs;
close (FILE);
# close (OUTFILE);
for (@lengthsToOutput) {
    $length = $_;
    next if ($wasDone{$length});
    # Do the stuff here
    open (OUTFILE, ">genome.ctg.fwd.pairs.$length");
    for ($i=0; $i<=$#ctgs; $i++) {
	$prefix = substr ($outputSeqs[$i], 0, $length);
	$suffix = substr ($outputSeqs[$i], -$length);
	print OUTFILE "$ctgs[$i] $prefix $suffix\n";
    }
    close (OUTFILE);
    # After file finished
    $wasDone{$length} = 1;
}

sub doOutputContigSeqs
{
    if(length($contigSeq) > $maxLength) {
	push (@outputSeqs, substr($contigSeq, 0, $maxLength) . " " . substr($contigSeq, -$maxLength)); }
#	print OUTFILE " ",substr($contigSeq, 0, $maxLength)," ",substr($contigSeq, -$maxLength),"\n"; }
    else {
	push (@outputSeqs, "$contigSeq $contigSeq"); }
#	print OUTFILE" $contigSeq $contigSeq\n"; }
}

sub findMax
{
    my (@lens) = @_;
    my ($maxLength, $i);

    $maxLength = $lens[0];
    for ($i=1; $i<=$#lens; $i++) {
	if ($maxLength < $lens[$i]) {
	    $maxLength = $lens[$i]; }
    }
    return ($maxLength);
}

