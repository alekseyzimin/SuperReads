#!/usr/bin/env perl
# (1)
$CeleraTerminatorDirectory = "/genome2/raid/tri/assembly/rhodobacter/runDir01/CA/9-terminator";
if ($#ARGV >= 0) {
    $CeleraTerminatorDirectory = $ARGV[0]; }
else {
    print STDERR "In ${0}: no Celera terminator directory specified. Bye!\n";
    exit (1); }
$file = "$CeleraTerminatorDirectory/genome.ctg.fasta";
open (OUTFILE, ">genome.ctg.fwd.pairs");
open (FILE, $file);
$line = <FILE>;
chomp ($line);
($ctg) = ($line =~ /^....(\S+)/);
print OUTFILE $ctg;
while ($line = <FILE>) {
    chomp ($line);
    if ($line =~ /^>/) {
	&doOutputContigSeqs;
	$contigSeq = "";
	($ctg) = ($line =~ /^....(\S+)/);
	print OUTFILE $ctg;
	next; }
    $contigSeq .= $line;
}
&doOutputContigSeqs;
close (FILE);
close (OUTFILE);

sub doOutputContigSeqs
{
    if(length($contigSeq) > 100) {
	print OUTFILE " ",substr($contigSeq, 0, 100)," ",substr($contigSeq, -100),"\n"; }
    else {
	print OUTFILE" $contigSeq $contigSeq\n"; }
}
