#!/usr/bin/env perl
# 2 args
# One arg specifies where the Celera terminator directory is
# The other arg is the length of the sequence at the kept at the ends of the contigs
for (@ARGV) {
    $arg = $_;
    if ($arg =~ /^\d+$/) {
	$seqLengthValue = $arg; }
    else {
	$CeleraTerminatorDirectory = $arg; }
}

$file1 = "genome.ctg.fwd.pairs";
$file2 = "$CeleraTerminatorDirectory/genome.posmap.ctgscf";
if ($seqLengthValue) {
    $file1 .= ".$seqLengthValue"; }
open(FILE, $file1);
while($line=<FILE>){
    chomp($line);
    @flds = split(" ",$line);
    $contig = $flds[0];
    $fwd_b{$contig} = $flds[1];
    $fwd_e{$contig} = $flds[2];
}
close(FILE);

$prefix = "cc";
$readnumber = 0;
open (FILE, $file2);

while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $contig = $flds[0];
    $scaffold = $flds[1];
    $contigOri = $flds[4];
    if ($scaffold eq $previousScaffold) {
	outputSequencePair ($previousContig, $previousContigOri, $contig, $contigOri); }
    else {
	$previousScaffold = $scaffold; }
    $previousContig = $contig;
    $previousContigOri = $contigOri;
}

sub outputSequencePair
{
    my ($previousContig, $previousContigOri, $contig, $contigOri) = @_;
    
    if($previousContigOri eq "f") {
	print ">$prefix",$readnumber," $previousContig f\n",$fwd_e{$previousContig},"\n"; }
    else {
	print ">$prefix",$readnumber," $previousContig r\n",reverse_complement($fwd_b{$previousContig}),"\n";
    }
    ++$readnumber;

    if ($contigOri eq "f") {
	print ">$prefix",$readnumber," $contig f\n",reverse_complement($fwd_b{$contig}),"\n"; }
    else {
	print ">$prefix",$readnumber," $contig r\n",$fwd_e{$contig},"\n"; }
    ++$readnumber;
}

sub reverse_complement
{
    my ($string) = @_;

    $string =~ tr/ACGTacgt/TGCAtgca/;
    $string = reverse ($string);
    return ($string);
}

