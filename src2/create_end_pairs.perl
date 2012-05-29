#!/usr/bin/env perl
#
$CeleraTerminatorDirectory = $ARGV[0];
$file1 = "genome.ctg.fwd.pairs";
$file2 = "$CeleraTerminatorDirectory/genome.posmap.ctgscf";
open(FILE, $file1);
while($line=<FILE>){
    chomp($line);
    @l=split(/\s+/,$line);
    $fwd_b{$l[0]}=$l[1];
    $fwd_e{$l[0]}=$l[2];
}
close(FILE);

$prefix="cc";
$readnumber=0;
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
    
    $readnumber++;
    if($contigOri eq "f"){
	print ">$prefix",$readnumber," $contig f\n",reverse_complement($fwd_b{$contig}),"\n"; }
    else {
	print ">$prefix",$readnumber," $contig r\n",$fwd_e{$contig},"\n"; }
    $readnumber++;
}

sub reverse_complement
{
    my $string=$_[0];
    my $rev_comp_sequence="";
    for(my $i=length($string);$i>=0;$i--)
    {
        if($if_qual==1)
        {
	    $rev_comp_sequence=$rev_comp_sequence.substr($string,$i,1);
        }
        else
        {
	    my $t=substr($string,$i,1);
	    $t=~tr/ACGTNacgtn/TGCANtgcan/;
	    $rev_comp_sequence=$rev_comp_sequence.$t;
        }
    }
    return($rev_comp_sequence);
}

