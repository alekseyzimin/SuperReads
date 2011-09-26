#!/usr/bin/env perl
#
open(FILE,"genome.ctg.fwd.pairs");
while($line=<FILE>){
chomp($line);
@l=split(/\s+/,$line);
$fwd_b{$l[0]}=$l[1];
$fwd_e{$l[0]}=$l[2];
}
close(FILE);

$prefix="cc";
$readnumber=0;
while($line=<STDIN>){
chomp($line);
@l=split(/\s+/,$line);
if($l[1] eq "f"){
print ">$prefix",$readnumber,"\n",$fwd_e{$l[0]},"\n";
}
else{
print ">$prefix",$readnumber,"\n",reverse_complement($fwd_b{$l[0]}),"\n";
}

$readnumber++;
if($l[3] eq "f"){
print ">$prefix",$readnumber,"\n",reverse_complement($fwd_b{$l[2]}),"\n";
}
else{
print ">$prefix",$readnumber,"\n",$fwd_e{$l[2]},"\n";
}
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

