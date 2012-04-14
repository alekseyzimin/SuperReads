#!/usr/bin/env perl
# SuperRead pipeline
# Copyright (C) 2012  Genome group at University of Maryland.
# 
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

