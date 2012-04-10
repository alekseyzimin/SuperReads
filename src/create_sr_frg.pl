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



#!/usr/bin/env perl
$rn="";
$shooting_index=0;
while($line=<STDIN>){
    if($line =~ /^>/){
	if(not($rn eq "")){
	    $l=length($seq);
	    $rev_seq=reverse_complement($seq);
	    $seq=$rev_seq lt $seq ? $rev_seq : $seq;
	    if($l<2048){
		print "$rn\n$seq\n";
	    }else{
		my @f=split(//,$seq);
		$k=0;
		$offset=1550;
		while(1){
		    print "$rn.",$k*$offset,"\n";
		    for($i=$k*$offset;($i<$k*$offset+2047&&$i<=$#f);$i++){
			print $f[$i];
		    }
		    $k++;
		    print "\n";
		    last if($i>$#f);
		}
	    }
	}
	chomp($line); 
	@l=split(/\s+/,$line);
	if(length($l[0])>100){
            print STDERR "SHOOTING_LONG_NAME_$shooting_index ",substr($l[0],1),"\n";
            $rn=">SHOOTING_LONG_NAME_".($shooting_index).":super-read";
            $shooting_index++;
        }else{
        $rn=$l[0].":super-read";
	}
	$seq="";
    }else{
	chomp($line); 
	$seq.=$line;
    }
}
#do not forget the last one!!!
$l=length($seq);
$rev_seq=reverse_complement($seq);
$seq=$rev_seq lt $seq ? $rev_seq : $seq;
if($l<2048){
    print "$rn\n$seq\n";
}else{
    my @f=split(//,$seq);
    $k=0;
    $offset=1550;
    while(1){
	print "$rn.",$k*$offset,"\n";
	for($i=$k*$offset;($i<$k*$offset+2047&&$i<=$#f);$i++){
	    print $f[$i];
	}
	$k++;
	print "\n";
	last if($i>$#f);
    }
}

sub reverse_complement{
    my $string=$_[0];
    my $rev_comp_sequence="";
    for(my $i=length($string);$i>=0;$i--){
        if($if_qual==1)        {
	    $rev_comp_sequence=$rev_comp_sequence.substr($string,$i,1);
        }else        {
	    my $t=substr($string,$i,1);
	    $t=~tr/ACGTNacgtn/TGCANtgcan/;
	    $rev_comp_sequence=$rev_comp_sequence.$t;
        }
    }
    return($rev_comp_sequence);
}


