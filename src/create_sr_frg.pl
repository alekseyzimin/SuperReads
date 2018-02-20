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

$rn="";
$max_len=$ARGV[0];
$shooting_index=0;
if($ARGV[1] eq ""){
$suffix="super-read";
}else{
$suffix=$ARGV[1];
}
while($line=<STDIN>){
    if($line =~ /^>/){
	if(not($rn eq "")){
	    $l=length($seq);
	    $rev_seq=reverse_complement($seq);
	    $seq=$rev_seq lt $seq ? $rev_seq : $seq;
	    if($l<$max_len){
		print "$rn\n$seq\n";
	    }else{
                $max_len_local=int($l/int($l/$max_len+1));
		$offset=int(($max_len_local-1)/2);
                $offset=$max_len_local-10000 if($max_len_local-$offset>10000);
	    	for($i=0;$i<$l;$i+=$offset){
			print "$rn.$max_len_local.$i\n",substr($seq,$i,$max_len_local),"\n";
		}
	    }
	}
	chomp($line); 
	@l=split(/\s+/,$line);
	if(length($l[0])>100){
            print STDERR "SR$shooting_index ",substr($l[0],1),"\n";
            $rn=">SR".($shooting_index).":".$suffix;
            $shooting_index++;
        }else{
        $rn=$l[0].":".$suffix;
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
            if($l<$max_len){
                print "$rn\n$seq\n";
            }else{
                $max_len_local=int($l/int($l/$max_len+1));
                $offset=int(($max_len_local-1)/2);
                $offset=$max_len_local-10000 if($max_len_local-$offset>10000);
                    for($i=0;$i<$l;$i+=$offset){
                        print "$rn.$max_len_local.$i\n",substr($seq,$i,$max_len_local),"\n";
                    }   
            }

sub reverse_complement {
my $seq=$_[0];
$seq=~tr/ACGTNacgtn/TGCANtgcan/;
$seq=reverse($seq);
return($seq);
}
