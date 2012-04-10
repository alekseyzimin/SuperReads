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
#assume that super reads are sorted by size, longest to shortest
$i=0;
$ll=0;
@sr=();
$timing=time();
@ku_sr_number=();
while($line=<STDIN>){
    chomp($line);
    @l=split(/\s+/,$line);
    $ll++;
    if($ll%100000==0){
	print STDERR "Processed $ll super reads, irreducible $i\n";
	$elapsed=time()-$timing;
	print STDERR "Processing ",int(500000/$elapsed)," super reads per second\n";
	$timing=time();
    }

#here we check the first and last k-unitig and find which super read the current one could be reduced to
#print "Processing sr $line\n";
    %candidate_sr_counts=();
    @f=split('_',$l[0]);
    $first_k_u=0;
    $ku=substr($f[$#f],0,length($f[$#f])-1);
    if(defined($ku_sr_number[$ku])){
           foreach $v(@{$ku_sr_number[$ku]}){
                $candidate_sr_counts{$v}=1;
                }
	}
        else{
            %candidate_sr_counts=();
            $first_k_u=scalar(@f);
        }

     if($first_k_u==0){
	$ku=substr($f[0],0,length($f[0])-1);
	if(defined($ku_sr_number[$ku])){
	    foreach $v(@{$ku_sr_number[$ku]}){
		if(defined($candidate_sr_counts{$v})){
		$candidate_sr_counts{$v}++;
		}
		}
	}
	else{
	    %candidate_sr_counts=();
	}
    }

    if(scalar(keys %candidate_sr_counts)>0){
	@candidate_sr=();
	foreach $v(keys %candidate_sr_counts){
#print "Candidate $v $sr[$v] count $candidate_sr_counts{$v}\n";
	    if($candidate_sr_counts{$v}>=2){
		push(@candidate_sr,$v);
#print "Final candidate $v $sr[$v]\n";
	    }
	}


	if(scalar(@candidate_sr)>0){
#print "Found ",scalar(keys %candidate_sr)," candidates for merging\n";
	    @candidate_sr_sorted = sort {$a>$b} (@candidate_sr);

#look for the match
	    $fwd_sr=$l[0];
            $rev_sr=reverse_sr(@f);
	    $flag=0;
 	    for($j=0;$j<scalar(@candidate_sr_sorted) && $j<20 ;$j++){
#print "trying candidate $sr[$candidate_sr_sorted[$j]], number $candidate_sr_sorted[$j]\n";
		if(index($sr[$candidate_sr_sorted[$j]],$fwd_sr)>-1){
		$flag=1;
		last;
		}
                if(index($sr[$candidate_sr_sorted[$j]],$rev_sr)>-1){
                $flag=1;
                last;
		}
	    }
	    if($flag==1){
		print "$fwd_sr $sr[$candidate_sr_sorted[$j]]\n";
		next;
	    }
	}
    }
#print "Maximal sr $line, i=$i\n";
    push(@sr,$l[0]);
    for($j=0;$j<scalar(@f);$j+=2){ 
	my $ku=substr($f[$j],0,length($f[$j])-1);
	push(@{$ku_sr_number[$ku]},$i);
    }
    $i++;
}

sub reverse_sr
{
    my @fsr=@_;
    my $new_sr;
    my $flag=0;
    my $dir;
    for(my $i=scalar(@fsr)-1;$i>-1;$i--){
    if($flag==0){
    if(substr($fsr[$i],length($fsr[$i])-1) eq "R"){
	$fsr[$i]=~s/R/F/;
    }
    else{
	$fsr[$i]=~s/F/R/;
    }
    $new_sr.=$fsr[$i];
}
else{
    $new_sr.="_".$fsr[$i]."_";
}
$flag=1-$flag;
}
return $new_sr;
}
