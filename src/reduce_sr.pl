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
	print STDERR "Processing ",int(100000/$elapsed)," super reads per second\n";
	$timing=time();
    }

#here we check the k-unitigs and find which super read the current one could be reduced to
#print "Processing sr $line\n";
    %candidate_sr_counts=();
    @f=split('_',$l[0]);
    for($j=0;$j<scalar(@f);$j+=2){
	my $ku=substr($f[$j],0,length($f[$j])-1);
	if(defined($ku_sr_number[$ku])){
	    foreach $v(@{$ku_sr_number[$ku]}){
		$candidate_sr_counts{$v}++;
		}
	}
	else{
	    %candidate_sr_counts=();
	    last;
	}
    }

    if(scalar(keys %candidate_sr_counts)>0){
	@candidate_sr=();
	$desired_count=int((scalar(@f)+1)/2+.001);
	foreach $v(keys %candidate_sr_counts){
#print "Candidate $v $sr[$v] count $candidate_sr_counts{$v}\n";
	    if($candidate_sr_counts{$v}>=$desired_count){
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
	    for($j=0;$j<scalar(@candidate_sr_sorted);$j++){
#print "trying candidate $sr[$candidate_sr_sorted[$j]], number $candidate_sr_sorted[$j]\n";
		last if(index($sr[$candidate_sr_sorted[$j]],$fwd_sr)>-1);
                last if(index($sr[$candidate_sr_sorted[$j]],$rev_sr)>-1);
	    }
	    if($j<scalar(@candidate_sr_sorted)){
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
