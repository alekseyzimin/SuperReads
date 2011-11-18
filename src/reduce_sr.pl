#!/usr/bin/env perl
#assume that super reads are sorted by size, longest to shortest
$i=0;
@sr=();
@ku_sr_number=();
while($line=<STDIN>){
chomp($line);
@l=split(/\s+/,$line);
#here we check the k-unitigs and find which super read the current one could be reduced to
#print "Processing sr $line\n";
%candidate_sr_counts=();
@f=split('_',$l[0]);
for($j=0;$j<scalar(@f);$j+=2){
my $ku=substr($f[$j],0,length($f[$j])-1);
if(defined($ku_sr_number[$ku])){
@c=split(' ',$ku_sr_number[$ku]);
#print "Preliminary candidates $ku_sr_number[$ku]\n";
foreach $v(@c){
#print "Preliminary candidate $sr[$v] $v\n";
$candidate_sr_counts{$v}++;
}
}
}
%candidate_sr=();
foreach $v(keys %candidate_sr_counts){
if($candidate_sr_counts{$v}>=int((scalar(@f)+1)/2+0.1)){
#print "Final candidate $v\n";
$candidate_sr{$v}=1;
}
}

if(scalar(keys %candidate_sr)>0){
#print "Found ",scalar(keys %candidate_sr)," candidates for merging\n";
@candidate_sr_sorted = sort {$a>$b} (keys %candidate_sr);
#look for the match
$fwd_sr=$l[0];
$rev_sr=reverse_sr($l[0]);
for($j=0;$j<scalar(@candidate_sr_sorted);$j++){
#print "trying candidate $sr[$candidate_sr_sorted[$j]], number $candidate_sr_sorted[$j]\n";
last if($sr[$candidate_sr_sorted[$j]] =~ /^($fwd_sr)/||$sr[$candidate_sr_sorted[$j]] =~ ("_".$fwd_sr) || $sr[$candidate_sr_sorted[$j]] =~ /^($rev_sr)/||$sr[$candidate_sr_sorted[$j]] =~ ("_".$rev_sr));
}
if($j<scalar(@candidate_sr_sorted)){
print "$fwd_sr $sr[$candidate_sr_sorted[$j]]\n";
next;
}
}
#print "Maximal sr $line, i=$i\n";
push(@sr,$l[0]);
for($j=0;$j<scalar(@f);$j+=2){ 
my $ku=substr($f[$j],0,length($f[$j])-1);
$ku_sr_number[$ku].="$i ";
}
$i++;
}

sub reverse_sr
{
my @fsr=split(/_/,$_[0]);
my $new_sr="";
my $flag=0;
my $dir;
for(my $i=scalar(@fsr)-1;$i>-1;$i--)
{
if($flag==0)
{
$dir="R";
$dir="F" if(substr($fsr[$i],length($fsr[$i])-1) eq "R");
$new_sr.=substr($fsr[$i],0,length($fsr[$i])-1).$dir;
}
else
{
$new_sr.="_".$fsr[$i]."_";
}
$flag=1-$flag;
}
return $new_sr;
}
