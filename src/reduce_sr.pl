#!/usr/bin/env perl
#
while($line=<STDIN>){
chomp($line);
@l=split(/\s+/,$line);
push (@sr,$l[0]);
$len{$l[0]}=$l[1];
}


for($j=0;$j<scalar(@sr);$j++){
push(@sr_rev,reverse_sr($sr[$j]));
$valid_indices{$j}=1;
}

for($i=1;$i<50;$i+=2){
#print "Size: $i\n";

for($j=0;$j<scalar(@sr);$j++){

next if(not(defined($valid_indices{$j})));
@s=split(/_/,$sr[$j]);
next if(scalar(@s)!=$i);
delete $valid_indices{$j};

#printf "Trying to merge $sr[$j]\n";
#print "fwd: $sr_tmp rev: $sr_tmp_rev\n";
foreach $k(keys %valid_indices){
#print "candidate $sr[$k]\n";
if($sr[$k] =~ /^($sr[$j])/||$sr[$k] =~ ("_".$sr[$j]) || $sr[$k] =~ /^($sr_rev[$j])/||$sr[$k] =~ ("_".$sr_rev[$j])){
#print "found match $sr[$j] $sr[$k]\n";
if((not defined($replace{$sr[$j]}))||$len{$replace{$sr[$j]}}<$len{$sr[$k]}){
$replace{$sr[$j]}=$sr[$k];
}
}
}
if(defined($replace{$sr[$j]})){
#print "merged $sr[$j] into $replace{$sr[$j]}\n";
push(@merged,$sr[$j]);
}
}
}

my $temp_ref;

for($j=0;$j<scalar(@merged);$j++){
$temp_ref=$merged[$j];
while(1){
if(defined($replace{$temp_ref})){
$temp_ref=$replace{$temp_ref};
}
else
{
last;
}
}
print "$merged[$j] $temp_ref\n";
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
