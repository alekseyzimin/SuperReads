#!/usr/bin/env perl
$rn="";
while($line=<STDIN>)
{
if($line =~ /^>/)
{
if(not($rn eq ""))
{
$l=length($seq);
if($l<2048)
{
print "$rn\n$seq\n";
}
else
{
my @f=split(//,$seq);
$k=0;
$offset=int(1000+rand(100));
while(1)
{
print "$rn.$k\n";
for($i=$k*$offset;($i<$k*$offset+2047&&$i<=$#f);$i++)
{
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
$rn=$l[0].":super-read";
$seq="";
}
else
{
chomp($line); 
$seq.=$line;
}
}
