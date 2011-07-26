#!/usr/bin/env perl
my $flag=0;

while($line=<STDIN>)
{
chomp($line);
@F=split(/\s+/,$line);
if($flag==0)
{
$sr_pair=$F[1];
if(int(substr($F[0],2))%2==0)
{
$insert=$F[0];
}
else
{
$insert=substr($F[0],0,2).int(substr($F[0],2))-1;
}
$flag=1;
}
elsif($flag==1)
{
$srh{"$sr_pair $F[1]"}.="$insert ";
$srh{"$F[1] $sr_pair"}.="$insert ";
push(@pairs,"$sr_pair $F[1]");
$flag=0;
}
}

foreach $v (@pairs)
{
@f=split(/ /,$v);
if(defined($srh{"$f[0] $f[1]"}) && defined($srh{"$f[1] $f[0]"}))
{
print "$v : ",$srh{"$f[0] $f[1]"},"\n";
}
}


