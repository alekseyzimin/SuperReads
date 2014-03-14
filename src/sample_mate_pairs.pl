#!/usr/bin/env perl
#
#
my $numer=$ARGV[0];
my $denom=$ARGV[1];
my $ssrand=$ARGV[2];
$ssrand=1 if($ssrand=="");
srand($ssrand);
die("Bad arguments $numer $denom") if($denom<=0||$numer<=0);

my $ratio=$numer/$denom;
if($ratio<1){
while($line=<STDIN>){
if(rand(1)<$ratio){
print $line;
$line=<STDIN>;
print $line;
$line=<STDIN>;
print $line;
$line=<STDIN>;
print $line;
}else{
$line=<STDIN>;
$line=<STDIN>;
$line=<STDIN>;
}
}
}else{
while($line=<STDIN>){
print $line;
}
}

