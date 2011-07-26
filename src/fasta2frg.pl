#!/usr/bin/perl

my $libId=$ARGV[0];

print STDOUT "{VER\n";
print STDOUT "ver:2\n";
print STDOUT "}\n";
print STDOUT "{LIB\n";
print STDOUT "act:A\n";
print STDOUT "acc:$libId\n";
print STDOUT "ori:I\n";
print STDOUT "mea:3000\n";
print STDOUT "std:300\n";
print STDOUT "src:\n";
print STDOUT ".\n";
print STDOUT "nft:1\n";
print STDOUT "fea:\n";
print STDOUT "doNotOverlapTrim=1\n";
print STDOUT ".\n";
print STDOUT "}\n";

while($line1=<STDIN>)
{
chomp($line1);
if($line1 =~ /^>/)
{
$readname1=substr($line1,1);
$line1=<STDIN>;
chomp($line1);
$sequence1=$line1;
next if(length($sequence1)<64);

        print STDOUT "{FRG\n";
        print STDOUT "act:A\n";
        print STDOUT "acc:$readname1\n";
        print STDOUT "rnd:0\n";
        print STDOUT "sta:G\n";
        print STDOUT "lib:$libId\n";
        print STDOUT "pla:0\n";
        print STDOUT "loc:0\n";
        print STDOUT "src:\n.\n";
        print STDOUT "seq:\n$sequence1\n.\n";
$sequence1 =~ tr/ACGTNacgtn/GGGGGGGGGG/;
        print STDOUT "qlt:\n$sequence1\n.\n";
        print STDOUT "hps:\n.\n";
        print STDOUT "clr:0,",length($sequence1),"\n";
        print STDOUT "}\n";

}
}
