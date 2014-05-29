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
$header=substr($line1,1);
@f=split(/\s+/,$header);
$readname1=$f[0];
$line1=<STDIN>;
chomp($line1);
$sequence1=$line1;
if(scalar(@f)==3){
$clr1=$f[1];
$clr2=$f[2];
}else{
$clr1=0;
$clr2=length($sequence1);
}
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
$sequence1 =~ tr/ACGTNacgtn/aaaaa99999/;# create fake quality scores
        print STDOUT "qlt:\n$sequence1\n.\n";
        print STDOUT "hps:\n.\n";
        print STDOUT "clv:$clr1,$clr2\n";
        print STDOUT "clr:$clr1,$clr2\n";
        print STDOUT "}\n";

}
}
