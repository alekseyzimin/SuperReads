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

#
#This program adds missing mates to the read fasta files "on the fly"

my $readnumberHold=-1;
my $prefixHold=0;
while($line=<STDIN>){
    chomp($line);
    @f=split(/\s+/,$line);
    $prefix=substr($f[0],1,2);
    $readnumber=substr($f[0],3);
    $editline="";
    for($i=1;$i<scalar(@f);$i++){
	$editline.="$f[$i] ";
	}

#print "DEBUG $prefix | $readnumber | $prefixHold | $readnumberHold\n";

    if($readnumber%2==0){#if the read is even we simply remember it
	if($readnumberHold!=-1){
	    print ">$prefixHold$readnumberHold $editlineHold\n$sequenceHold\n>$prefixHold",$readnumberHold+1,"\nN\n";
	}
	$prefixHold=$prefix;
	$readnumberHold=int($readnumber);
        $editlineHold=$editline;
	$line=<STDIN>;
	chomp($line);
	$sequenceHold=$line;
    }
    elsif($readnumberHold==-1){#the previous even read is missing
	print ">$prefix",$readnumber-1,"\nN\n$line\n";
	$line=<STDIN>;
	print $line;
    }
    elsif($readnumber-1!=$readnumberHold){#previous mate is missing odd and current is missing even
	print ">$prefixHold$readnumberHold $editlineHold\n$sequenceHold\n>$prefixHold",$readnumberHold+1,"\nN\n>$prefix",$readnumber-1,"\nN\n$line\n";
	$line=<STDIN>;
	print $line;
	$prefixHold="";
	$readnumberHold=-1;
    }
    elsif($readnumber-1==$readnumberHold){
	print ">$prefixHold$readnumberHold $editlineHold\n$sequenceHold\n$line\n";
	$line=<STDIN>;
	print $line;
	$prefixHold="";
	$readnumberHold=-1;
    }
    else{
	print "$line\n";
	die("error reading input file");
    }
}

