#!/usr/bin/env perl
#
#This program adds missing mates to the read fasta files "on the fly"

my $readnumberHold=-1;
my $prefixHold=0;
while($line=<STDIN>){
    chomp($line);
    @f=split(/\s+/,$line);
    $prefix=substr($f[0],1,2);
    $readnumber=substr($f[0],3);

#print "DEBUG $prefix | $readnumber | $prefixHold | $readnumberHold\n";

    if($readnumber%2==0){#if the read is even we simply remember it
	if($readnumberHold!=-1){
	    print ">$prefixHold",$readnumberHold,"\n$sequenceHold\n>$prefixHold",$readnumberHold+1,"\nN\n";
	}
	$prefixHold=$prefix;
	$readnumberHold=int($readnumber);
	$line=<STDIN>;
	chomp($line);
	$sequenceHold=$line;
    }
    elsif($readnumberHold==-1){#the prevois even read is missing
	print ">$prefix",$readnumber-1,"\nN\n$line\n";
	$line=<STDIN>;
	print $line;
    }
    elsif($readnumber-1!=$readnumberHold){#previous mate is missing odd and current is missing even
	print ">$prefixHold$readnumberHold\n$sequenceHold\n>$prefixHold",$readnumberHold+1,"\nN\n>$prefix",$readnumber-1,"\nN\n$line\n";
	$line=<STDIN>;
	print $line;
	$prefixHold="";
	$readnumberHold=-1;
    }
    elsif($readnumber-1==$readnumberHold){
	print ">$prefixHold$readnumberHold\n$sequenceHold\n$line\n";
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

