#!/bin/bash
#$1 path to the stores
#$2 prefix of the stores
#$3 file with ALL jump library UIDs
#$4 threshold for chimeric sizes (750)

echo "Filtering libraries for chimerism...";

#dump all read uids
gatekeeper -dumpfragments -tabular $1/$2.gkpStore |awk '{print $1" "$3}' > $2.uidMuid

#dump tigStore and create layout of unitigs for jumping libraries ONLY
tigStore -g $1/$2.gkpStore -t $1/$2.tigStore 1 -U -d layout | awk '{
if($1 ~ /^u/){
	unitig=$2;
}else if($1 ~ /^F/){
	print $5" "unitig" "$(NF-1)" "$NF;
	}	
}' | perl -e '{
open(FILE,$ARGV[1]);
while($line=<FILE>){
	chomp($line);
	$h{$line}=1;
}

open(FILE, $ARGV[0]);
while($line=<FILE>){
	@l=split(/\s+/,$line);
	push(@uid,$l[0]);
}

while($line=<STDIN>){
	chomp($line);
	@l=split(/\s+/,$line);
	if(defined($h{$uid[$l[0]]})){
		$prefix=substr($uid[$l[0]],0,2); 
		$readnumber=int(substr($uid[$l[0]],2)); 
                $pairnumber=$readnumber;
                $pairnumber-- if($pairnumber%2); 
                if($l[2]>$l[3]){
		$ori="R";
		}else{
                $ori="F";
		}
		print "$prefix$readnumber $l[1] $l[2] $ori $prefix$pairnumber\n";
		}
}
}' $2.uidMuid $3  >layout.txt  

#figure out which jumping mates are non-junction 
sort -k5,5 -S 10% layout.txt | filter_alt.pl innie > genome.chimeric.uid

