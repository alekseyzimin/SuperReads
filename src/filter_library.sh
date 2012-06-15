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
		$readnumber-- if($readnumber %2!=0);
		print "$prefix $readnumber $l[1] $l[2] $l[3]\n";
		}
}
}' $2.uidMuid $3  >layout.txt  

#figure out which jumping mates are non-junction (chimeric)
awk '{
print $5" "$3" "$1" "$2;
}' layout.txt |sort -k2,4 -S 10%|uniq -D -f 1 |awk '
BEGIN{flag=0}
{
if(flag==1){
	if($1>pos){
		if($1-pos<'$4')
			print $3""$4"\n"$3""$4+1;
	}else{
		if(pos-$1<'$4')
			print $3""$4"\n"$3""$4+1;
	}
}else{
	pos=$1;
}
flag=1-flag;
}'  >$2.chimeric.uid

#figure out redundant jumping mates, using both beginnings and ends of reads
awk '{print $1$2" "$3" "$5}' layout.txt |sort -k2,3 -S 20% |uniq -D -f 1 |awk '{print $2" "$3" "$1}' |sort -k3,3 -S 10%|uniq -D -f 2|awk '
BEGIN{flag=0}
{
if(flag==1){
	index1=int(substr(c1_1,1,length(c1_1)-1))*20000+c1_2;
	index2=int(substr($1,1,length($1)-1))*20000+$2;
	if(index1>index2){
		print c1_1" "$1" "c1_2" "$2" "c;
	}else{
		print $1" "c1_1" "$2" "c1_2" "c;
	}
}
c=$3;
c1_1=$1;
c1_2=$2;
flag=1-flag;
}'|perl -ane '{
chomp;
$range=1;
$code=0;
for($i=-$range;$i<=$range;$i++){
	for($j=-$range;$j<=$range;$j++){
		$code++ if(defined($h{"$F[0] $F[1] ".($F[2]+$i)." ".($F[3]+$j)}));
	}
}
if($code==0){
	$h{"$F[0] $F[1] $F[2] $F[3]"}=1;
}else{
print "$F[4]\n",substr($F[4],0,2),int(substr($F[4],2))+1,"\n";
}
}' > $2.redundant.uid

#figure out redundant jumping mates, using both beginnings and ends of reads
awk '{print $1$2" "$3" "$4}' layout.txt |sort -k2,3 -S 20% |uniq -D -f 1 |awk '{print $2" "$3" "$1}' |sort -k3,3 -S 10%|uniq -D -f 2|awk '
BEGIN{flag=0}
{
if(flag==1){ 
        index1=int(substr(c1_1,1,length(c1_1)-1))*20000+c1_2;
        index2=int(substr($1,1,length($1)-1))*20000+$2;
        if(index1>index2){
                print c1_1" "$1" "c1_2" "$2" "c;
        }else{
                print $1" "c1_1" "$2" "c1_2" "c;
        }
}
c=$3;
c1_1=$1;
c1_2=$2;
flag=1-flag;
}'|perl -ane '{
chomp;
$range=1;
$code=0;
for($i=-$range;$i<=$range;$i++){
        for($j=-$range;$j<=$range;$j++){
                $code++ if(defined($h{"$F[0] $F[1] ".($F[2]+$i)." ".($F[3]+$j)}));
        }
}
if($code==0){
        $h{"$F[0] $F[1] $F[2] $F[3]"}=1;
}else{
print "$F[4]\n",substr($F[4],0,2),int(substr($F[4],2))+1,"\n";
}
}' > $2.redundant.uid


echo -n "Non-junction reads "
wc -l $2.chimeric.uid
echo -n "Redundant reads "
wc -l $2.redundant.uid
