#!/bin/bash
#$1 -- PREFIX
rm -f f_*
rm -f *.out
rm -f *.broken

grep -H FAILED `ls ../5-consensus/${1}_*.err ../5-consensus/${1}_*.success |awk -F '.' '{print $3}' |sort |uniq -u | awk '{print ".."$1".err"}'` > unitig_failures
awk '{split($1,a,"_"); print int(a[2])" "$3}' unitig_failures|uniq| awk '{print "tigStore -g ../'$1'.gkpStore -t ../'$1'.tigStore 1 -up "$1" -d layout -u "$2" > f_unitig"$2}' > extract_layouts.sh
bash ./extract_layouts.sh

awk 'BEGIN{unum=0}{if($3 !=unum){print $3" "substr($9,1,length($9)-1); unum=$3}}' unitig_failures > unitig_split_points
perl -e '{
open(FILE,"unitig_split_points"); 
while($line=<FILE>){
	chomp($line);
	@f=split(/\s+/,$line);
	$h{$f[0]}=$f[1];
	}
foreach $v (keys %h){
	$count1=0;
	$count2=0;
	open(FILE,"f_unitig$v");
	while($line=<FILE>){
		if($line=~/^data.unitig_coverage_stat/){
		chomp($line);
		@f1=split(/\s+/,$line); 
		$cov=$f1[1];
		next
		};
		next if(not($line =~ /^FRG/)); 
		@f2=split(/\s+/,$line);
		if($f2[4] eq $h{$v}){
			$count2++;
		}elsif($count2==0){
			$count1++;
		}else{
			$count2++;
		}
		}
	close(FILE);
	open(FILE,"f_unitig$v");
	open(OUTFILE,">f_unitig$v.fixed");
	print OUTFILE "\nunitig $v\nlen 0\ncns\nqlt\ndata.unitig_coverage_stat $cov\ndata.unitig_microhet_prob 1.000000\ndata.unitig_status        X\ndata.unitig_unique_rept   X\ndata.contig_status        U\ndata.num_frags            $count1\ndata.num_unitigs          0\n";
	$count=-1;
	while($line=<FILE>){
		next if(not($line =~ /^FRG/));
		$count ++;
		if($count ==$count1){
		print OUTFILE "\nunitig -1\nlen 0\ncns\nqlt\ndata.unitig_coverage_stat $cov\ndata.unitig_microhet_prob 1.000000\ndata.unitig_status        X\ndata.unitig_unique_rept   X\ndata.contig_status        U\ndata.num_frags            $count2\ndata.num_unitigs          0\n";
		} 
	print OUTFILE $line;
	}close(OUTFILE);
}
}' 

awk '{split($1,a,"_"); print int(a[2])" "$3}' unitig_failures|uniq| awk '{print "tigStore -g ../'$1'.gkpStore -t ../'$1'.tigStore 1 -up "$1" -R f_unitig"$2".fixed.final"}' > replace_layouts.sh
for f in $(ls *.fixed);do
  awk 'BEGIN{flag=0}{if($1 ~ /^unitig/){if($2==-1){flag=0;} else{flag=1}} if(flag==1){print $0}}' $f > ${f}.final
done
bash ./replace_layouts.sh

awk '{split($1,a,"_"); print int(a[2])" "$3}' unitig_failures|uniq| awk '{print "utgcns -g ../'$1'.gkpStore -t ../'$1'.tigStore 1 "$1}' |uniq > recompute_consensus.sh
bash ./recompute_consensus.sh 1>utgcns.err 2>&1

cat *.fixed | awk 'BEGIN{flag=0}{if($1 ~ /^unitig/){if($2==-1){flag=1;} else{flag=0;}} if(flag==1 && $0 ~ /^FRG/){print "frg iid "$5" mateiid 0"}}' > gkp.edits.msg
gatekeeper --edit gkp.edits.msg ../$1.gkpStore 1>gkp.out 2>&1
awk '{if($5>0){print "frg iid "$5" mateiid 0"}}' gkp.out > gkp.edits1.msg
gatekeeper --edit gkp.edits1.msg ../$1.gkpStore 1>gkp1.out 2>&1

touch ..//5-consensus/consensus.success
