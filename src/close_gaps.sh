#!/bin/sh
export LD_LIBRARY_PATH=/genome8/raid/alekseyz/test-SR-CA/MSR-CA-1.3.2/bin/../lib:$LD_LIBRARY_PATH
export PATH=/genome8/raid/alekseyz/test-SR-CA/MSR-CA-1.3.2/CA/Linux-amd64/bin/:/genome8/raid/alekseyz/test-SR-CA/MSR-CA-1.3.2/bin/:/genome8/raid/alekseyz/test-SR-CA/MSR-CA-1.3.2/bin/:$PATH

cat CA/9-terminator/genome.ctg.fasta | perl -ane 'if($F[0] =~ /^>/){chomp; print "\n$F[0]\n"}else{chomp;print}' | perl -ane 'BEGIN{$f=1}{if($f==0){print substr($F[0],4);}else{if(length($F[0])>100){print " ",substr($F[0],0,100)," ",substr($F[0],length($F[0])-100),"\n";}else{print " $F[0] $F[0]\n";}}$f=1-$f;}' >genome.ctg.fwd.pairs; 

awk '{print $1" "$5" "$2}' CA/9-terminator/genome.posmap.ctgscf |uniq -D -f 2 | awk 'BEGIN{prev_c="";prev_ori="";prev_sc=""}{if(prev_c ==""){prev_c=$1;prev_ori=$2;prev_sc=$3}else{if(prev_sc == $3){print prev_c" "prev_ori" "$1" "$2};prev_c=$1;prev_ori=$2;}prev_sc=$3}'| create_end_pairs.pl > contig_end_pairs.fa

awk '{print $1" "$5" "$2}' CA/9-terminator/genome.posmap.ctgscf |uniq -D -f 2 | awk 'BEGIN{prev_c="";prev_ori="";prev_sc=""}{if(prev_c ==""){prev_c=$1;prev_ori=$2;prev_sc=$3}else{if(prev_sc == $3){print prev_c" "prev_ori" "$1" "$2};prev_c=$1;prev_ori=$2;}prev_sc=$3}' | awk 'BEGIN{rn=0;}{print "cc"rn" cc"++rn" "$0;rn++}' > read_pair_contig_pair.txt

jellyfish count -m 31 -t 16 -C -r -s 200000000 -o k_u_hash pe.cor.fa sj.cor.fa
create_k_unitigs --cont-on-low  --low-stretch=30 -C -t 16 -o close_gap_k_u -m 2 -M 2  k_u_hash_0
rm close_gap_k_u_*.counts
cat close_gap_k_u_*.fa | awk 'BEGIN{n=0}{if($1~/^>/){print ">"n" "$2" "$3" "$4;n++;}else{print $0}}' > close_gap_k_u.all.fa
rm close_gap_k_u_*.fa
createSuperReadsForDirectory.perl  -mikedebug -force-join -default-mean 500 -default-stdev 200 -join-shooting  -minreadsinsuperread 1 -kunitigsfile close_gap_k_u.all.fa -join-mates -l 31 -s 200000000 -t 16 -mkudisr 0 work4 contig_end_pairs.fa 1>super.out 2>super.out

awk '{if(int(substr($1,3))%2==0){print $3" "$2" "$1;}else{print $3" "$2" "substr($1,1,2)""int(substr($1,3))-1}}' work4/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt|uniq -d -f 1 | awk '{print $3}' > closed_gap_pairs.txt

echo -n "closed gaps ";wc -l closed_gap_pairs.txt

nucmer -l 15 -c 21 -b 2 contig_end_pairs.fa /genome9/raid/alekseyz/GAGE/rhodobacter/max_contig_size/original/all.fa 
echo -n "closed gaps with reads "; show-coords -lcH out.delta | awk '{rn=int(substr($18,3));if(rn%2==0)print $0" "rn" "rn" "$19; else print $0" "rn" "rn-1" "$19}' | sort -k21,22|uniq -D -f 20 |uniq -u -f 19 > gaps_closed.coords;awk '{print "cc"$(NF-1)}'  gaps_closed.coords |uniq -c |awk '{if($1>3)print $2}' >closed_gap_pairs.reads.txt;wc -l closed_gap_pairs.reads.txt

cat closed_gap_pairs.txt closed_gap_pairs.reads.txt >closed_gap_pairs.all.txt

cat CA/9-terminator/genome.ctg.fasta| perl -ane 'if($F[0] =~ /^>/){chomp; print "\n",substr($F[0],4)," "}else{chomp;print}' |perl -e '{open(FILE,"read_pair_contig_pair.txt");while($line=<FILE>){chomp($line);@l=split(/\s+/,$line);$ctg{$l[0]}=$l[2];$ctg{$l[1]}=$l[4];}open(FILE,"closed_gap_pairs.all.txt");while($line=<FILE>){chomp($line);$closed{$ctg{$line}}=1;}$seq="";$name="";while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);next if($l[0] eq "");if(defined($closed{$l[0]})){$name.="$l[0]_";$seq.=$l[1];}else{if(not($name eq "")){print ">$name$l[0]\n$seq$l[1]\n";$name="";$seq="";}else{print ">$l[0]\n$l[1]\n";}}}}' > merged_ctg.fa

~/myprogs/compute_n50_stats.pl  merged_ctg.fa
