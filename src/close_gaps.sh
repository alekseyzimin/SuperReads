#!/bin/bash
cat CA/9-terminator/genome.ctg.fasta| perl -ane 'if($F[0] =~ /^>/){chomp; print "\n",substr($F[0],4)," "}else{chomp;print}' | awk '{print $1" "substr($2,1,150)" "substr($2,length($2)-150)}' >genome.ctg.fwd.pairs; 
cat CA/9-terminator/genome.ctg.fasta| reverse_complement | perl -ane 'if($F[0] =~ /^>/){chomp; print "\n",substr($F[0],4)," "}else{chomp;print}' | awk '{print $1" "substr($2,1,150)" "substr($2,length($2)-150)}' >genome.ctg.rev.pairs; 

awk '{print $1" "$5" "$2}' CA/9-terminator/genome.posmap.ctgscf |uniq -D -f 2 | awk 'BEGIN{prev_c="";prev_ori="";prev_sc=""}{if(prev_c ==""){prev_c=$1;prev_ori=$2;prev_sc=$3}else{if(prev_sc == $3){print prev_c" "prev_ori" "$1" "$2};prev_c=$1;prev_ori=$2;}prev_sc=$3}'| perl -e '{open(FILE,"genome.ctg.fwd.pairs");while($line=<FILE>){chomp($line);@l=split(/\s+/,$line);$fwd_b{$l[0]}=$l[1];;$fwd_e{$l[0]}=$l[2];}close(FILE);open(FILE,"genome.ctg.rev.pairs");while($line=<FILE>){chomp($line);@l=split(/\s+/,$line);$rev_b{$l[0]}=$l[1];;$rev_e{$l[0]}=$l[2];}close(FILE);$prefix="cc";$readnumber=0;while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);if($l[1] eq "f"){print ">$prefix",$readnumber,"\n$fwd_e{$l[0]}\n";}else{print ">$prefix",$readnumber,"\n$rev_e{$l[0]}\n";}$readnumber++;if($l[3] eq "f"){print ">$prefix",$readnumber,"\n$fwd_b{$l[0]}\n";}else{print ">$prefix",$readnumber,"\n$rev_b{$l[0]}\n";}$readnumber++;}}' > contig_end_pairs.fa

awk '{print $1" "$5" "$2}' CA/9-terminator/genome.posmap.ctgscf |uniq -D -f 2 | awk 'BEGIN{prev_c="";prev_ori="";prev_sc=""}{if(prev_c ==""){prev_c=$1;prev_ori=$2;prev_sc=$3}else{if(prev_sc == $3){print prev_c" "prev_ori" "$1" "$2};prev_c=$1;prev_ori=$2;}prev_sc=$3}' | awk 'BEGIN{rn=0;}{print "cc"rn" cc"++rn" "$0;rn++}' > read_pair_contig_pair.txt

create_k_unitigs -C -t 48 -o close_gap_k_u -m 2 -M 2 --cont-on-low --low-stretch=21 k_u_hash_0
rm close_gap_k_u_*.counts
cat close_gap_k_u_*.fa > close_gap_k_u.all.fa
rm close_gap_k_u_*.fa
createSuperReadsForDirectory.perl -minreadsinsuperread 1 -kunitigsfile close_gap_k_u.all.fa -join-mates -l 31 -s 200000000 -t 16 -mkudisr 0 work4 contig_end_pairs.fa 1>super.out 2>super.err
awk '{print $1" "$2}'  work4/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt |sort -k2,2|uniq -d -f 1 | awk '{print $1}' > closed_gap_pairs.txt

echo -n "closed gaps ";wc -l closed_gap_pairs.txt

cat CA/9-terminator/genome.ctg.fasta| perl -ane 'if($F[0] =~ /^>/){chomp; print "\n",substr($F[0],4)," "}else{chomp;print}' |perl -e '{open(FILE,"read_pair_contig_pair.txt");while($line=<FILE>){chomp($line);@l=split(/\s+/,$line);$ctg{$l[0]}=$l[2];$ctg{$l[1]}=$l[4];}open(FILE,"closed_gap_pairs.txt");while($line=<FILE>){chomp($line);$closed{$ctg{$line}}=1;}$seq="";$name="";while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);next if($l[0] eq "");if(defined($closed{$l[0]})){$name.="$l[0]_";$seq.=$l[1];}else{if(not($name eq "")){print ">$name$l[0]\n$seq$l[1]\n";$name="";$seq="";}else{print ">$l[0]\n$l[1]\n";}}}}' > merged_ctg.fa

