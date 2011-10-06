#!/usr/bin/env perl
#This script reads the the config file and runs the full super reads pipeline
#config file format:
#
#PATHS
#CA_PATH= path to runCA
#JELLYFISH_PATH= path to jellyfish
#SR_PATH=  path to super reads routine installation
#END
#
#DATA
#PE= id mean stdev fastq_forward fastq_reverse
#JUMP= id mean stdev fastq_forward fastq_reverse 
#OTHER=frag1.frg
#END
#
#PARAMETERS
#WINDOW=10
#MAX_ERR_PER_WINDOW=3
#TRIM_PARAM=2
#NUM_THREADS=24
#JF_SIZE=2000000000
#EXTEND_JUMP_READS=0
#END
#

#parsing config file
my $JELLYFISH_PATH="";
my $SR_PATH="";
my $CA_PATH="";
my $CA_PARAMETERS=" cgwErrorRate=0.15 utgErrorRate=0.03 merylMemory=8192 ovlMemory=4GB ";
my $WINDOW=10;
my $MAX_ERR_PER_WINDOW=3;
my $TRIM_PARAM=2;
my $NUM_THREADS=2;
my $EXTEND_JUMP_READS=0;
my $JF_SIZE=100000000;
my $in_paths=0;
my $in_data=0;
my $in_parameters=0;
my @pe_info_array=();
my @jump_info_array=();
my @other_info_array=();
my %used_library_ids;


open(FILE,$ARGV[0]);
while($line=<FILE>)
{
chomp($line);
if($line =~ /^PATHS/)
{
if($in_parameters==1)
{
die("error in config file: mixed PARAMETERS and PATHS");
}
if($in_data==1)
{
die("error in config file: mixed PATHS and DATA");
}
if($in_paths==1)
{
die("duplicate PATHS header"); 
}
$in_paths=1;
next;
}
elsif($line =~ /^DATA/)
{
if($in_parameters==1)
{
die("error in config file: mixed PARAMETERS and DATA");
}
if($in_paths==1)
{
die("error in config file: mixed PATHS and DATA");
}
if($in_data==1)
{
die("duplicate PATHS header");
}
$in_data=1;
next;
}
elsif($line =~ /^PARAMETERS/)
{
if($in_data==1)
{
die("error in config file: mixed PARAMETERS and DATA");
}
if($in_paths==1)
{
die("error in config file: mixed PARAMETERS and PATHS");
}
if($in_parameters==1)
{
die("duplicate PARAMETERS header");
}
$in_parameters=1;
next;
}

if($in_parameters==1)
{
if($line =~ /^END/)
{
$in_parameters=0;
next;
}
elsif($line =~ /^EXTEND_JUMP_READS/)
{
@f=split(/=/,$line);
$f[1]=~s/^\s+//;
$f[1]=~s/\s+$//;
$EXTEND_JUMP_READS=int($f[1]);
die("bad value for EXTEND_JUMP_READS") if($EXTEND_JUMP_READS!=1&&$EXTEND_JUMP_READS!=0);
next;
}
elsif($line =~ /^WINDOW/)
{
@f=split(/=/,$line);
$f[1]=~s/^\s+//;
$f[1]=~s/\s+$//;
$WINDOW=int($f[1]);
die("bad value for WINDOW") if($WINDOW<1);
next;
}
elsif($line =~ /^CA_PARAMETERS/)
{
@f=split(/=/,$line);
$CA_PARAMETERS=$f[1];
for($i=2;$i<scalar(@f);$i++)
{

$CA_PARAMETERS.="=$f[$i]";
}
die("bad value for CA_PARAMETERS") if($CA_PARAMETERS eq "");
next;
}
elsif($line =~ /^MAX_ERR_PER_WINDOW/)
{
@f=split(/=/,$line);
$f[1]=~s/^\s+//;
$f[1]=~s/\s+$//;
$MAX_ERR_PER_WINDOW=int($f[1]);
die("bad value for MAX_ERR_PER_WINDOW") if($MAX_ERR_PER_WINDOW<1);
next;
}
elsif($line =~ /^TRIM_PARAM/)
{
@f=split(/=/,$line);
$f[1]=~s/^\s+//;
$f[1]=~s/\s+$//;
$TRIM_PARAM=int($f[1]);
die("bad value for TRIM_PARAM") if($TRIM_PARAM<1);
next;
}
elsif($line =~ /^NUM_THREADS/)
{
@f=split(/=/,$line);
$f[1]=~s/^\s+//;
$f[1]=~s/\s+$//;
$NUM_THREADS=int($f[1]);
die("bad value for NUM_THREADS") if($NUM_THREADS<1);
next;
}
elsif($line =~ /^JF_SIZE/)
{
@f=split(/=/,$line);
$f[1]=~s/^\s+//;
$f[1]=~s/\s+$//;
$JF_SIZE=int($f[1]);
die("bad value for JF_SIZE") if($JF_SIZE<100000);
next;
}
else
{
next;
}
}


if($in_paths==1)
{
if($line =~ /^END/)
{
$in_paths=0;
next;
}
elsif($line =~ /^JELLYFISH_PATH/)
{
@f=split(/=/,$line);
$JELLYFISH_PATH=$f[1];
$JELLYFISH_PATH=~s/^\s+//;
$JELLYFISH_PATH=~s/\s+$//;
$JELLYFISH_PATH.="/" if(not($JELLYFISH_PATH eq ""));
next;
}
elsif($line =~ /^CA_PATH/)
{
@f=split(/=/,$line);
$CA_PATH=$f[1];
$CA_PATH=~s/^\s+//;
$CA_PATH=~s/\s+$//;
$CA_PATH.="/" if(not($JELLYFISH_PATH eq ""));
next;
}
elsif($line =~ /^SR_PATH/)
{
@f=split(/=/,$line);
$SR_PATH=$f[1];
$SR_PATH=~s/^\s+//;
$SR_PATH=~s/\s+$//;
$SR_PATH.="/" if(not($SR_PATH eq ""));
next;
}
else
{
next;
}
}
elsif($in_data==1)
{
if($line =~ /^END/)
{
if(scalar(@pe_info_array)==0 && scalar(@jump_info_array)==0 && scalar(@other_info_array)==0)
{
die("incomplete DATA section");
}
$in_data=0;
next;
}
elsif($line =~ /^PE/)
{
$pe_info_line="";
@f=split(/=/,$line);
$f[1]=~s/^\s+//;
@f1=split(/\s+/,$f[1]);
die("improper id for PE library $f1[0]") if(not(length($f1[0])==2));
die("duplicate id for PE library $f1[0]") if(defined($used_library_ids{$f1[0]}));
$pe_info_line.="$f1[0] ";
$used_library_ids{$f1[0]}=1;
die("improper mean for PE library $f1[0]") if(not(int($f1[1])>0));
$pe_info_line.="$f1[1] ";
die("improper stdev for PE library $f1[0]") if(not(int($f1[2])>0));
$pe_info_line.="$f1[2] ";
die("missing forward file for PE library $f1[0]") if(not(-e $f1[3]));
$pe_info_line.="$f1[3] ";
if(length($f1[4])>0)
{
die("missing reverse file for PE library $f1[0]") if(not(-e $f1[4]));
$pe_info_line.="$f1[4] ";
}
else
{
$pe_info_line.="$f1[3] ";
}
push(@pe_info_array,$pe_info_line);
}
elsif($line =~ /^JUMP/)
{
$jump_info_line="";
@f=split(/=/,$line);
$f[1]=~s/^\s+//;
@f1=split(/\s+/,$f[1]);
die("improper id for JUMP library $f1[0]") if(not(length($f1[0])==2));
die("duplicate id for JUMP library $f1[0]") if(defined($used_library_ids{$f1[0]}));
$jump_info_line.="$f1[0] ";
$used_library_ids{$f1[0]}=1;
die("improper mean for JUMP library $f1[0]") if(not(int($f1[1])>0));
$jump_info_line.="$f1[1] ";
die("improper stdev for JUMP library $f1[0]") if(not(int($f1[2])>0));
$jump_info_line.="$f1[2] ";
die("missing forward file for JUMP library $f1[0]") if(not(-e $f1[3]));
$jump_info_line.="$f1[3] ";
die("missing reverse file for JUMP library $f1[0]") if(not(-e $f1[4]));
$jump_info_line.="$f1[4] ";
push(@jump_info_array,$jump_info_line);
}
elsif($line=~ /OTHER/)
{
@f=split(/=/,$line);
$f[1]=~s/^\s+//;
$f[1]=~s/\s+$//;
if($f[1] =~/\.frg$/ && (-e $f[1]))
{
push(@other_info_array,$f[1]);
}
else
{
die("missing or incorrect frg file for OTHER library");
}
}
}
}
die("no read data files specified") if(scalar(@pe_info_array)+scalar(@jump_info_array)+scalar(@other_info_array)==0);

print "PATHS:\nCA_PATH = $CA_PATH\nJELLYFISH_PATH = $JELLYFISH_PATH\nSR_PATH = $SR_PATH\n\n";
print "Verifying PATHS...\n";
if(-e "${JELLYFISH_PATH}jellyfish")
{
print "jellyfish ok\n";
}
else
{
die("jellyfish not found at ${JELLYFISH_PATH}");
}
if(-e "${CA_PATH}runCA")
{
print "runCA ok\n";
}
else
{
die("runCA not found at ${CA_PATH}");
}
if(-e "${SR_PATH}createSuperReadsForDirectory.perl")
{
print "createSuperReadsForDirectory.perl ok\n";
}
else
{
die("createSuperReadsForDirectory.perl not found at ${SR_PATH}");
}

print "\nPARAMETERS:\nWINDOW = $WINDOW\nMAX_ERR_PER_WINDOW = $MAX_ERR_PER_WINDOW\nTRIM_PARAM = $TRIM_PARAM\nNUM_THREADS = $NUM_THREADS\nJF_SIZE = $JF_SIZE\nEXTEND_JUMP_READS = $EXTEND_JUMP_READS\nCA_PARAMETERS = $CA_PARAMETERS\n\n";
print "DATA:\nPE:\n";
foreach $v(@pe_info_array)
{
print "$v\n";
}
print "JUMP:\n";
foreach $v(@jump_info_array)
{
print "$v\n";
}
print "OTHER:\n";
foreach $v(@other_info_array)
{
print "$v\n";
}

die ("no PE data specified") if(scalar(@pe_info_array)==0);

print("\ncreating script file for the actions...");

open(FILE,">assemble.sh");
print FILE "#!/bin/bash\n\n";
print FILE "export LD_LIBRARY_PATH=$JELLYFISH_PATH../lib:\$LD_LIBRARY_PATH\n";
print FILE "export PATH=$CA_PATH:$SR_PATH:$JELLYFISH_PATH:\$PATH\n";
$i=0;
$list_pe_files="";
$list_jump_files="";
##################################################renaming reads####################################################
print FILE "echo -n 'processing PE library reads ';date;\n";

foreach $v(@pe_info_array)
{
@f=split(/\s+/,$v);
$list_pe_files.="$f[0].renamed.fastq ";
next if(-e "$f[0].renamed.fastq");
if($f[3] eq $f[4])
{
print FILE "cat $f[3] | perl -e '{\$library=\$ARGV[0];\$readnumber=0;while(\$line=<STDIN>){if(\$line=~ /^@/){\$line=<STDIN>;chomp(\$line);\@seq=split(/\\s+/,\$line);\$line=<STDIN>;\$line=<STDIN>;\@qlt=split(/\\s+/,\$line);print \"\@\",\"\$library\$readnumber\\n\$seq[0]\\n+\\n\$qlt[0]\\n\";\$readnumber+=2;}}}' $f[0] > $f[0].renamed.fastq &\nPID$i=\$!\n";
}
else
{
print FILE "paste $f[3] $f[4] | perl -e '{\$library=\$ARGV[0];\$readnumber=0;while(\$line=<STDIN>){if(\$line=~ /^@/){\$line=<STDIN>;chomp(\$line);\@seq=split(/\\s+/,\$line);\$line=<STDIN>;\$line=<STDIN>;\@qlt=split(/\\s+/,\$line);print \"\@\",\"\$library\$readnumber\\n\$seq[0]\\n+\\n\$qlt[0]\\n\";\$readnumber++;print \"\@\",\"\$library\$readnumber\\n\$seq[1]\\n+\\n\$qlt[1]\\n\";\$readnumber++;}}}' $f[0] > $f[0].renamed.fastq &\nPID$i=\$!\n";
}
$i++;
}
print FILE "wait ";
for(my $j=0;$j<$i;$j++)
{
print FILE "\$PID$j ";
}
print FILE "\n";

foreach $v(@pe_info_array)
{
@f=split(/\s+/,$v);
print FILE "MAX_$f[0]=`tail $f[0].renamed.fastq |grep '^\@$f[0]' |tail -n 1|awk '{print substr(\$1,4)}'`\n"
}


if(scalar(@jump_info_array)>0)
{
print FILE "echo -n 'processing JUMP library reads ';date;\n";
$i=0;
foreach $v(@jump_info_array)
{
@f=split(/\s+/,$v);
$list_jump_files.="$f[0].renamed.fastq ";
next if(-e "$f[0].renamed.fastq");
if($f[3] eq $f[4])
{
die("duplicate jump library $f[0] files");
}
else
{
print FILE "paste $f[3] $f[4] | perl -e '{\$library=\$ARGV[0];\$readnumber=0;while(\$line=<STDIN>){if(\$line=~ /^@/){\$line=<STDIN>;chomp(\$line);\@seq=split(/\\s+/,\$line);\$line=<STDIN>;\$line=<STDIN>;\@qlt=split(/\\s+/,\$line);print \"\@\",\"\$library\$readnumber\\n\$seq[0]\\n+\\n\$qlt[0]\\n\";\$readnumber++;print \"\@\",\"\$library\$readnumber\\n\$seq[1]\\n+\\n\$qlt[1]\\n\";\$readnumber++;}}}' $f[0] > $f[0].renamed.fastq &\nPID$i=\$!\n";
}
$i++;
}
print FILE "wait ";
for(my $j=0;$j<$i;$j++)
{
print FILE "\$PID$j ";
}
print FILE "\n";

foreach $v(@jump_info_array)
{
@f=split(/\s+/,$v);
print FILE "MAX_$f[0]=`tail $f[0].renamed.fastq |grep '^\@$f[0]' |tail -n 1|awk '{print substr(\$1,4)}'`\n"
}
}
######################################################done renaming reads####################################################################

########################################################Jellyfish#################################

print FILE "echo -n 'running Jellyfish ';date;\n";
print FILE "Q2_CHAR=`cat $list_pe_files |head -n 4000 | awk 'BEGIN{flag=0}{if(\$0 ~ /^\+/){flag=1}else if(flag==1){print \$0;flag=0}}'  | perl -ne 'BEGIN{\$q0_char=\"\@\";}{chomp;\@f=split \"\";foreach \$v(\@f){if(ord(\$v)<ord(\$q0_char)){\$q0_char=\$v;}}}END{print chr(ord(\$q0_char)+2)}'`\n";
print FILE "echo Q2_CHAR: \$Q2_CHAR\n";

print FILE "\ncat $list_pe_files | filter_illumina_data_quals '\$Q2_CHAR'  > pe_trim.fa\n" if(not(-e "pe_trim.fa"));

#for our error correctin we run jellyfish twice: once on all pe bases and once on pe bases with quality >2 
print FILE "jellyfish count -t $NUM_THREADS -C -r -o pe_trim -s $JF_SIZE -m 24 pe_trim.fa\n"  if(not(-e "pe_trim_0"));
print FILE "jellyfish count -t $NUM_THREADS -C -r -o pe_all -s $JF_SIZE -m 24 $list_pe_files\n"  if(not(-e "pe_all_0"));

#check if the JF_SIZE was big enough:  we want to end up with a single raw database for pe_all and pe_trim
print FILE "if [[ -e pe_trim_1 || -e pe_all_1 ]];then\n";
print FILE "echo \"Increase JF_SIZE\"\n";
print FILE "exit\n";
print FILE "fi\n";
########################################################done Jellyfish#################################
#initialize list of frg files
$list_of_frg_files="";
print FILE "\n\n\n";
######################################################error correct PE#############################################################
if(not(-e "pe.cor.fa"))
{
print FILE "echo -n 'error correct PE ';date;\n";
print FILE "\nerror_correct_reads -d pe_trim_0 -d pe_all_0 -C -m 1 -s 1 -g 1 -a 3 -t $NUM_THREADS -w $WINDOW -e $MAX_ERR_PER_WINDOW $list_pe_files\n";
print FILE "cat error_corrected_*.fa  | homo_trim $TRIM_PARAM > pe.cor.fa\n";
print FILE "rm error_corrected_*\n";
}
#compute average PE read length -- we will need this for Astat later
print FILE "PE_AVG_READ_LENGTH=`head -n 1000000 pe.cor.fa |tail -n 500000| grep -v '^>' | awk 'BEGIN{n=0;m=0;}{m+=length(\$0);n++;}END{print int(m/n)}'`\n";
print FILE "echo \"Average PE read length after error correction: \$PE_AVG_READ_LENGTH\"\n"; 
######################################################done error correcct PE##########################################################
print FILE "\n";
########################################################error correct, super reads for JUMP#################################
if(scalar(@jump_info_array)>0)
{
if(not(-e "sj.cor.fa"))
{
print FILE "echo -n 'error correct JUMP ';date;\n";
print FILE "\nerror_correct_reads -d pe_trim_0 -d pe_all_0 -C -m 1 -s 1 -g 2 -a 3 -t $NUM_THREADS -w $WINDOW -e $MAX_ERR_PER_WINDOW $list_jump_files\n";
print FILE "cat error_corrected_*.fa  | homo_trim $TRIM_PARAM > sj.cor.fa\n";
print FILE "rm error_corrected_*\n";
}
#############################################################done error correct JUMP#############################################
#
print FILE "\n\n\n";
#######################################################build k-unitigs##############################################################
print FILE "\n";
if(not(-e "guillaumeKUnitigsAtLeast32bases_all.fasta"))
{
print FILE "jellyfish count -m 31 -t $NUM_THREADS -C -r -s $JF_SIZE -o k_u_hash pe.cor.fa sj.cor.fa\n";
print FILE "create_k_unitigs -C -t $NUM_THREADS  -m 2 -M 2 -l 31 -o k_unitigs k_u_hash_0 1> /dev/null 2>&1\n";
print FILE "cat k_unitigs_*.fa > guillaumeKUnitigsAtLeast32bases_all.fasta\n";
print FILE "rm k_unitigs_*.fa  k_unitigs_*.counts\n";
}
}
else
{
if(not(-e "guillaumeKUnitigsAtLeast32bases_all.fasta"))
{
print FILE "jellyfish count -m 31 -t $NUM_THREADS -C -r -s $JF_SIZE -o k_u_hash pe.cor.fa\n";
print FILE "create_k_unitigs -C -t $NUM_THREADS  -m 2 -M 2 -l 31 -o k_unitigs k_u_hash_0 1> /dev/null 2>&1\n";
print FILE "cat k_unitigs_*.fa > guillaumeKUnitigsAtLeast32bases_all.fasta\n";
print FILE "rm k_unitigs_*.fa  k_unitigs_*.counts\n";
}
}
print FILE "ESTIMATED_GENOME_SIZE=`perl -ane '{if(\$F[0]=~/^>/){print \"\\n\"}else{print \$F[0]}}'  guillaumeKUnitigsAtLeast32bases_all.fasta | wc|awk '{print int(\$3)-(30*(int(\$1)-1));}'`\n";
print FILE "echo \"Estimated genome size: \$ESTIMATED_GENOME_SIZE\"\n";
######################################################done k-unitigs#################################################################
print FILE "\n\n\n";
########################################super reads for jump#######################################################################
if(scalar(@jump_info_array)>0)
{
print FILE "echo -n 'filtering JUMP ';date;\n";

#creating super reads. for filtering
print FILE "createSuperReadsForDirectory.perl -minreadsinsuperread 2 -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.fasta -l 31 -s $JF_SIZE -t $NUM_THREADS -M 2 -m 2 -join-mates -join-shooting -mkudisr 0 work2 sj.cor.fa 1> super2.err 2>&1\n" if(not(-e "work2"));;
print FILE "\n";

#now, using read positions in super reads, we find out which mates got joined -- these are the ones that do not have the biotin in the middle, call them chimeric
if(not(-e "chimeric_sj.txt"))
{
print FILE "awk '{if(int(substr(\$1,3))%2==0){print \$3\" \"\$2\" \"\$1;}else{print \$3\" \"\$2\" \"substr(\$1,1,2)\"\"int(substr(\$1,3))-1}}' work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt |uniq -D -f 1| awk 'BEGIN{insert=\"\";}{if(\$3!=insert){start=\$1;insert=\$3}else{if(start>\$1){print insert\" \"start-\$1}else{print insert\" \"\$1-start}}}' | perl -ane '{if(\$F[1]<750&&\$F[1]>0){print STDOUT \"\$F[0]\\n\",substr(\$F[0],0,2),int(substr(\$F[0],2))+1,\"\\n\";}}' 1> chimeric_sj.txt \n" if(not(-e "chimeric_sj.txt"));
print FILE "\n";
print FILE "awk '{if(int(substr(\$1,3))%2==0){print \$4\" \"\$2\" \"\$1;}else{print \$4\" \"\$2\" \"substr(\$1,1,2)\"\"int(substr(\$1,3))-1}}' work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt |uniq -d | perl -ane '{print STDOUT \"\$F[2]\\n\",substr(\$F[2],0,2),int(substr(\$F[2],2))+1,\"\\n\";}' 1>> chimeric_sj.txt \n";
print FILE "\n";
}
#we also do initial redundancy filtering here, based on positions of reads in suoer reads
print FILE "cat work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt | awk '{if(int(substr(\$1,3))%2==0){print \$2\" \"\$3\" \"\$1}else{print \$2\" \"\$3\" \"substr(\$1,1,2)\"\"int(substr(\$1,3))-1}}' |uniq -D -f 2|awk 'BEGIN{flag=0}{if(flag==1){index1=int(substr(c1_1,1,length(c1_1)-1))*20000+c1_2;index2=int(substr(\$1,1,length(\$1)-1))*20000+\$2;if(index1>index2){print c1_1\" \"\$1\" \"c1_2\" \"\$2\" \"c}else{print \$1\" \"c1_1\" \"\$2\" \"c1_2\" \"c}}c=\$3;c1_1=\$1;c1_2=\$2;flag=1-flag;}'|perl -ane '{chomp;\$range=2;\$code=0;for(\$i=-\$range;\$i<=\$range;\$i++){for(\$j=-\$range;\$j<=\$range;\$j++){\$code++ if(defined(\$h{\"\$F[0] \$F[1] \".(\$F[2]+\$i).\" \".(\$F[3]+\$j)}))}}if(\$code==0){\$h{\"\$F[0] \$F[1] \$F[2] \$F[3]\"}=1}else{print \"\$F[4]\\n\",substr(\$F[4],0,2),int(substr(\$F[4],2))+1,\"\\n\"}}' > redundant_sj.txt\n" if(not(-e "redundant_sj.txt")); 

print FILE "echo 'Chimeric/Redundant jump reads:';wc -l  chimeric_sj.txt redundant_sj.txt;\n";

#remove all chimeric and all redundant reads from sj.cor.fa
print FILE "extractreads.pl <(cat chimeric_sj.txt redundant_sj.txt | perl -e '{while(\$line=<STDIN>){chomp(\$line);\$h{\$line}=1}open(FILE,\$ARGV[0]);while(\$line=<FILE>){chomp(\$line);print \$line,\"\\n\" if(not(defined(\$h{\$line})));}}' <(awk '{print \$1}' work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt)) sj.cor.fa 1 |reverse_complement > sj.cor.clean.fa\n" if(not(-e "sj.cor.clean.fa"));

#we extend the filtered and reverse complemented jump reads if asked -- the reason to do that is that sometimes the jump library reads are shorter than 64 bases and CA cannot use them
if($EXTEND_JUMP_READS==1)
{
print FILE "createSuperReadsForDirectory.perl -minreadsinsuperread 1 -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.fasta -l 31 -s $JF_SIZE -t $NUM_THREADS -M 2 -m 2 -jumplibraryreads -mkudisr 0 work3 sj.cor.clean.fa 1>super3.err 2>&1\n" if(not(-e "work3"));
print FILE "ln -sf work3/superReadSequences.jumpLibrary.fasta sj.cor.ext.fa\n";
}
else
{
print FILE "ln -sf sj.cor.clean.fa sj.cor.ext.fa\n";
} 

print FILE "\n";

#here we create the frg files for CA from the jump libraries: each jump library will contribute one jump frg file and one additional frg file of linking information from "chimers"
print FILE "echo -n 'creating FRG files ';date;\n";

print FILE "rm -rf compute_jump_coverage.txt\n";

for($i=0;$i<scalar(@jump_info_array);$i++)
{
@f=split(/\s+/,$jump_info_array[$i]);
$list_of_frg_files.="$f[0].cor.clean.frg ";
print FILE "echo -n \"$f[1] \" >> compute_jump_coverage.txt\n";
print FILE "grep -A 1 '^>$f[0]' sj.cor.ext.fa > $f[0].tmp\n";
print FILE "error_corrected2frg $f[0] $f[1] $f[2] \$MAX_$f[0] $f[0].tmp |tee $f[0].cor.clean.frg | grep '^{LKG' |wc -l >> compute_jump_coverage.txt\n";
print FILE "rm $f[0].tmp\n";
}
print FILE "JUMP_BASES_COVERED=`awk 'BEGIN{b=0}{b+=\$1*\$2;}END{print b}' compute_jump_coverage.txt`\n";
}
#here we reduce jump library coverage: we know the genome size (from k-unitigs) and JUMP_BASES_COVERED contains total jump library coverage :)
print FILE "perl -e '{\$cov='\$JUMP_BASES_COVERED'/'\$ESTIMATED_GENOME_SIZE'; print \"JUMP insert coverage: \$cov\\n\"; \$optimal_cov=60;if(\$cov>\$optimal_cov){print \"Reducing JUMP insert coverage from \$cov to \$optimal_cov\\n\";\$prob_coeff=\$optimal_cov/\$cov;open(FILE,\"gkp.edits.msg\");while(\$line=<FILE>){chomp(\$line);\@f=split(/\\s+/,\$line);\$deleted{\$f[2]}=1;}close(FILE); open(FILE,\"sj.cor.clean.fa\");while(\$line=<FILE>){next if(not(\$line =~ /^>/));chomp(\$line);if(int(substr(\$line,3))%2==0) {print STDERR substr(\$line,1),\"\\n\" if(rand(1)>\$prob_coeff);}}}}' 2>mates_to_break.txt\n";

for($i=0;$i<scalar(@jump_info_array);$i++)
{
@f=split(/\s+/,$jump_info_array[$i]);
###if one of the reads 
print FILE "mv $f[0].cor.clean.frg $f[0].cor.clean.frg.tmp\n";
print FILE "perl -e '{open(FILE,\"mates_to_break.txt\");while(\$line=<FILE>){chomp(\$line);\$h{\$line}=1;}while(\$line=<STDIN>){if(\$line=~/^\\{LKG/){\$flag=0;\@lines=();while(\$line=<STDIN>){last if(\$line=~/^\\}/);push(\@lines,\$line);chomp(\$line);if(defined(\$h{substr(\$line,4)})){\$flag++}}if(\$flag==0){print \"{LKG\\n\",\@lines,\"}\\n\";}}else{print \$line;}}}' < $f[0].cor.clean.frg.tmp > $f[0].cor.clean.frg\n"; 
print FILE "rm $f[0].cor.clean.frg.tmp\n";
}

###############################################################done with JUMP library################################################################
print FILE "\n\n\n";
########################################################error correct, super reads for PE#################################
print FILE "echo -n 'computing super reads from PE ';date;\n";

#we create super reads from PE
print FILE "createSuperReadsForDirectory.perl -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.fasta -l 31 -s $JF_SIZE -t $NUM_THREADS -M 2 -m 2 -join-mates -join-shooting -mkudisr 0 work1 pe.cor.fa 1>super1.err 2>&1\n" if(not(-e "work1"));
print FILE "\n";

#now we extract those PE mates that did not end up in the same super read -- we call them linking mates, they will be useful for scaffolding
print FILE "extractreads.pl <( awk '{if(int(substr(\$1,3))%2==0){print \$1\" \"\$2\" \"\$1;}else{print \$1\" \"\$2\" \"substr(\$1,1,2)\"\"int(substr(\$1,3))-1}}' work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt |uniq -D -f 2 | uniq -u -f 1 | awk '{print \$1}' )  pe.cor.fa 1 > pe.linking.fa\n" if(not(-e "pe.linking.fa"));
print FILE "\n";

#create frg files for PE data
foreach $v(@pe_info_array)
{
@f=split(/\s+/,$v);
$list_of_frg_files.="$f[0].linking.frg ";
if(not(-e "$f[0].linking.frg"))
{
print FILE "grep -A 1 '^>$f[0]' pe.linking.fa > $f[0].tmp\n";
print FILE "error_corrected2frg $f[0] $f[1] $f[2] \$MAX_$f[0] $f[0].tmp > $f[0].linking.frg\n";
print FILE "rm $f[0].tmp\n";
}
}
print FILE "echo -n 'Linking PE reads ';\ncat ??.linking.frg |grep '^{FRG' |wc -l;\n";
print FILE "\n";
#create frg file for super reads
print FILE "awk 'BEGIN{f=1}{if(f==0){print l\" \"\$1\" \"length(\$1)}else{l=\$1}f=1-f;}' work1/superReadSequences.fasta |sort -grk3,3 -S 20%| awk '{print \$1\"\\n\"\$2}' | create_sr_frg.pl | fasta2frg.pl sr >  superReadSequences_shr.frg\n" if(not(-e "superReadSequences_shr.frg"));
print FILE "\n\n\n";
################################################################done with PE library################################################################
foreach $v(@other_info_array)
{
$list_of_frg_files.="$v ";
}
print FILE "\n\n\n";
###############################################################Celera Assembler################################################################

print FILE "\necho -n 'Celera Assembler ';date;\n";

print FILE "rm -rf CA\n";

#data filtering
if(scalar(@jump_info_array)>0)
{
#we run CA fully for small data sets (under 10M reads) and only up to the unitigs for large genomes (over 10M reads)

print FILE "let TOTAL_READS=`wc -l pe.cor.fa| awk '{print \$1}'`\n";

print FILE "runCA jellyfishHashSize=$JF_SIZE utgErrorRate=0.03 merylMemory=8192 ovlMemory=4GB stopAfter=unitigger ovlMerThreshold=300 bogBreakAtIntersections=0 doOverlapBasedTrimming=0 unitigger=bog bogBadMateDepth=1000000 -p genome -d CA merylThreads=$NUM_THREADS frgCorrThreads=1 frgCorrConcurrency=$NUM_THREADS cnsConcurrency=$NUM_THREADS ovlCorrConcurrency=$NUM_THREADS ovlConcurrency=$NUM_THREADS ovlThreads=1 superReadSequences_shr.frg $list_of_frg_files  1> runCA0.out 2>&1\n\n";

#here we filter libraries for chimerism and redundancy
#we also reduce the insert coverage by jump libraries if necessary: no more than 100x insert coverage by all libraries

print FILE "if [[ ! -e CA/4-unitigger/unitigger.err ]];then\n";
print FILE "echo \"CA failed, check output under CA/ and runCA0.out\"\n";
print FILE "exit\n";
print FILE "fi\n";

print FILE "cd CA/\nmv 4-unitigger 4-unitigger-filter\ncd 4-unitigger-filter\ngrep '^>' ../../sj.cor.clean.fa |awk '{print substr(\$1,2)}' > sj.uid\nfilter_library.sh ../ genome sj.uid 750\n";

#we should not check for redundancy on the extended jump reads -- it will wipe them all out
if($EXTEND_JUMP_READS==0)
{
print FILE "cat genome.redundant.uid |awk '{print \"frg uid \"\$1\" isdeleted 1\"}' > gkp.edits.msg\n";
print FILE "cat genome.chimeric.uid |awk '{print \"frg uid \"\$1\" mateiid 0\"}' >> gkp.edits.msg\n";

}
else
{
print FILE "cat genome.chimeric.uid |awk '{print \"frg uid \"\$1\" mateiid 0\"}'  > gkp.edits.msg\n";
}
print FILE "echo -n \"Deleted reads due to redundancy/chimerism: \"\nwc -l gkp.edits.msg\n";
print FILE "gatekeeper --edit gkp.edits.msg ../genome.gkpStore 1>gatekeeper.err 2>&1\n";
print FILE "cd ../\nrm -rf *.tigStore\nrm -rf *.ovlStore\nrm -rf 0-* 1-* 2-* 3-*\ncd ../\n\n";
print FILE "\n";
print FILE "if [[ \$TOTAL_READS -lt 20000000 ]];then \n";#run the assember to completion, for better filtering
print FILE "runCA $CA_PARAMETERS jellyfishHashSize=$JF_SIZE stopAfter=consensusAfterUnitigger doOverlapBasedTrimming=0 unitigger=bog -p genome -d CA merylThreads=$NUM_THREADS frgCorrThreads=1 frgCorrConcurrency=$NUM_THREADS cnsConcurrency=$NUM_THREADS ovlCorrConcurrency=$NUM_THREADS ovlConcurrency=$NUM_THREADS ovlThreads=1 superReadSequences_shr.frg $list_of_frg_files   1> runCA1.out 2>&1\n";
#now we check if the unitig consensus which is sometimes problematic, failed, and fix the unitigs
print FILE "if [[ -e \"CA/5-consensus/consensus.success\" ]];then\n";
print FILE "echo \"unitig consensus OK\"\n";
print FILE "else\n";
print FILE "echo \"fixing unitig consensus...\"\n";
print FILE "mkdir CA/fix_unitig_consensus\n";
print FILE "cd CA/fix_unitig_consensus\n";
print FILE "cp `which fix_unitigs.sh` .\n";
print FILE "./fix_unitigs.sh genome \n";
print FILE "cd ../../\n";
print FILE "fi\n";
#we now recompute the A-stat for the unitigs based on positions of PE reads in the super-reads
print FILE "recompute_astat_superreads.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt\n";
print FILE "runCA $CA_PARAMETERS jellyfishHashSize=$JF_SIZE unitigger=bog -p genome -d CA cnsConcurrency=$NUM_THREADS computeInsertSize=0 doExtendClearRanges=0 1>runCA2.out 2>&1\n";
print FILE "if [[ -e \"CA/9-terminator/genome.qc\" ]];then\n";
print FILE "echo \"CA success, checking for chimeric mates\"\n";
print FILE "else\n";
print FILE "echo \"CA failed, check output under CA/ and runCA2.out\"\n";
print FILE "exit\n";
print FILE "fi\n";
print FILE "cd CA/4-unitigger-filter\n";
print FILE "echo -n 'Found additional chimeric mates: '\n";
print FILE "grep badOuttie ../9-terminator/genome.posmap.mates | perl -ane 'BEGIN{open(FILE,\"sj.uid\");while(\$line=<FILE>){chomp(\$line); \$h{\$line}=1;}}{print \"frg uid \$F[0] isdeleted 1\\nfrg uid \$F[1] isdeleted 1\\n\" if(defined(\$h{\$F[0]}));}' |tee gkp.edits.final.msg | wc -l\n";
print FILE "gatekeeper --edit gkp.edits.final.msg ../genome.gkpStore 1>gatekeeper.err 2>&1\n"; 
print FILE "cd ../\nrm -rf *.tigStore\nrm -rf *.ovlStore\nrm -rf 0-* 1-* 2-* 3-* 4-unitigger 5-* 6-* 7-* 8-* 9-*\ncd ../\n\n";
print FILE "fi\n";
print FILE "\n\n";
}
#this if statement is here because if OTHER frg is specified, we will have to do OBT, it will slow us down, but it has to be done :(
if(scalar(@other_info_array)>0)
{
print FILE "runCA $CA_PARAMETERS jellyfishHashSize=$JF_SIZE stopAfter=consensusAfterUnitigger doOverlapBasedTrimming=1 unitigger=bog -p genome -d CA merylThreads=$NUM_THREADS frgCorrThreads=1 frgCorrConcurrency=$NUM_THREADS cnsConcurrency=$NUM_THREADS ovlCorrConcurrency=$NUM_THREADS ovlConcurrency=$NUM_THREADS ovlThreads=1 superReadSequences_shr.frg $list_of_frg_files   1> runCA1.out 2>&1\n";
}
else
{
print FILE "runCA $CA_PARAMETERS jellyfishHashSize=$JF_SIZE stopAfter=consensusAfterUnitigger doOverlapBasedTrimming=0 unitigger=bog -p genome -d CA merylThreads=$NUM_THREADS frgCorrThreads=1 frgCorrConcurrency=$NUM_THREADS cnsConcurrency=$NUM_THREADS ovlCorrConcurrency=$NUM_THREADS ovlConcurrency=$NUM_THREADS ovlThreads=1 superReadSequences_shr.frg $list_of_frg_files   1> runCA1.out 2>&1\n";
}
#now we check if the unitig consensus which is sometimes problematic, failed, and fix the unitigs
print FILE "if [[ -e \"CA/5-consensus/consensus.success\" ]];then\n";
print FILE "echo \"unitig consensus OK\"\n";
print FILE "else\n";
print FILE "echo \"fixing unitig consensus...\"\n";
print FILE "mkdir CA/fix_unitig_consensus\n";
print FILE "cd CA/fix_unitig_consensus\n";
print FILE "cp `which fix_unitigs.sh` .\n";
print FILE "./fix_unitigs.sh genome \n";
print FILE "cd ../../\n";
print FILE "fi\n";

#we now recompute the A-stat for the unitigs based on positions of PE reads in the super-reads
print FILE "recompute_astat_superreads.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt\n";

#and we continue into the scaffolder... we do ECR only if OTHER data is specified
if(scalar(@other_info_array)>0)
{
print FILE "runCA $CA_PARAMETERS unitigger=bog -p genome -d CA cnsConcurrency=$NUM_THREADS computeInsertSize=0 doExtendClearRanges=1 1>runCA2.out 2>&1\n";
}
else
{
print FILE "runCA $CA_PARAMETERS unitigger=bog -p genome -d CA cnsConcurrency=$NUM_THREADS computeInsertSize=0 doExtendClearRanges=0 1>runCA2.out 2>&1\n";
}

print FILE "if [[ -e \"CA/9-terminator/genome.qc\" ]];then\n";
print FILE "echo \"CA success\"\n";
print FILE "else\n";
print FILE "echo \"CA failed, check output under CA/ and runCA2.out\"\n";
print FILE "exit\n";
print FILE "fi\n";

#here we close gaps in scaffolds:  we use create_k_unitigs allowing to continue on count 1 sequence and then generate fake reads from the 
#end sequences of contigs that are next to each other in scaffolds, and then use super reads software to close the gaps

#Done !!!! Hoorayyyy!!! :)
print FILE "echo -n 'All done ';date;\n";


close(FILE);
system("chmod 0755 assemble.sh");
print "done.\nexecute assemble.sh to run the assembly\n";


