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
my $CA_PARAMETERS=" utgErrorRate=0.03 merylMemory=8192 ovlMemory=4GB ";
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

print "\nPARAMETERS:\nWINDOW = $WINDOW\nMAX_ERR_PER_WINDOW = $MAX_ERR_PER_WINDOW\nTRIM_PARAM = $TRIM_PARAM\nNUM_THREADS = $NUM_THREADS\nJF_SIZE = $JF_SIZE\nEXTEND_JUMP_READS = $EXTEND_JUMP_READS\n\n";
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
print FILE "echo Q2_CHAR = \$Q2_CHAR\n";
print FILE "\ncat $list_pe_files | filter_illumina_data_quals '\$Q2_CHAR'  > pe_trim.fa\n" if(not(-e "pe_trim.fa"));
print FILE "jellyfish count -t $NUM_THREADS -C -r -o pe_trim -s $JF_SIZE -m 24 pe_trim.fa\n"  if(not(-e "pe_trim_0"));
print FILE "jellyfish count -t $NUM_THREADS -C -r -o pe_all -s $JF_SIZE -m 24 $list_pe_files\n"  if(not(-e "pe_all_0"));

print FILE "if [[ -e pe_trim_1 || -e pe_all_1 ]];then\n";
print FILE "echo \"Increase JF_SIZE\"\n";
print FILE "exit\n";
print FILE "fi\n";
########################################################done Jellyfish#################################

########################################################error correct, super reads for PE#################################
if(not(-e "pe.cor.fa"))
{
print FILE "echo -n 'error correct PE ';date;\n";
print FILE "\nerror_correct_reads -d pe_trim_0 -d pe_all_0 -C -m 1 -s 1 -g 1 -t $NUM_THREADS -w $WINDOW -e $MAX_ERR_PER_WINDOW $list_pe_files\n";
print FILE "cat error_corrected_*.fa  | homo_trim $TRIM_PARAM > pe.cor.fa\n";
print FILE "rm error_corrected_*\n";
}

print FILE "echo -n 'computing super reads from PE ';date;\n";

print FILE "\n";
print FILE "createSuperReadsForDirectory.perl -l 31 -s $JF_SIZE -t $NUM_THREADS -M 2 -m 2 -join-mates -mkudisr 0 work1 pe.cor.fa 1> super.out 2>super.err\n" if(not(-e "work1"));

print FILE "\n";
print FILE "extractreads.pl <( awk '{if(int(substr(\$1,3))%2==0){print \$1\" \"\$2\" \"\$1;}else{print \$1\" \"\$2\" \"substr(\$1,1,2)\"\"int(substr(\$1,3))-1}}' work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt |sort -k3,3 -S 10% |uniq -D -f 2 | uniq -u -f 1 | awk '{print \$1}' )  pe.cor.fa 1 > pe.linking.fa\n" if(not(-e "pe.linking.fa"));
print FILE "\n";
$list_of_frg_files="";
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
print FILE "create_sr_frg.pl < work1/superReadSequences.fasta | fasta2frg.pl sr >  superReadSequences_shr.frg\n" if(not(-e "superReadSequences_shr.frg"));
print FILE "\n";

###############################################################done with PE library################################################################

########################################################error correct, super reads for JUMP#################################
if(scalar(@jump_info_array)>0)
{
if(not(-e "sj.cor.fa"))
{
print FILE "echo -n 'error correct JUMP ';date;\n";
print FILE "\nerror_correct_reads -d pe_trim_0 -d pe_all_0 -C -m 1 -s 1 -g 2 -t $NUM_THREADS -w $WINDOW -e $MAX_ERR_PER_WINDOW $list_jump_files\n";
print FILE "cat error_corrected_*.fa  | homo_trim $TRIM_PARAM > sj.cor.fa\n";
print FILE "rm error_corrected_*\n";
}

print FILE "echo -n 'filtering JUMP ';date;\n";

print FILE "createSuperReadsForDirectory.perl -minreadsinsuperread 1 -kunitigsfile work1/guillaumeKUnitigsAtLeast32bases_all.fasta -l 31 -s $JF_SIZE -t $NUM_THREADS -M 2 -m 2 -join-mates -mkudisr 0 work2 sj.cor.fa 1> super1.out 2>super1.err\n" if(not(-e "work2"));;

print FILE "\n";

print FILE "awk '{if(int(substr(\$1,3))%2==0){print \$3\" \"\$2\" \"\$1;}else{print \$3\" \"\$2\" \"substr(\$1,1,2)\"\"int(substr(\$1,3))-1}}' work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt |sort -k3,3 -S 10% |uniq -D -f 1| awk 'BEGIN{insert=\"\";}{if(\$3!=insert){start=\$1;insert=\$3}else{if(start>\$1){print insert\" \"start-\$1}else{print insert\" \"\$1-start}}}' | awk '{if(\$2<1000){print \$1\"\\n\"substr(\$1,0,2)\"\"int(substr(\$1,3))+1}}' > chimeric_sj.txt\n" if(not(-e "chimeric_sj.txt"));

print FILE "echo -n 'Chimeric jump reads ';wc -l  chimeric_sj.txt;\n";

print FILE "extractreads.pl <(cat chimeric_sj.txt <(awk '{print \$1}' work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt) |sort -S 10% |uniq -u) sj.cor.fa 1 |reverse_complement > sj.cor.clean.fa\n" if(not(-e "sj.cor.clean.fa"));

if($EXTEND_JUMP_READS==1)
{
print FILE "createSuperReadsForDirectory.perl -minreadsinsuperread 1 -kunitigsfile work1/guillaumeKUnitigsAtLeast32bases_all.fasta -l 31 -s $JF_SIZE -t $NUM_THREADS -M 2 -m 2 -join-mates -jumplibraryreads -mkudisr 0 work2 sj.cor.clean.fa 1> super2.out 2>super2.err\n" if(not(-e "work3"));;
print FILE "ln -sf work3/superReadSequences.jumpLibrary.fasta sj.cor.ext.fa\n";
}
else
{
print FILE "ln -sf sj.cor.clean.fa sj.cor.ext.fa\n";
} 

print FILE "\n";

print FILE "echo -n 'creating FRG files ';date;\n";

foreach $v(@jump_info_array)
{
@f=split(/\s+/,$v);
$list_of_frg_files.="$f[0].cor.clean.frg ";
if(not(-e "$f[0].cor.clean.frg"))
{
print FILE "grep -A 1 '^>$f[0]' sj.cor.ext.fa > $f[0].tmp\n";
print FILE "error_corrected2frg $f[0] $f[1] $f[2] \$MAX_$f[0] $f[0].tmp > $f[0].cor.clean.frg\n";
print FILE "rm $f[0].tmp\n";
}
}
}
###############################################################done with JUMP library################################################################
foreach $v(@other_info_array)
{
$list_of_frg_files.="$v ";
}

###############################################################Celera Assembler################################################################
print FILE "echo -n 'Celera Assembler ';date;\n";

print FILE "rm -rf CA\n";

if(scalar(@jump_info_array)>0)
{
print FILE "runCA $CA_PARAMETERS stopAfter=unitigger ovlMerThreshold=300 bogBreakAtIntersections=0 doOverlapBasedTrimming=0 unitigger=bog bogBadMateDepth=1000000 -p genome -d CA merylThreads=$NUM_THREADS frgCorrThreads=1 frgCorrConcurrency=$NUM_THREADS cnsConcurrency=$NUM_THREADS ovlCorrConcurrency=$NUM_THREADS ovlConcurrency=$NUM_THREADS ovlThreads=1 superReadSequences_shr.frg $list_of_frg_files  1> runCA0.out 2>&1\n\n";

print FILE "cd CA/\nmv 4-unitigger 4-unitigger-filter\ncd 4-unitigger-filter\ngrep '^>' ../../sj.cor.fa |awk '{print substr(\$1,2)}' > sj.uid\nfilter_library.sh ../ genome sj.uid 1000 \n";
print FILE "cat genome.redundant.uid genome.chimeric.uid |awk '{print \"frg uid \"\$1\" isdeleted 1\"}' > gkp.edits.msg\n";
print FILE "gatekeeper --edit gkp.edits.msg ../genome.gkpStore 1>gatekeeper.err 2>&1\n";
print FILE "cd ../\nrm -rf *.tigStore\nrm -rf *.ovlStore\nrm -rf 0-overlaptrim-overlap\nrm -rf 0-overlapper\nrm -rf 1-* 2-* 3-*\ncd ../\n\n";
}

print FILE "runCA $CA_PARAMETERS stopAfter=consensusAfterUnitigger doOverlapBasedTrimming=0 unitigger=bog -p genome -d CA merylThreads=$NUM_THREADS frgCorrThreads=1 frgCorrConcurrency=$NUM_THREADS cnsConcurrency=$NUM_THREADS ovlCorrConcurrency=$NUM_THREADS ovlConcurrency=$NUM_THREADS ovlThreads=1 superReadSequences_shr.frg $list_of_frg_files   1> runCA1.out 2>&1\n";
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
print FILE "recompute_astat_superreads.sh genome CA 100 work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt\n";
print FILE "runCA unitigger=bog -p genome -d CA cnsConcurrency=$NUM_THREADS computeInsertSize=0 doExtendClearRanges=0 cgwErrorRate=0.15 1>runCA2.out 2>&1\n";

print FILE "echo -n 'All done ';date;\n";

close(FILE);
system("chmod 0755 assemble.sh");
print "done.\nexecute assemble.sh to run the assembly\n";

