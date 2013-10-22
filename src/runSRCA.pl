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
#KMER_COUNT_THRESHOLD=1
#NUM_THREADS=24
#JF_SIZE=2000000000
#EXTEND_JUMP_READS=0
#END
#

# use strict;
# use warnings;
use Env qw(@PATH @LD_LIBRARY_PATH);
use File::Spec;

#parsing config file
my $KMER="auto";
my $JELLYFISH_PATH="";
my $SR_PATH="";
my $CA_PATH="";
my $CA_PARAMETERS="";
my $WINDOW=10;
my $MAX_ERR_PER_WINDOW=3;
my $KMER_COUNT_THRESHOLD=1;
my $KMER_RELIABLE_THRESHOLD=3;
my $TRIM_PARAM=2;
my $NUM_THREADS=2;
my $EXTEND_JUMP_READS=1;
my $JF_SIZE=100000000;
my $in_paths=0;
my $in_data=0;
my $in_parameters=0;
my @pe_info_array=();
my @jump_info_array=();
my @other_info_array=();
my %used_library_ids;
my $LIMIT_JUMP_COVERAGE=500;
my $USE_LINKING_MATES=1;
my $homo_trim_string=" | homo_trim  $TRIM_PARAM ";
my $CLOSE_GAPS=1;

if(not(-e $ARGV[0])||$ARGV[0] eq ""|| $ARGV[0] eq "-h"){
    print "USAGE: runSRCA.pl <config_file>\n";
    exit;
}

my $config_file = $ARGV[0];
open(FILE, $config_file);
while($line=<FILE>){
    chomp($line);
    next if($line =~ /^\#/);
    if($line =~ /^PATHS/){
	if($in_parameters==1){
	    die("error in config file: mixed PARAMETERS and PATHS");
	}
	if($in_data==1){
	    die("error in config file: mixed PATHS and DATA");
	}
	if($in_paths==1){
	    die("duplicate PATHS header"); 
	}
	$in_paths=1;
	next;
    }
    elsif($line =~ /^DATA/){
	if($in_parameters==1){
	    die("error in config file: mixed PARAMETERS and DATA");
	}
	if($in_paths==1){
	    die("error in config file: mixed PATHS and DATA");
	}
	if($in_data==1){
	    die("duplicate PATHS header");
	}
	$in_data=1;
	next;
    }
    elsif($line =~ /^PARAMETERS/){
	if($in_data==1){
	    die("error in config file: mixed PARAMETERS and DATA");
	}
	if($in_paths==1){
	    die("error in config file: mixed PARAMETERS and PATHS");
	}
	if($in_parameters==1){
	    die("duplicate PARAMETERS header");
	}
	$in_parameters=1;
	next;
    }

    if($in_parameters==1){
	if($line =~ /^END/){
	    $in_parameters=0;
	    next;
	}
	elsif($line =~ /^EXTEND_JUMP_READS/){
	    @f=split(/=/,$line);
	    $f[1]=~s/^\s+//;
	    $f[1]=~s/\s+$//;
	    $EXTEND_JUMP_READS=int($f[1]);
	    die("bad value for EXTEND_JUMP_READS, enter 0 or 1") if($EXTEND_JUMP_READS!=1 && $EXTEND_JUMP_READS!=0);
	    next;
	}
        elsif($line =~ /^DO_HOMOPOLYMER_TRIM/){
            @f=split(/=/,$line);
            $f[1]=~s/^\s+//;
            $f[1]=~s/\s+$//;
            $DO_HOMOPOLYMER_TRIM=int($f[1]);
	    $homo_trim_string="" if($DO_HOMOPOLYMER_TRIM==0);    
            die("bad value for DO_HOMOPOLYMER_TRIM, enter 0 or 1") if($DO_HOMOPOLYMER_TRIM!=1 && $DO_HOMOPOLYMER_TRIM!=0);
            next;
        }
	elsif($line =~ /^CA_PARAMETERS/){
	    @f=split(/=/,$line);
	    $CA_PARAMETERS=$f[1];
	    for($i=2;$i<scalar(@f);$i++){
		$CA_PARAMETERS.="=$f[$i]";
	    }
	    die("bad value for CA_PARAMETERS") if($CA_PARAMETERS eq "");
	    next;
	}
        elsif($line =~ /^LIMIT_JUMP_COVERAGE/){
            @f=split(/=/,$line);
            $f[1]=~s/^\s+//;
            $f[1]=~s/\s+$//;
            $LIMIT_JUMP_COVERAGE=int($f[1]);
            die("bad value for LIMIT_JUMP_COVERAGE, enter int > 1") if($LIMIT_JUMP_COVERAGE<1);
            next;
        }
        elsif($line =~ /^GRAPH_KMER_SIZE/){
            @f=split(/=/,$line);
            $f[1]=~s/^\s+//;
            $f[1]=~s/\s+$//;
            $KMER=$f[1];
            die("bad value for GRAPH_KMER_SIZE, enter auto or number > 15 and < 151") if(not($KMER eq "auto") && ($KMER<15 || $KMER>151));
            next;
        }
        elsif($line =~ /^USE_LINKING_MATES/){
            @f=split(/=/,$line);
            $f[1]=~s/^\s+//;
            $f[1]=~s/\s+$//;
            $USE_LINKING_MATES=int($f[1]);
            die("bad value for USE_LINKING_MATES, enter 0 or 1") if($USE_LINKING_MATES!=1 && $USE_LINKING_MATES!=0);
            next;
        }
	elsif($line =~ /^KMER_COUNT_THRESHOLD/){
	    @f=split(/=/,$line);
	    $f[1]=~s/^\s+//;
	    $f[1]=~s/\s+$//;
	    $KMER_COUNT_THRESHOLD=int($f[1]);
	    $KMER_RELIABLE_THRESHOLD=3*$KMER_COUNT_THRESHOLD;
	    die("bad value for KMER_COUNT_THRESHOLD") if($KMER_COUNT_THRESHOLD<1);
	    next;
	}
	elsif($line =~ /^NUM_THREADS/){
	    @f=split(/=/,$line);
	    $f[1]=~s/^\s+//;
	    $f[1]=~s/\s+$//;
	    $NUM_THREADS=int($f[1]);
	    die("bad value for NUM_THREADS") if($NUM_THREADS<1);
	    next;
	}
	elsif($line =~ /^JF_SIZE/){
	    @f=split(/=/,$line);
	    $f[1]=~s/^\s+//;
	    $f[1]=~s/\s+$//;
	    $JF_SIZE=int($f[1]);
	    die("bad value for JF_SIZE, enter int > 100000") if($JF_SIZE<100000);
	    next;
	}
	else{
	    next;
	}
    }


    if($in_paths==1){
	if($line =~ /^END/){
	    $in_paths=0;
	    next;
	}
	elsif($line =~ /^JELLYFISH_PATH/){
	    @f=split(/=/,$line);
	    $JELLYFISH_PATH=$f[1];
	    $JELLYFISH_PATH=~s/^\s+//;
	    $JELLYFISH_PATH=~s/\s+$//;
	    $JELLYFISH_PATH.="/" if(not($JELLYFISH_PATH eq ""));
	    next;
	}
	elsif($line =~ /^CA_PATH/){
	    @f=split(/=/,$line);
	    $CA_PATH=$f[1];
	    $CA_PATH=~s/^\s+//;
	    $CA_PATH=~s/\s+$//;
	    $CA_PATH.="/" if(not($CA_PATH eq ""));
	    next;
	}
	elsif($line =~ /^SR_PATH/){
	    @f=split(/=/,$line);
	    $SR_PATH=$f[1];
	    $SR_PATH=~s/^\s+//;
	    $SR_PATH=~s/\s+$//;
	    $SR_PATH.="/" if(not($SR_PATH eq ""));
	    next;
	}
	else{
	    next;
	}
    }
    elsif($in_data==1){
	if($line =~ /^END/){
	    if(scalar(@pe_info_array)==0 && scalar(@jump_info_array)==0 && scalar(@other_info_array)==0){
		die("incomplete DATA section");
	    }
	    $in_data=0;
	    next;
	}
	elsif($line =~ /^PE/){
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
	    if(length($f1[4])>0){
		die("missing reverse file for PE library $f1[0]") if(not(-e $f1[4]));
		$pe_info_line.="$f1[4] ";
	    }
	    else{
		$pe_info_line.="$f1[3] ";
	    }
	    push(@pe_info_array,$pe_info_line);
	}
	elsif($line =~ /^JUMP/){
	    $jump_info_line="";
	    @f=split(/=/,$line);
	    $f[1]=~s/^\s+//;
	    @f1=split(/\s+/,$f[1]);
	    die("improper id for JUMP library $f1[0]") if(not(length($f1[0])==2));
	    die("duplicate id for JUMP library $f1[0]") if(defined($used_library_ids{$f1[0]}));
	    $jump_info_line.="$f1[0] ";
	    $used_library_ids{$f1[0]}=1;
	    die("improper mean for JUMP library $f1[0]") if(int($f1[1])==0);
	    $jump_info_line.="$f1[1] ";
	    die("improper stdev for JUMP library $f1[0]") if(not(int($f1[2])>0));
	    $jump_info_line.="$f1[2] ";
	    die("missing forward file for JUMP library $f1[0]") if(not(-e $f1[3]));
	    $jump_info_line.="$f1[3] ";
	    die("missing reverse file for JUMP library $f1[0]") if(not(-e $f1[4]));
	    $jump_info_line.="$f1[4] ";
	    push(@jump_info_array,$jump_info_line);
	}
	elsif($line=~ /OTHER/){
	    @f=split(/=/,$line);
	    $f[1]=~s/^\s+//;
	    $f[1]=~s/\s+$//;
	    if($f[1] =~/\.frg$/ && (-e $f[1])){
		push(@other_info_array,$f[1]);
	    }
	    else{
		die("missing or incorrect frg file for OTHER library");
	    }
	}
    }
}
die("no read data files specified") if(scalar(@pe_info_array)+scalar(@jump_info_array)+scalar(@other_info_array)==0);

print "PATHS:\nCA_PATH = $CA_PATH\nJELLYFISH_PATH = $JELLYFISH_PATH\nSR_PATH = $SR_PATH\n\n";
print "Verifying PATHS...\n";
unshift(@PATH, $JELLYFISH_PATH, $SR_PATH, $CA_PATH);
unshift(@LD_LIBRARY_PATH, "${JELLYFISH_PATH}/../lib");

print("@PATH\n");
sub check_exec {
  for my $e (@_) {
    system("$e --help >/dev/null 2>&1") == 0 or die "$e not found or failed to run";
    print "$e OK\n";
  }
}
check_exec "jellyfish", "runCA", "createSuperReadsForDirectory.perl";

print "\nPARAMETERS:\nGRAPH_KMER_SIZE = $KMER\nUSE_LINKING_MATES = $USE_LINKING_MATES\nKMER_COUNT_THRESHOLD = $KMER_COUNT_THRESHOLD\nLIMIT_JUMP_COVERAGE = $LIMIT_JUMP_COVERAGE\nNUM_THREADS = $NUM_THREADS\nJF_SIZE = $JF_SIZE\nEXTEND_JUMP_READS = $EXTEND_JUMP_READS\nCA_PARAMETERS = $CA_PARAMETERS\nDO_HOMOPOLYMER_TRIM = $DO_HOMOPOLYMER_TRIM\n\n";
print "DATA:\nPE:\n";
foreach $v(@pe_info_array){
    print "$v\n";
}
print "JUMP:\n";
foreach $v(@jump_info_array){
    print "$v\n";
}
print "OTHER:\n";
foreach $v(@other_info_array){
    print "$v\n";
}

die ("no PE data specified") if(scalar(@pe_info_array)==0);

print("\ncreating script file for the actions...");

my $config_abs_path = File::Spec->rel2abs($config_file);
my $cmd_abs_path = File::Spec->rel2abs(__FILE__);

open(FILE,">assemble.sh");
print FILE <<"EOS";
#!/bin/bash

#assemble.sh generated by runSRCA.pl
CONFIG_PATH="$config_abs_path"
CMD_PATH="$cmd_abs_path"

# Parse command line switches
usage() {
  echo "Usage: \$0 [-r] [-c]"
  exit 1
}
while getopts ":rc" o; do
  case "\${o}" in
    c)
    echo "configuration file is '\$CONFIG_PATH'"
    ;;
    r)
    echo "Rerunning configuration"
    exec perl "\$CMD_PATH" "\$CONFIG_PATH"
    ;;
    *)
    usage
    ;;
  esac
done

set +e
# Set some paths and prime system to save environment variables
save () {
  (echo -n "$1=\""; eval "echo -n \"\$$1\""; echo "\"") >> environment.sh
}

rm -f environment.sh; touch environment.sh
export LD_LIBRARY_PATH="$JELLYFISH_PATH/../lib:\$LD_LIBRARY_PATH"
save LD_LIBRARY_PATH
export PATH="$CA_PATH:$SR_PATH:$JELLYFISH_PATH:\$PATH"
save PATH
EOS

$i=0;
$list_pe_files="";
$list_jump_files="";
$list_of_frg_files="";
$rerun_pe=0;
$rerun_sj=0;

my $tmplist;
foreach $v(@other_info_array){
    $list_of_frg_files.="$v ";
    $tmplist.="$v ";
}

###renaming reads###

print FILE "echo -n 'processing PE library reads ';date;\n";
print FILE "rm -rf meanAndStdevByPrefix.pe.txt\n";
foreach $v(@pe_info_array){
    @f=split(/\s+/,$v);
    $list_pe_files.="$f[0].renamed.fastq ";
    print FILE "echo '$f[0] $f[1] $f[2]' >> meanAndStdevByPrefix.pe.txt\n";
    next if(-e "$f[0].renamed.fastq");
    $rerun_pe=1;
    $rerun_sj=1;
    print FILE "rename_filter_fastq.pl $f[0] $f[3] $f[4] > $f[0].renamed.fastq &\nPID$i=\$!\n";   
    $i++;
}
print FILE "wait ";
for(my $j=0;$j<$i;$j++){
    print FILE "\$PID$j ";
}
print FILE "\n";

print FILE "rm -rf meanAndStdevByPrefix.sj.txt\n";
if(scalar(@jump_info_array)>0){
    print FILE "echo -n 'processing JUMP library reads ';date;\n";
    $i=0;
    foreach $v(@jump_info_array){
	@f=split(/\s+/,$v);
	print FILE "echo '$f[0] 500 100' >> meanAndStdevByPrefix.sj.txt\n";
	$list_jump_files.="$f[0].renamed.fastq ";
	next if(-e "$f[0].renamed.fastq");
        $rerun_sj=1;
	die("duplicate jump library $f[0] files") if($f[3] eq $f[4]);
	print FILE "rename_filter_fastq.pl $f[0] $f[3] $f[4] > $f[0].renamed.fastq &\nPID$i=\$!\n";
	$i++;
    }
    print FILE "wait ";
    for(my $j=0;$j<$i;$j++){
	print FILE "\$PID$j ";
    }
    print FILE "\n";
}
###done renaming reads###
print FILE "\n";
###compute minimum and average PE read length and gc content, and kmer size###
print FILE "head -q -n 2  $list_pe_files | grep --text -v '^\@' > pe_data.tmp\n";
print FILE "PE_AVG_READ_LENGTH=`awk '{n+=length(\$1);m++;}END{print int(n/m)}' pe_data.tmp`\n";
print FILE "save PE_AVG_READ_LENGTH\n";
print FILE "echo \"Average PE read length \$PE_AVG_READ_LENGTH\"\n";
if(uc($KMER) eq "AUTO"){
#here we have to estimate gc content and recompute kmer length for the graph
     print FILE "KMER=`for f in $list_pe_files;do head -n 80000 \$f |tail -n 40000;done | perl -e 'while(\$line=<STDIN>){\$line=<STDIN>;chomp(\$line);push(\@lines,\$line);\$line=<STDIN>;\$line=<STDIN>}\$min_len=100000;\$base_count=0;foreach \$l(\@lines){\$base_count+=length(\$l);if(length(\$l)<\$min_len){\$min_len=length(\$l)} \@f=split(\"\",\$l);foreach \$base(\@f){if(uc(\$base) eq \"G\" || uc(\$base) eq \"C\"){\$gc_count++}}} \$gc_ratio=\$gc_count/\$base_count;\$kmer=0;if(\$gc_ratio<0.5){\$kmer=int(\$min_len*.7);}elsif(\$gc_ratio>=0.5 && \$gc_ratio<0.6){\$kmer=int(\$min_len*.5);}else{\$kmer=int(\$min_len*.33);} \$kmer=31 if(\$kmer<31); print \$kmer'`\n";
     print FILE "save KMER\n";
     print FILE "echo \"choosing kmer size of \$KMER for the graph\"\n";
}else{
     print FILE "KMER=$KMER\n";
}

if(scalar(@jump_info_array)>0){    
    print FILE "KMER_J=`for f in $list_jump_files;do head -n 80000 \$f |tail -n 40000;done | perl -e 'while(\$line=<STDIN>){\$line=<STDIN>;chomp(\$line);push(\@lines,\$line);\$line=<STDIN>;\$line=<STDIN>}\$min_len=100000;\$base_count=0;foreach \$l(\@lines){\$base_count+=length(\$l);if(length(\$l)<\$min_len){\$min_len=length(\$l)} \@f=split(\"\",\$l);foreach \$base(\@f){if(uc(\$base) eq \"G\" || uc(\$base) eq \"C\"){\$gc_count++}}} \$gc_ratio=\$gc_count/\$base_count;\$kmer=0;if(\$gc_ratio<0.5){\$kmer=int(\$min_len*.7);}elsif(\$gc_ratio>=0.5 && \$gc_ratio<0.6){\$kmer=int(\$min_len*.5);}else{\$kmer=int(\$min_len*.33);} \$kmer=31 if(\$kmer<31);  print \$kmer'`\n";
    print FILE "save KMER_J\n";
}else{
   print FILE "KMER_J=\$KMER\n";
}
###done###
###Jellyfish###
print FILE "echo -n 'running Jellyfish ';date;\n";
print FILE "MIN_Q_CHAR=`cat $list_pe_files |head -n 50000 | awk 'BEGIN{flag=0}{if(\$0 ~ /^\+/){flag=1}else if(flag==1){print \$0;flag=0}}'  | perl -ne 'BEGIN{\$q0_char=\"\@\";}{chomp;\@f=split \"\";foreach \$v(\@f){if(ord(\$v)<ord(\$q0_char)){\$q0_char=\$v;}}}END{\$ans=ord(\$q0_char);if(\$ans<64){print \"33\\n\"}else{print \"64\\n\"}}'`\n";
print FILE "save MIN_Q_CHAR\n";
print FILE "echo MIN_Q_CHAR: \$MIN_Q_CHAR\n";

#for our error correction we run jellyfish twice: once on all pe bases and once on pe bases with quality >2 
print FILE "JF_SIZE=`ls -l *.fastq | awk '{n+=\$5}END{s=int(n/46); if(s>$JF_SIZE)print s;else print \"$JF_SIZE\";}'`\n";
print FILE "save JF_SIZE\n";
print FILE "perl -e '{if(int('\$JF_SIZE')>$JF_SIZE){print \"WARNING: JF_SIZE set too low, increasing JF_SIZE to at least '\$JF_SIZE', this automatic increase may be not enough!\\n\"}}'\n"; 
if(not(-e "combined_0") || $rerun_pe==1){
    print FILE "jellyfish count -t $NUM_THREADS -p 126 -C -r -o pe_trim -s \$JF_SIZE -m 24 --quality-start \$MIN_Q_CHAR --min-quality 5 $list_pe_files\n";
    print FILE "jellyfish count -t $NUM_THREADS -p 126 -C -r -o pe_all -s \$JF_SIZE -m 24 $list_pe_files\n";

    print FILE "CUTOFF=`jellyfish histo -t 48 pe_trim_0 | 
        perl -ane '{
        if(\$F[0]>1){
                \$n+=\$F[0]*\$F[1];
                \$m+=\$F[1];
                }
        }END{
        \$er=0.01; 
        \$c=\$n/\$m; 
        for(\$x=2;\$x<1000;\$x++){
                if(\$er/sqrt(2*3.1415927*\$x)*((\$er*\$c/3/\$x)**\$x)*exp(-\$er*\$c/3+\$x)<0.000001){
                        print (\$x+1);
                        last;
                }
        }
        }'`\n";

    print FILE "save CUTOFF\n";
    print FILE "echo \"Error correction Poisson cutoff = \$CUTOFF\"\n";
    print FILE "echo \$CUTOFF > cutoff.txt\n";
#check if the JF_SIZE was big enough:  we want to end up with a single raw database for pe_all and pe_trim
    print FILE "if [[ -e pe_trim_1 || -e pe_all_1 ]];then\n";
    print FILE "echo \"Increase JF_SIZE in config file, the recommendation is to set this to genome_size*coverage/2\"\n";
    print FILE "rm -f pe_trim_? pe_all_?\n";
    print FILE "exit\n";
    print FILE "fi\n";

    print FILE "quorum_combine_jf_dbs -m 1 pe_trim_0 pe_all_0 -o combined\n";
    print FILE "rm -f pe_trim_? pe_all_?\n";
    $rerun_pe=1;
    $rerun_sj=1;
}

###done Jellyfish###
print FILE "\n";
###error correct PE###

if(not(-e "pe.cor.fa")||$rerun_pe==1){
    print FILE "echo -n 'error correct PE ';date;\n";
    print FILE "cat combined_0 > /dev/null\n";
    print FILE "\nquorum_error_correct_reads   -p `cat cutoff.txt` --contaminant=$SR_PATH/../share/adapter_0 -d combined_0 -c 2 -C -m $KMER_COUNT_THRESHOLD -s 1 -g 1 -a $KMER_RELIABLE_THRESHOLD -t $NUM_THREADS -w $WINDOW -e $MAX_ERR_PER_WINDOW $list_pe_files 2>error_correct.log $homo_trim_string | add_missing_mates.pl > pe.cor.fa\n";
    $rerun_pe=1;
}

###done error correct PE###
print FILE "\n";
###error correct JUMP###

if(scalar(@jump_info_array)>0){
    if(not(-e "sj.cor.fa")||$rerun_sj==1){
	print FILE "echo -n 'error correct JUMP ';date;\n";
        print FILE "cat combined_0 > /dev/null\n";
	print FILE "\nquorum_error_correct_reads -p `cat cutoff.txt` --contaminant=$SR_PATH/../share/adapter_0 -d combined_0 -c 2 -C -m $KMER_COUNT_THRESHOLD -s 1 -g 2 -a $KMER_RELIABLE_THRESHOLD -t $NUM_THREADS -w $WINDOW -e $MAX_ERR_PER_WINDOW $list_jump_files 2>error_correct.log $homo_trim_string | add_missing_mates.pl > sj.cor.fa\n";
        $rerun_sj=1;
    }
}

###done error correct JUMP###
print FILE "\n";
###estimate genome size###

if(scalar(@jump_info_array)>0){
    $k_u_arg="pe.cor.fa sj.cor.fa";
}else{
    $k_u_arg="pe.cor.fa";
}

if(not(-e "k_u_hash_0")||$rerun_pe==1||$rerun_sj==1){
    print FILE "jellyfish-2.0 count -m 31 -t $NUM_THREADS -C -s \$JF_SIZE -o k_u_hash_0 $k_u_arg\n";
}


print FILE "ESTIMATED_GENOME_SIZE=`jellyfish-2.0 histo -t $NUM_THREADS -h 1 k_u_hash_0 | tail -n 1 |awk '{print \$2}'`\n";
print FILE "save ESTIMATED_GENOME_SIZE\n";
print FILE "echo \"Estimated genome size: \$ESTIMATED_GENOME_SIZE\"\n";

####done estimate genome size###
print FILE "\n";
###build k-unitigs###

if(not(-e "guillaumeKUnitigsAtLeast32bases_all.fasta")||$rerun_pe==1||$rerun_sj==1){
    print FILE "create_k_unitigs_large_k -c \$((\$KMER-1)) -t $NUM_THREADS -m \$KMER -n \$ESTIMATED_GENOME_SIZE -l \$KMER -f 0.000001 $k_u_arg  | grep --text -v '^>' | perl -ane '{\$seq=\$F[0]; \$F[0]=~tr/ACTGacgt/TGACtgac/;\$revseq=reverse(\$F[0]); \$h{(\$seq ge \$revseq)?\$seq:\$revseq}=1;}END{\$n=0;foreach \$k(keys \%h){print \">\",\$n++,\" length:\",length(\$k),\"\\n\$k\\n\"}}' > guillaumeKUnitigsAtLeast32bases_all.fasta\n";
    print FILE "if [[ \$KMER -eq \$KMER_J ]];then\n";
    print FILE "ln -s guillaumeKUnitigsAtLeast32bases_all.fasta guillaumeKUnitigsAtLeast32bases_all.jump.fasta\n";
    print FILE "else\n";
    print FILE "create_k_unitigs_large_k -c \$((\$KMER_J-1)) -t $NUM_THREADS -m \$KMER_J -n \$ESTIMATED_GENOME_SIZE -l \$KMER_J -f 0.000001 $k_u_arg  | grep --text -v '^>' | perl -ane '{\$seq=\$F[0]; \$F[0]=~tr/ACTGacgt/TGACtgac/;\$revseq=reverse(\$F[0]); \$h{(\$seq ge \$revseq)?\$seq:\$revseq}=1;}END{\$n=0;foreach \$k(keys \%h){print \">\",\$n++,\" length:\",length(\$k),\"\\n\$k\\n\"}}' > guillaumeKUnitigsAtLeast32bases_all.jump.fasta\n";
    print FILE "fi\n";
    $rerun_pe=1;
    $rerun_sj=1;
}

###done build k-unitigs###
print FILE "\n";
###super reads and filtering for jump###

if( not(-d "CA") || $rerun_pe || $rerun_sj ){
    if(scalar(@jump_info_array)>0){
	print FILE "echo -n 'filtering JUMP ';date;\n";

#creating super reads for filtering
	if($rerun_pe==1||$rerun_sj==1||not(-e "work2")){
	    print FILE "rm -rf work2\n";
	    $rerun_sj=1;
	}
	print FILE "createSuperReadsForDirectory.perl  -maxnodes 100000 -minreadsinsuperread 1 -l \$KMER_J -join-aggressive 1 -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.sj.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.jump.fasta -t $NUM_THREADS -mikedebug work2 sj.cor.fa 1> super2.err 2>&1\n";

#check if the super reads pipeline finished successfully
	print FILE "if [[ ! -e work2/superReads.success ]];then\n";
	print FILE "echo \"Super reads failed, check super2.err and files in ./work2/\"\n";
	print FILE "exit\n";
	print FILE "fi\n";

#now, using read positions in super reads, we find out which mates got joined -- these are the ones that do not have the biotin in the middle, call them chimeric
	if(not(-e "chimeric_sj.txt")||$rerun_pe==1||$rerun_sj==1){
	    print FILE "filter_alt.pl outtie < work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt >  chimeric_sj.txt \n";
	    $rerun_sj=1;
	}

#we also do initial redundancy filtering here, based on positions of reads in suoer reads
	if(not(-e "redundant_sj.txt")||$rerun_pe==1||$rerun_sj==1){
	    print FILE "filter_redundancy.pl 2 < work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt > redundant_sj.txt\n";
	    $rerun_sj=1;
	} 

	print FILE "echo 'Chimeric/Redundant jump reads:';wc -l  chimeric_sj.txt redundant_sj.txt;\n";

#remove all chimeric and all redundant reads from sj.cor.fa
	if(not(-e "sj.cor.clean.rev.fa")||$rerun_pe==1||$rerun_sj==1){
	    print FILE "extractreads.pl <(cat chimeric_sj.txt redundant_sj.txt | perl -e '{
		while(\$line=<STDIN>){
		chomp(\$line);
		\$h{\$line}=1
		}
		open(FILE,\$ARGV[0]);
		while(\$line=<FILE>){
		chomp(\$line);
		print \$line,\"\\n\" if(not(defined(\$h{\$line})));
		}
		}' <(awk '{
			prefix=substr(\$1,1,2); 
			readnumber=int(substr(\$1,3));  
			if(readnumber\%2==0){
				last_readnumber=readnumber; 
				last_prefix=prefix;
			}else{
				if(last_readnumber==readnumber-1 && last_prefix==prefix){
					print prefix\"\"last_readnumber\"\\n\"prefix\"\"readnumber;
				}
			}
			}' work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt)) sj.cor.fa 1 > sj.cor.clean.fa\n";
	    $rerun_sj=1;
	    print FILE "rm -rf sj.cor.clean.rev.fa\n";
	    for($i=0;$i<scalar(@jump_info_array);$i++){
        		@f=split(/\s+/,$jump_info_array[$i]);
	            	my $if_innie="";
        	        $if_innie=" | reverse_complement " if($f[1]>0);
            		print FILE "grep --text -A 1 '^>$f[0]' sj.cor.clean.fa | grep --text -v '^\\-\\-' $if_innie >> sj.cor.clean.rev.fa\n";
        		}
#here we perform another round of filtering bad mates
	print FILE "cat sj.cor.clean.rev.fa | putReadsIntoGroupsBasedOnSuperReads --super-read-sequence-file work2/superReadSequences.fasta --read-placements-file work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt > sj.cor.clean.rev.fa.bak\n";
	print FILE "rm -f sj.cor.clean.rev.fa\n";
	print FILE "findReversePointingJumpingReads.perl -s \$JF_SIZE --Celera-terminator-directory . --jumping-library-read-file sj.cor.clean.rev.fa.bak --reads-file pe.cor.fa --output-directory work4 --min-kmer-len 23 --max-kmer-len 80 --num-threads $NUM_THREADS --maxnodes 1000 --reduce-read-set-kmer-size 27 --max-reads-in-memory 100000000 --faux-insert-mean 500 --faux-insert-stdev 100 --num-joins-per-directory 101 1>findReversePointingJumpingReads.err 2>&1 \n";
	print FILE "extractreads_not.pl work4/output.txt sj.cor.clean.rev.fa.bak 1 > sj.cor.clean.rev.fa\n";
	print FILE "echo Found extra chimeric mates: \n";
	print FILE "wc -l work4/output.txt\n";
	print FILE "rm -rf work4/readFile.??? work4/workReadsVsFaux work4/workFauxVsFaux\n";
	}

#here we extend the jumping library reads if they are too short
        if(not(-e "sj.cor.ext.fa")||$rerun_pe==1||$rerun_sj==1){      
        $rerun_sj=1;
	if($EXTEND_JUMP_READS==1){
        print FILE "createSuperReadsForDirectory.perl -jumplibraryreads -minreadsinsuperread 1 -l \$KMER_J -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.sj.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.jump.fasta -t $NUM_THREADS -mikedebug work3 sj.cor.clean.rev.fa 1> super2.err 2>&1\n";

#check if the super reads pipeline finished successfully
        print FILE "if [[ ! -e work3/superReads.success ]];then\n";
        print FILE "echo \"Super reads failed, check super2.err and files in ./work2/\"\n";
        print FILE "exit\n";
        print FILE "fi\n";

        print FILE "ln -sf work3/superReadSequences.jumpLibrary.fasta sj.cor.ext.fa\n";
	}else{
	print FILE "ln -sf sj.cor.clean.rev.fa sj.cor.ext.fa\n";
	}
	}

#here we create the frg files for CA from the jump libraries: each jump library will contribute one jump frg file and one additional frg file of linking information from "chimers"
	print FILE "echo -n 'creating FRG files ';date;\n";
	print FILE "rm -rf compute_jump_coverage.txt\n";
	
	for($i=0;$i<scalar(@jump_info_array);$i++){
	    @f=split(/\s+/,$jump_info_array[$i]);
	    print FILE "echo -n \"$f[1] \" >> compute_jump_coverage.txt\n";
	    print FILE "grep --text -A 1 '^>$f[0]' sj.cor.ext.fa | grep --text -v '^\\-\\-' > $f[0].tmp\n";
	    print FILE "error_corrected2frg $f[0] ",abs($f[1])," $f[2] 2000000000 $f[0].tmp | grep --text '^{LKG' |wc -l >> compute_jump_coverage.txt\n";
        }
	print FILE "JUMP_BASES_COVERED=`awk 'BEGIN{b=0}{b+=\$1*\$2;}END{print b}' compute_jump_coverage.txt`\n";
        print FILE "save JUMP_BASES_COVERED\n";

#here we reduce jump library coverage: we know the genome size (from k-unitigs) and JUMP_BASES_COVERED contains total jump library coverage :)
	print FILE "perl -e '{srand(1);\$cov=int('\$JUMP_BASES_COVERED'/'\$ESTIMATED_GENOME_SIZE'); print \"JUMP insert coverage: \$cov\\n\"; \$optimal_cov=$LIMIT_JUMP_COVERAGE;if(\$cov>\$optimal_cov){print \"Reducing JUMP insert coverage from \$cov to \$optimal_cov\\n\";\$prob_coeff=\$optimal_cov/\$cov;open(FILE,\"gkp.edits.msg\");while(\$line=<FILE>){chomp(\$line);\@f=split(/\\s+/,\$line);\$deleted{\$f[2]}=1;}close(FILE); open(FILE,\"sj.cor.clean.fa\");while(\$line=<FILE>){next if(not(\$line =~ /^>/));chomp(\$line);if(int(substr(\$line,3))%2==0) {print STDERR substr(\$line,1),\"\\n\" if(rand(1)>\$prob_coeff);}}}}' 2>mates_to_break.txt\n";
	print FILE "extractreads_not.pl mates_to_break.txt sj.cor.ext.fa 1 >  sj.cor.ext.reduced.fa\n";
	for($i=0;$i<scalar(@jump_info_array);$i++){
	    @f=split(/\s+/,$jump_info_array[$i]);
	    $list_of_frg_files.="$f[0].cor.clean.frg ";
	    print FILE "grep --text -A 1 '^>$f[0]' sj.cor.ext.reduced.fa | grep --text -v '^\\-\\-' > $f[0].tmp\n";
	    print FILE "error_corrected2frg $f[0] ",abs($f[1])," $f[2] 2000000000 $f[0].tmp > $f[0].cor.clean.frg\n";
	    print FILE "rm -f $f[0].tmp\n";
	}
    }

###done with super reads and filtering for jump###
    print FILE "\n";
##super reads for PE###

    print FILE "echo -n 'computing super reads from PE ';date;\n";

#create super reads from PE    
    if($rerun_pe==1|| not(-e "work1")){
	print FILE "rm -rf work1\n";
	$rerun_pe=1;
    }

    print FILE "createSuperReadsForDirectory.perl -l \$KMER -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.pe.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.fasta -t $NUM_THREADS -mikedebug work1 pe.cor.fa 1> super1.err 2>&1\n";

#check if the super reads pipeline finished successfully
    print FILE "if [[ ! -e work1/superReads.success ]];then\n";
    print FILE "echo \"Super reads failed, check super1.err and files in ./work1/\"\n";
    print FILE "exit\n";
    print FILE "fi\n";

    if($USE_LINKING_MATES==1){
#now we extract those PE mates that did not end up in the same super read -- we call them linking mates, they will be useful for scaffolding
	if(not(-e "pe.linking.fa")||$rerun_pe==1){
	    print FILE "extractreads.pl <( awk 'BEGIN{last_readnumber=-1;last_super_read=\"\"}{readnumber=int(substr(\$1,3));if(readnumber%2>0){readnumber--}super_read=\$2;if(readnumber==last_readnumber){if(super_read!=last_super_read){print read;print \$1;}}else{read=\$1;last_super_read=\$2}last_readnumber=readnumber}' work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt )  pe.cor.fa 1 > pe.linking.fa\n";
	    $rerun_pe=1;
	}

#create frg files for PE data
	foreach $v(@pe_info_array){
	    @f=split(/\s+/,$v);
	    $list_of_frg_files.="$f[0].linking.frg ";
	    if(not(-e "$f[0].linking.frg")||$rerun_pe==1){
		print FILE "grep --text -A 1 '^>$f[0]' pe.linking.fa | grep --text -v '^\\-\\-' > $f[0].tmp\n";
		print FILE "error_corrected2frg $f[0] $f[1] $f[2] 2000000000 $f[0].tmp > $f[0].linking.frg\n";
		print FILE "rm $f[0].tmp\n";
	    }
	}
	print FILE "echo -n 'Linking PE reads ';\ncat ??.linking.frg |grep --text '^{FRG' |wc -l;\n";
    }

#create frg file for super reads
    if(not(-e "superReadSequences_shr.frg")||$rerun_pe==1){
	print FILE "cat work1/superReadSequences.fasta | create_sr_frg.pl 2>renamed_sr.txt | fasta2frg.pl sr >  superReadSequences_shr.frg\n";
    }

###done with super reads for PE###
    print FILE "\n";
}

###Celera Assembler###
if(not(-e "CA/9-terminator/genome.qc")|| $rerun_pe || $rerun_sj){
    print FILE "\necho -n 'Celera Assembler ';date;\n";
    if($rerun_sj==1||$rerun_pe==1){
	print FILE "rm -rf CA\n";
    }

#this if statement is here because if OTHER frg is specified, we will have to do OBT+ECR, it will slow us down, but it has to be done :(
        if(scalar(@other_info_array)>0){
            $ovlMerSize=22;
            $other_parameters="doFragmentCorrection=1 doOverlapBasedTrimming=1 doExtendClearRanges=2 ovlMerSize=22";
        }else{
            $ovlMerSize=30;
            $other_parameters="doFragmentCorrection=0 doOverlapBasedTrimming=0 doExtendClearRanges=0 ovlMerSize=30";
        }

    if(not(-d "CA/7-0-CGW")|| $rerun_pe || $rerun_sj){
###figure out the optimal parameters for CA###
	print FILE "TOTAL_READS=`cat  *.frg |grep --text '^{FRG'|wc -l`\n";
        print FILE "save TOTAL_READS\n";
	print FILE "ovlRefBlockSize=`perl -e '\$s=int('\$TOTAL_READS'/8); if(\$s>100000){print \$s}else{print \"100000\"}'`\n";
        print FILE "save ovlRefBlockSize\n";
	print FILE "ovlHashBlockSize=`perl -e '\$s=int('\$TOTAL_READS'/80); if(\$s>10000){print \$s}else{print \"10000\"}'`\n";
        print FILE "save ovlHashBlockSize\n";
	print FILE "ovlCorrBatchSize=\$ovlHashBlockSize\n";
        print FILE "save ovlCorrBatchSize\n";
####done figuring out CA parameters###

#estimating mer threshold for overlapper to cover 90% of all distinct k-mers
	if(not(-d "CA/genome.ovlStore")|| $rerun_pe || $rerun_sj){
	    print FILE "ovlMerThreshold=`jellyfish-2.0 histo -t $NUM_THREADS k_u_hash_0 | awk '{thresh=75;if(\$1>1) {dist+=\$2;if(dist>int(\"'\$ESTIMATED_GENOME_SIZE'\")*0.98&&flag==0){if(\$1>thresh) thresh=\$1;flag=1}}}END{print thresh}'`\n";
	    print FILE "echo ovlMerThreshold=\$ovlMerThreshold\n\n";
	}

#filter out non-junction fragments

	if(scalar(@jump_info_array)>0){
	    if(not(-e "CA/4-unitigger-filter/gkp.edits.msg") || $rerun_pe || $rerun_sj){
		print FILE "runCA  gkpFixInsertSizes=0 jellyfishHashSize=\$JF_SIZE ovlRefBlockSize=\$ovlRefBlockSize ovlHashBlockSize=\$ovlHashBlockSize ovlCorrBatchSize=\$ovlCorrBatchSize utgErrorRate=0.03 merylMemory=8192 ovlMemory=4GB stopAfter=unitigger ovlMerThreshold=\$ovlMerThreshold bogBreakAtIntersections=0 unitigger=bog bogBadMateDepth=1000000 -p genome -d CA merylThreads=$NUM_THREADS frgCorrThreads=1 frgCorrConcurrency=$NUM_THREADS cnsConcurrency=$NUM_THREADS ovlCorrConcurrency=$NUM_THREADS ovlConcurrency=$NUM_THREADS ovlThreads=1 $other_parameters superReadSequences_shr.frg $list_of_frg_files  1> runCA0.out 2>&1\n\n";


		print FILE "if [[ ! -e CA/4-unitigger/unitigger.err ]];then\n";
		print FILE "echo \"CA failed, check output under CA/ and runCA0.out\"\n";
		print FILE "exit\n";
		print FILE "fi\n";

		print FILE "cd CA/\nmv 4-unitigger 4-unitigger-filter\ncd 4-unitigger-filter\ngrep --text '^>' ../../sj.cor.ext.reduced.fa |awk '{print substr(\$1,2)}' > sj.uid\nfilter_library.sh ../ genome sj.uid\n";
		print FILE "cat genome.chimeric.uid |awk '{print \"frg uid \"\$1\" mateiid 0\"}'  > gkp.edits.msg\n";
		print FILE "echo -n \"Found additional non-junction reads: \"\nwc -l gkp.edits.msg\n";
		print FILE "gatekeeper --edit gkp.edits.msg ../genome.gkpStore 1>gatekeeper.err 2>&1\n";
		print FILE "cd ../\nrm -rf *.tigStore\ncd ../\n\n";
		print FILE "\n";
	    }
	    print FILE "if [[ -e \"CA/4-unitigger-filter/gkp.edits.msg\" ]];then\n";
	    print FILE "echo \"Filter success\"\n";
	    print FILE "else\n";
	    print FILE "echo \"Filtering failed, check output under CA/4-unitigger-filter/ and runCA0.out\"\n";
	    print FILE "exit\n";
	    print FILE "fi\n";
	}

	if(not(-e "CA/5-consensus/consensus.success")|| $rerun_pe || $rerun_sj){
	    print FILE "runCA ovlMerThreshold=\$ovlMerThreshold gkpFixInsertSizes=0 $CA_PARAMETERS jellyfishHashSize=\$JF_SIZE ovlRefBlockSize=\$ovlRefBlockSize ovlHashBlockSize=\$ovlHashBlockSize ovlCorrBatchSize=\$ovlCorrBatchSize stopAfter=consensusAfterUnitigger unitigger=bog -p genome -d CA merylThreads=$NUM_THREADS frgCorrThreads=1 frgCorrConcurrency=$NUM_THREADS cnsConcurrency=$NUM_THREADS ovlCorrConcurrency=$NUM_THREADS ovlConcurrency=$NUM_THREADS ovlThreads=1 $other_parameters superReadSequences_shr.frg $list_of_frg_files   1> runCA1.out 2>&1\n";

	    print FILE "if [[ -e \"CA/4-unitigger/unitigger.err\" ]];then\n";
	    print FILE "echo \"Overlap/unitig success\"\n";
	    print FILE "else\n";
	    print FILE "echo \"Overlap/unitig failed, check output under CA/ and runCA1.out\"\n";
	    print FILE "exit\n";
	    print FILE "fi\n";

#we now recompute the A-stat for the unitigs based on positions of PE reads in the super-reads
	    print FILE "recompute_astat_superreads.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt\n";

#here we filter for repetitive kmers in the unique unitigs
	    print FILE "NUM_SUPER_READS=`cat superReadSequences_shr.frg $tmplist | grep --text '^{FRG' |wc -l`\n";
            print FILE "save NUM_SUPER_READS\n";
	    print FILE "cd CA\n";
	    print FILE "tigStore -g genome.gkpStore -t genome.tigStore 2 -d layout -U | tr -d '-' | awk 'BEGIN{print \">unique unitigs\"}{if(\$1 == \"cns\"){seq=\$2}else if(\$1 == \"data.unitig_coverage_stat\" && \$2>=5){print seq\"N\"}}' | jellyfish-2.0 count -L 2 -C -m $ovlMerSize -s \$ESTIMATED_GENOME_SIZE -t $NUM_THREADS -o unitig_mers /dev/fd/0\n";
	    print FILE "cat <(overlapStore -b 1 -e \$NUM_SUPER_READS -d genome.ovlStore  | awk '{if(\$1<'\$NUM_SUPER_READS' && \$2<'\$NUM_SUPER_READS') print \$0}'|filter_overlap_file -t $NUM_THREADS <(gatekeeper -dumpfastaseq genome.gkpStore ) unitig_mers /dev/fd/0) <(overlapStore -d genome.ovlStore | awk '{if(\$1>='\$NUM_SUPER_READS' || \$2>='\$NUM_SUPER_READS') print \$1\" \"\$2\" \"\$3\" \"\$4\" \"\$5\" \"\$6\" \"\$7}')  |convertOverlap -b -ovl > overlaps.ovb\n";
	    print FILE "rm -rf 4-unitigger 5-consensus genome.tigStore genome.ovlStore\n";
	    print FILE "overlapStore -c genome.ovlStore -M 4096 -t $NUM_THREADS -g genome.gkpStore overlaps.ovb 1>overlapstore.err 2>&1\n";
	    print FILE "cd ..\n";
	    
#and now we rerun the assembler
	    print FILE "runCA ovlMerThreshold=\$ovlMerThreshold gkpFixInsertSizes=0 $CA_PARAMETERS jellyfishHashSize=\$JF_SIZE ovlRefBlockSize=\$ovlRefBlockSize ovlHashBlockSize=\$ovlHashBlockSize ovlCorrBatchSize=\$ovlCorrBatchSize stopAfter=consensusAfterUnitigger unitigger=bog -p genome -d CA merylThreads=$NUM_THREADS frgCorrThreads=1 frgCorrConcurrency=$NUM_THREADS cnsConcurrency=$NUM_THREADS ovlCorrConcurrency=$NUM_THREADS ovlConcurrency=$NUM_THREADS ovlThreads=1 $other_parameters superReadSequences_shr.frg $list_of_frg_files   1> runCA1.out 2>&1\n";

	    print FILE "if [[ -e \"CA/4-unitigger/unitigger.err\" ]];then\n";
	    print FILE "echo \"Overlap/unitig success\"\n";
	    print FILE "else\n";
	    print FILE "echo \"Overlap/unitig failed, check output under CA/ and runCA1.out\"\n";
	    print FILE "exit\n";
	    print FILE "fi\n";
#now we check if the unitig consensus which is sometimes problematic, failed, and fix the unitigs

	    print FILE "if [[ -e \"CA/5-consensus/consensus.success\" ]];then\n";
	    print FILE "echo \"Unitig consensus success\"\n";
	    print FILE "else\n";
	    print FILE "echo \"Fixing unitig consensus...\"\n";
	    print FILE "mkdir CA/fix_unitig_consensus\n";
	    print FILE "cd CA/fix_unitig_consensus\n";
	    print FILE "cp `which fix_unitigs.sh` .\n";
	    print FILE "./fix_unitigs.sh genome \n";
	    print FILE "cd ../../\n";
	    print FILE "fi\n";

#we now recompute the A-stat for the unitigs based on positions of PE reads in the super-reads
	    print FILE "recompute_astat_superreads.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt\n";
	}
    }

#and we continue into the scaffolder...
    print FILE "runCA $CA_PARAMETERS unitigger=bog -p genome -d CA cnsConcurrency=$NUM_THREADS computeInsertSize=0 $other_parameters 1>runCA2.out 2>&1\n";

    print FILE "if [[ -e \"CA/9-terminator/genome.qc\" ]];then\n";
    print FILE "echo \"CA success\"\n";
    print FILE "else\n";
    print FILE "echo \"CA failed, check output under CA/ and runCA2.out\"\n";
    print FILE "exit\n";
    print FILE "fi\n";
}
#here we close gaps in scaffolds:  we use create_k_unitigs allowing to continue on count 1 sequence and then generate fake reads from the 
#end sequences of contigs that are next to each other in scaffolds, and then use super reads software to close the gaps for k=17...31
if($CLOSE_GAPS){

    print FILE "echo -n 'Gap closing ';date;\n";
    my $reads_argument="";
    @f=split(" ",$list_pe_files);
    foreach $v(@f){
	$reads_argument.="--reads-file $v ";
    }
    @f=split(" ",$list_jump_files);
    foreach $v(@f){
        $reads_argument.="--reads-file $v ";
    }

    print FILE "closeGapsLocally.perl -s $JF_SIZE --Celera-terminator-directory CA/9-terminator $reads_argument --output-directory CA/10-gapclose --min-kmer-len 17 --max-kmer-len \$((\$PE_AVG_READ_LENGTH-5)) --num-threads $NUM_THREADS --contig-length-for-joining \$((\$PE_AVG_READ_LENGTH-1)) --contig-length-for-fishing 200 --reduce-read-set-kmer-size 21 1>gapClose.err 2>&1\n";
    print FILE "if [[ -e \"CA/10-gapclose/genome.ctg.fasta\" ]];then\n";
    print FILE "echo \"Gap close success. Output sequence is in CA/10-gapclose/genome.\{ctg,scf\}.fasta\"\n";
    print FILE "else\n";
    print FILE "echo \"Gap close failed, you can still use pre-gap close files under CA/9-terminator/. Check gapClose.err for problems.\"\n";
    print FILE "exit\n";
    print FILE "fi\n";
}

###Done !!!! Hoorayyyy!!! :)###
print FILE "echo -n 'All done ';date;\n";

close(FILE);
system("chmod 0755 assemble.sh");
print "done.\nexecute assemble.sh to run the assembly\n";

