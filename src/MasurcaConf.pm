package MasurcaConf;

use strict;
use warnings;

our (@ISA, @EXPORT, @EXPORT_OK);

BEGIN {
  require Exporter;
  @ISA = qw(Exporter);
  @EXPORT_OK = qw(&fail $default_config $config_file);
}

my $default_values = {
  # Parameters
  EXTEND_JUMP_READS       => 0,
  USE_GRID                => 0,
  GRID_QUEUE              => "all.q",
  GRID_BATCH_SIZE         => 5000000000,
  GRID_ENGINE             => "SGE",
  LHE_COVERAGE            => 25,
  JF_SIZE                 => 100000000,
  KMER                    => "auto",
  KMER_COUNT_THRESHOLD    => 1,
  KMER_RELIABLE_THRESHOLD => 3,
  TRIM_PARAM              => 3,
  RELIABLE_Q_PARAM	  => 40,
  NUM_THREADS             => 2,
  NUM_CNS_THREADS         => 2,
  LIMIT_JUMP_COVERAGE     => 300,
  MEGA_READS_ONE_PASS     => "",
  USE_LINKING_MATES       => 0,
  DO_HOMOPOLYMER_TRIM     => 0,
  CLOSE_GAPS              => 1,
  NO_MMAP                 => 1,
  STOP_AFTER_SR           => 0,
  CA_PARAMETERS           => "",
  SOAP_ASSEMBLY           => 0,
  FLYE_ASSEMBLY           => 0,

  # Data
  PE_INFO       => [],
  JUMP_INFO     => [],
  OTHER_INFO    => [],
  MOLECULO_INFO => [],
  PACBIO_INFO   => [],
  NANOPORE_INFO => [],
  NANOPORE_RNA_INFO => [],
  REF_INFO => [],
  POLISH_INFO => [],
};

our $config_file;

sub can_read {
  my ($file) = @_;

  my $res = open(my $io, "<", $file);
  close($io);
  return $res;
}

sub default_config {
  return <<'EOS';
# example configuration file 

# DATA is specified as type {PE,JUMP,OTHER,PACBIO} and 5 fields:
# 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads
# 5)fastq(.gz)_rev_reads. The PE reads are always assumed to be
# innies, i.e. --->.<---, and JUMP are assumed to be outties
# <---.--->. If there are any jump libraries that are innies, such as
# longjump, specify them as JUMP and specify NEGATIVE mean. Reverse reads
# are optional for PE libraries and mandatory for JUMP libraries. Any
# OTHER sequence data (454, Sanger, Ion torrent, etc) must be first
# converted into Celera Assembler compatible .frg files (see
# http://wgs-assembler.sourceforge.com)
DATA
#Illumina paired end reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
#if single-end, do not specify <reverse_reads>
#If mean/stdev are unknown use 500 and 50 -- these are safe values that will work for most runs
#MUST HAVE Illumina paired end reads to use MaSuRCA
PE= pe 500 50  /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq
#Illumina mate pair reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
#JUMP= sh 3600 200  /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
#pacbio OR nanopore reads must be in a single fasta or fastq file with absolute path, can be gzipped
#if you have both types of reads supply them both as NANOPORE type
#PACBIO=/FULL_PATH/pacbio.fa
#NANOPORE=/FULL_PATH/nanopore.fa
#Legacy reads (Sanger, 454, etc) in one frg file, concatenate your frg files into one if you have many
#OTHER=/FULL_PATH/file.frg
#synteny-assisted assembly, concatenate all reference genomes into one reference.fa; works for Illumina-only data
#REFERENCE=/FULL_PATH/nanopore.fa
END

PARAMETERS
#PLEASE READ all comments to essential parameters below, and set the parameters according to your project
#set this to 1 if your Illumina mate pair (jumping) library reads are shorter than 100bp
EXTEND_JUMP_READS=0
#this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE = auto
#set this to 1 for all Illumina-only assemblies
#set this to 0 if you have more than 15x coverage by long reads (Pacbio or Nanopore) or any other long reads/mate pairs (Illumina MP, Sanger, 454, etc)
USE_LINKING_MATES = 0
#specifies whether to run the assembly on the grid
USE_GRID=0
#specifies grid engine to use SGE or SLURM
GRID_ENGINE=SGE
#specifies queue (for SGE) or partition (for SLURM) to use when running on the grid MANDATORY
GRID_QUEUE=all.q
#batch size in the amount of long read sequence for each batch on the grid
GRID_BATCH_SIZE=500000000
#use at most this much coverage by the longest Pacbio or Nanopore reads, discard the rest of the reads
#can increase this to 30 or 35 if your long reads reads have N50<7000bp
LHE_COVERAGE=25
#this parameter is useful if you have too many Illumina jumping library reads. Typically set it to 60 for bacteria and 300 for the other organisms 
LIMIT_JUMP_COVERAGE = 300
#these are the additional parameters to Celera Assembler; do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
#CABOG ASSEMBLY ONLY: set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS =  cgwErrorRate=0.15
#CABOG ASSEMBLY ONLY: whether to attempt to close gaps in scaffolds with Illumina  or long read data
CLOSE_GAPS=1
#number of cpus to use, set this to the number of CPUs/threads per node you will be using
NUM_THREADS = 32
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*20
JF_SIZE = 200000000
#ILLUMINA ONLY. Set this to 1 to use SOAPdenovo contigging/scaffolding module.  
#Assembly will be worse but will run faster. Useful for very large (>=8Gbp) genomes from Illumina-only data
SOAP_ASSEMBLY=0
#If you are doing Hybrid Illumina paired end + Nanopore/PacBio assembly ONLY (no Illumina mate pairs or OTHER frg files).  
#Set this to 1 to use Flye assembler for final assembly of corrected mega-reads.  
#A lot faster than CABOG, AND QUALITY IS THE SAME OR BETTER.   
#DO NOT use if you have less than 20x coverage by long reads.
FLYE_ASSEMBLY=1
END
EOS
}

# Like die, but don't print backtrace like information
sub fail {
  my ($msg, $line) = @_;
  print(STDERR "Error");
  if(defined($config_file) &&  defined($line)) {
    print(STDERR " line $_[1] of configuration file '$config_file':\n");
  } else {
    print(STDERR ": ");
  }
  chomp($msg);
  print(STDERR $msg, "\n");
  exit(1);
}


# Match a line of the form "KEY = VALUE". Returns undef if it does not match
sub read_param {
  my ($line) = @_;
  fail("Can't parse line '$line'", $.) unless $line =~ /^\s*(\w+)\s*=\s*(.*?)\s*$/;
  my ($k, $v) = ($1, $2);
  chomp($v);
  return ($k, $v);
}

my %used_library_ids;


sub parse_parameters {
  my ($key, $param, $res) = @_;

  if($key eq "EXTEND_JUMP_READS"){
    fail("bad value for EXTEND_JUMP_READS, enter 0 or 1", $.) unless($param =~ /^[01]$/);
    $$res{EXTEND_JUMP_READS} = int($param);
  } elsif($key eq "DO_HOMOPOLYMER_TRIM"){
    fail("bad value for DO_HOMOPOLYMER_TRIM, enter 0 or 1", $.) unless($param =~ /^[01]$/);
    $$res{DO_HOMOPOLYMER_TRIM} = int($param);
  } elsif($key eq "TRIM_PARAM") {
    fail("bad value for TRIM_PARAM, it should be a positive integer", $.) unless $param =~ /^\d*$/;
    $$res{TRIM_PARAM} = length($param) > 0 ? int($param) : 2;
  } elsif($key eq "RELIABLE_Q_PARAM") {
    fail("bad value for RELIABLE_Q_PARAM, it should be a positive integer", $.) unless $param =~ /^\d*$/;
    $$res{RELIABLE_Q_PARAM} = length($param) > 0 ? int($param) : 35;
  } elsif($key eq "CLOSE_GAPS"){
    fail("bad value for CLOSE_GAPS, enter 0 or 1", $.) unless($param =~ /^[01]$/);
    $$res{CLOSE_GAPS} = int($param);
  } elsif($key eq "MEGA_READS_ONE_PASS"){
    fail("bad value for MEGA_READS_ONE_PASS, enter 0 or 1", $.) unless($param =~ /^[01]$/);
    $$res{MEGA_READS_ONE_PASS} = "--onepass" if(not($param==0)) ;
  } elsif($key eq "USE_GRID"){
    fail("bad value for USE_GRID, enter 0 or 1", $.) unless($param =~ /^[01]$/);
    $$res{USE_GRID} = int($param);
  } elsif($key eq "LHE_COVERAGE"){
    fail("bad value for LHE_COVERAGE, it should be a positive integer", $.) unless($param =~ /^\d*$/);
    $$res{LHE_COVERAGE} = int($param);
  } elsif($key eq "CA_PARAMETERS"){
    fail("bad value for CA_PARAMETERS", $.) if($param eq "");
    $$res{CA_PARAMETERS} = $param;
  } elsif($key eq "GRID_QUEUE"){
    fail("bad value for GRID_QUEUE", $.) if($param eq "");
    $$res{GRID_QUEUE} = $param;
  } elsif($key eq "GRID_BATCH_SIZE"){
    fail("bad value for GRID_BATCH_SIZE", $.) unless($param =~ /^\d*$/);
    $$res{GRID_BATCH_SIZE} = $param;
  } elsif($key eq "GRID_ENGINE"){
    fail("bad value for GRID_ENGINE, only SGE or SLURM are allowed", $.) unless($param eq "SGE" || $param eq "SLURM" || $param eq "MANUAL");
    $$res{GRID_ENGINE} = $param; 
  } elsif($key eq "LIMIT_JUMP_COVERAGE"){
    fail("bad value for LIMIT_JUMP_COVERAGE, enter a number > 1", $.) if($param<=1);
    $$res{LIMIT_JUMP_COVERAGE} = int($param);
  } elsif($key eq "GRAPH_KMER_SIZE"){
    fail("bad value for GRAPH_KMER_SIZE, enter auto or number >= 15 and <= 151", $.) if(not($param eq "auto") && ($param<15 || $param>151));
    $$res{KMER} = $param;
  } elsif($key eq "USE_LINKING_MATES"){
    fail("bad value for USE_LINKING_MATES, enter 0 or 1", $.) unless $param =~ /^[01]$/;
    $$res{USE_LINKING_MATES} = int($param);
  } elsif($key eq "KMER_COUNT_THRESHOLD"){
    fail("bad value for KMER_COUNT_THRESHOLD. Enter a number >= 1", $.) if($param<1);
    $$res{KMER_COUNT_THRESHOLD} = int($param);
    $$res{KMER_RELIABLE_THRESHOLD} = 3*int($param)
  } elsif($key eq "NUM_THREADS"){
    fail("bad value for NUM_THREADS. Enter a number >= 1", $.) if($param<1);
    $$res{NUM_THREADS} = int($param);
    $$res{NUM_CNS_THREADS} = int($param/4)+1;
  } elsif($key eq "JF_SIZE"){
    fail("bad value for JF_SIZE, enter a number >= 100000", $.) if($param<100000);
    $$res{JF_SIZE} = int($param);
  } elsif($key eq "NO_MMAP"){
    fail("bad value for NO_MMAP, enter 0 or 1", $.) unless($param =~ /^[01]$/);
    $$res{NO_MMAP} = int($param);
  } elsif($key eq "STOP_AFTER_SUPERREADS"){
    fail("bad value for STOP_AFTER_SUPERREADS, enter 0 or 1", $.) unless($param =~ /^[01]$/);
    $$res{STOP_AFTER_SR} = int($param);
  } elsif($key eq "SOAP_ASSEMBLY") {
    fail("bad value for SOAP_ASSEMBLY, enter 0 or 1", $.) unless($param =~ /^[01]$/);
    $$res{SOAP_ASSEMBLY} = int($param);   
  } elsif($key eq "FLYE_ASSEMBLY") {
    fail("bad value for FLYE_ASSEMBLY, enter 0 or 1", $.) unless($param =~ /^[01]$/);
    $$res{FLYE_ASSEMBLY} = int($param);
  } else {
    return 1;
  }
  return 0;
}

sub parse_data {
  my ($key, $param, $res) = @_;
  
  if($key eq "PE"){
    my @f=split(" ", $param);
    fail("improper id for PE library '$f[0]'. It should be two character long (like 'p0')", $.) if(not(length($f[0])==2));
    fail("duplicate id for PE library '$f[0]'", $.) if(defined($used_library_ids{$f[0]}));
    $used_library_ids{$f[0]}=1;
    fail("improper mean '$f[1]' for PE library '$f[0]'. It must be a positive number", $.) unless(int($f[1])>0);
    fail("improper stdev '$f[2]' for PE library '$f[0]'. It must be a positive number", $.) unless(int($f[2])>0);
    can_read($f[3]) or fail("invalid forward file for PE library '$f[0]': '$f[3]' $!", $.);
    if(defined($f[4])){
      can_read($f[4]) or fail("invalid reverse file for PE library '$f[0]': '$f[4]' $!", $.);
    } else {
      push(@f, $f[3]);
    }
    push(@{$$res{PE_INFO}}, \@f);
  } elsif($key eq "JUMP"){
    my @f = split(" ", $param);
    fail("improper id for JUMP library '$f[0]'. It should be two character long (like 'j1')", $.) if(not(length($f[0])==2));
    fail("duplicate id for JUMP library '$f[0]'", $.) if(defined($used_library_ids{$f[0]}));
    $used_library_ids{$f[0]}=1;
    fail("improper mean '$f[1]' for JUMP library '$f[0]'. It must be a number strictly greater or less than zero", $.) if(int($f[1])==0);
    fail("improper stdev '$f[2]' for JUMP library '$f[0]'. It must be a positive number", $.) unless(int($f[2])>0);
    can_read($f[3]) or fail("invalid forward file for JUMP library '$f[0]': '$f[3]' $!", $.);
    can_read($f[4]) or fail("invalid reverse file for JUMP library '$f[0]': '$f[4]' $!", $.);
    push(@{$$res{JUMP_INFO}}, \@f);
  } elsif($key eq "OTHER"){
    fail("incorrect frg file name '$param'. It must end in '.frg'", $.) unless($param =~/\.frg$/);
    can_read($param) or fail("invalid frg file for OTHER: '$param' $!", $.);
    push(@{$$res{OTHER_INFO}}, $param);
  } elsif($key eq "MOLECULO") {
    can_read($param) or fail("invalid file for MOLECULO: '$param' $!", $.);
    push(@{$$res{MOLECULO_INFO}}, $param);
  } elsif($key eq "NANOPORE") {
    can_read($param) or fail("invalid file for NANOPORE: '$param' $!", $.);
    push(@{$$res{NANOPORE_INFO}}, $param);
  } elsif($key eq "NANOPORE_RNA") {
    can_read($param) or fail("invalid file for NANOPORE_RNA: '$param' $!", $.);
    push(@{$$res{NANOPORE_RNA_INFO}}, $param);
  } elsif($key eq "REFERENCE") {
    can_read($param) or fail("invalid file for REFERENCE: '$param' $!", $.);
    push(@{$$res{REF_INFO}}, $param);
  } elsif($key eq "POLISH") {
    can_read($param) or fail("invalid file for POLISH: '$param' $!", $.);
    push(@{$$res{POLISH_INFO}}, $param);
  } elsif($key eq "PACBIO") {
    can_read($param) or fail("invalid file for PACBIO: '$param' $!", $.);
    push(@{$$res{PACBIO_INFO}}, $param);
  } else {
    return 1;
  }
  return 0;
}

sub parse {
  ($config_file) = @_;
  my ($in_parameters, $in_data) = (0, 0);
  
  open(FILE, "<", $config_file) or fail("Can't open config file '$config_file': $!");
  my %res = %$default_values;
  while(my $line=<FILE>){
    chomp($line);
    next if($line =~ /^\s*(#|$)/);
    if($line =~ /^DATA\s*$/){
      fail("error in config file: mixed PARAMETERS and DATA", $.) if $in_parameters;
      fail("duplicate DATA header", $.) if $in_data;
      $in_data = 1;
      next;
    } elsif($line =~ /^PARAMETERS\s*$/){
      fail("error in config file: mixed PARAMETERS and DATA", $.) if $in_data;
      fail("duplicate PARAMETERS header", $.) if $in_parameters;
      $in_parameters=1;
      next;
    } elsif($line =~ /^PATHS\s*$/) {
      warn("PATHS section is obsolete. You should remove it from your configuration file. Skipping to next section...");
      while(<FILE> !~ /^END\s*/) { }
      next;
    } elsif($line =~ /^END\s*$/) {
      fail("Unexpected 'END' keyword. Not in PARAMETERS or DATA section", $.) if(!($in_parameters || $in_data));
      $in_data = $in_parameters = 0;
      next;
    }

    my ($key, $param) = read_param($line);
    my $error = 1;
    $error = parse_parameters($key, $param, \%res) if($in_parameters==1);
    $error = parse_data($key, $param, \%res) if($in_data==1);
    fail("Invalid line '$line'", $.) if $error;
  }
  return %res;
}

1;
