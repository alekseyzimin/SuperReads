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
  LHE_COVERAGE            => 30,
  JF_SIZE                 => 100000000,
  KMER                    => "auto",
  KMER_COUNT_THRESHOLD    => 1,
  KMER_RELIABLE_THRESHOLD => 3,
  TRIM_PARAM              => 3,
  RELIABLE_Q_PARAM	  => 40,
  NUM_THREADS             => 2,
  NUM_CNS_THREADS         => 2,
  LIMIT_JUMP_COVERAGE     => 300,
  USE_LINKING_MATES       => 0,
  DO_HOMOPOLYMER_TRIM     => 0,
  CLOSE_GAPS              => 1,
  NO_MMAP                 => 1,
  STOP_AFTER_SR           => 0,
  CA_PARAMETERS           => "",
  SOAP_ASSEMBLY           => 0,

  # Data
  PE_INFO       => [],
  JUMP_INFO     => [],
  OTHER_INFO    => [],
  MOLECULO_INFO => [],
  PACBIO_INFO   => [],
  NANOPORE_INFO => [],
  REF_INFO => [],
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
PE= pe 180 20  /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq
JUMP= sh 3600 200  /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
#pacbio reads must be in a single fasta file! make sure you provide absolute path
PACBIO=/FULL_PATH/pacbio.fa
OTHER=/FULL_PATH/file.frg
END

PARAMETERS
#this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE = auto
#set this to 1 for all Illumina-only assemblies
#set this to 1 if you have less than 20x long reads (454, Sanger, Pacbio) and less than 50x CLONE coverage by Illumina, Sanger or 454 mate pairs
#otherwise keep at 0
USE_LINKING_MATES = 0
#this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms 
LIMIT_JUMP_COVERAGE = 300
#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
#set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS =  cgwErrorRate=0.15
#minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if Illumina coverage >100
KMER_COUNT_THRESHOLD = 1
#whether to attempt to close gaps in scaffolds with Illumina data
CLOSE_GAPS=1
#auto-detected number of cpus to use
NUM_THREADS = 16
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage
JF_SIZE = 200000000
#set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>5Gbp) genomes
SOAP_ASSEMBLY=0
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
  } elsif($key eq "REFERENCE") {
    can_read($param) or fail("invalid file for REFERENCE: '$param' $!", $.);
    push(@{$$res{REF_INFO}}, $param);
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
