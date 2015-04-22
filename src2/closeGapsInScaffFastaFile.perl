#!/usr/bin/env perl
# This program is used to close gaps in a scaffold sequence file.
# Example invocation:
# closeGapsInScaffFastaFile.perl --scaffold-fasta-file 9-terminator/genome.scf.fa --min-kmer-len 17 --max-ker-len 31 --work-directory . --num-threads 16 --output-directory outputDir --reads-file pe.cor.fa --reads-file sj.cor.fa
#
# There are no args, only flags (mostly) with arguments. They are as follows:
#
# Required flags:
# --scaffold-fasta-file filename : file containing the scaffold sequences
#
# All the other switches are the switches understood by closeGapsLocally.

use strict;
use warnings;
use File::Basename;
use POSIX;

my $exeDir = dirname ($0);
my $scaffoldFastaFile;

sub processArgs {
  my @nargv;

  for(my $i = 0; $i < @ARGV; $i++) {
    if($ARGV[$i] eq "--scaffold-fasta-file") {
      $scaffoldFastaFile = $ARGV[$i + 1];
      $i++;
    } else {
      push(@nargv, $ARGV[$i]);
    }
  }
  return @nargv;
}

sub system_error {
  my ($msg) = @_;
  print(STDERR $msg, ": ") if $msg;
  my $e = 1;
  my $st = ${^CHILD_ERROR_NATIVE};
  if(WIFEXITED($st)) {
    $e = WEXITSTATUS($st);
    printf(STDERR "Exit status is %d\n", $e);
  } elsif(WIFSIGNALED($st)) {
    my $s = WTERMSIG($st);
    print(STDERR "Killed by signal $s\n");
  } else {
    print(STDERR "Unknown error\n");
  }
  exit($e);
}

sub reportUsage {
  open (FILE, $0);
  my $line = <FILE>;
  while ($line = <FILE>) {
    last unless ($line =~ /^\#/);
    chomp ($line);
    my ($line2) = ($line =~ /^..(.*)$/);
    print(STDERR $line2 || "", "\n");
  }
  close (FILE);
  exit (1);
}

my @nargv = processArgs;
reportUsage if(!defined($scaffoldFastaFile));

# Also generates genome.posmap.ctgscf
# Also generates genome.asm
my $cmd = "$exeDir/splitFileAtNs '$scaffoldFastaFile' > genome.ctg.fasta";
print "$cmd\n";
system ($cmd) == 0 or
   system_error("Splitting scaffold file at Ns failed");

my @cmd = ("$exeDir/closeGapsLocally.perl",  @nargv, "--Celera-terminator-directory", ".");
print(join(" ", @cmd), "\n");
system (@cmd) == 0 or
   system_error("closeGapsLocally failed");
