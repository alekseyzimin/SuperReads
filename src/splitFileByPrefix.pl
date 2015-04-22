#!/usr/bin/env perl

my $index=0;
my @fh=();
my %findex=();

foreach my $v (@ARGV) {
            local *OUTFILE;
            open(OUTFILE,">$v.cor.clean.fa");
            push(@fh,*OUTFILE);
            $findex{$v}=$index;
            ++$index;
}

my $current_index=0;
while($line=<STDIN>){
$current_index=$findex{substr($line,1,2)} if($line=~ /^>/);
print {$fh[$current_index]} $line;
$line=<STDIN>;
print {$fh[$current_index]} $line;
}
foreach my $v (@fh) {
close($v);
}

