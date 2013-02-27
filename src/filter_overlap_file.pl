#!/usr/bin/env perl


my $gkpStore=$ARGV[0];
my $repeatKmerFile=$ARGV[1];
my $kmer=$ARGV[2];


open(FILE,$repeatKmerFile);
while($line=<FILE>){
    next if($line =~ /^>/);
    chomp($line);
    $repKMers{$line}=1;
    $repKMers{reverse_complement($line)}=1;
}
close(FILE);

my @frags=();
open(FILE,"gatekeeper -dumpfastaseq $gkpStore  | ");
while($line=<FILE>){
    if($line =~ /^>/){
        chomp($line);
        my @f=split(' ',$line);
	my @ff=split(/,/,$f[0]);
        $iid=$ff[1];
    }else{
        chomp($line);
        $frags[$iid].=$line;
    }
}
close(FILE);

while($line=<STDIN>){
    chomp($line);
    $line =~ s/^\s+//;
    my @f=split(' ',$line);
    next if($f[0]>$f[1]);
    my %kmers1=();
    my %kmers2=();
    my $read1=$frags[$f[0]];
    my $read2;
    if($f[2] eq "I"){
    $read2=reverse_complement($frags[$f[1]]);
    }else{
    $read2=$frags[$f[1]];
    }
    if($f[3]>=0){
	$start1=$f[3];
	$start2=0;
    }else{
	$start1=0;
	$start2=-$f[3];
    }
    if($f[4]>=0){
	$end1=length($read1);
	$end2=length($read2)-$f[4];
    }else{
	$end1=length($read1)+$f[4];
	$end2=length($read2);
    }
    #print "$line\ncoordinates $start1 $end1 $start2 $end2 ",length($read1)," ",length($read2),"\n";
    for($i=$start1;$i<=$end1-$kmer;$i++){
	$kmers1{substr($read1,$i,$kmer)}=1;
    }
    for($i=$start2;$i<=$end2-$kmer;$i++){
	my $local_kmer=substr($read2,$i,$kmer);
	$kmers2{$local_kmer}=1 if(defined($kmers1{$local_kmer}));
    }
    foreach $k(keys %kmers2){
	if(not(defined($repKMers{$k}))){
	    print $line,"\n";
	    last;
	}
    }
} 

sub reverse_complement{
    my $str=$_[0];
    $str =~ tr/acgtACGTNn/tgcaTGCANn/;
    $str = reverse ($str);
    return ($str);
}


