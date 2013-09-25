#!/usr/bin/env perl

my $library = $ARGV[0];
my $infile1 = $ARGV[1];
my $infile2 = $ARGV[2];

if($infile1 eq $infile2 || $infile2 eq ""){
    open(INFILE1,"zcat -cf $infile1 | "); 
    $readnumber=0;
    $seq1=""; 
    $qlt1="";
    while(defined($line1=<INFILE1>)){
	if($line1=~ /^@/){ 
	    if(defined($line1=<INFILE1>)){
		chomp($line1);
		$seq1=$line1;
	    }else{
		last;
	    }
	    if(defined($line1=<INFILE1>)){
		if($line1=~ /^+/){
		}else{
		    last;
		}
	    }
	    if(defined($line1=<INFILE1>)){
		chomp($line1);
		$qlt1=$line1;
	    }else{
		last;
	    }
	    if($seq1 !~ /[^ACGTN]/){
		print "\@","$library$readnumber\n$seq1\n+\n$qlt1\n";
	    }
	    $readnumber+=2;
	}
    }
}else{
    open(INFILE1,"zcat -cf $infile1 | ");open(INFILE2,"zcat -cf $infile2 | "); 
    $readnumber=0; 
    $seq1=""; 
    $seq2=""; 
    $qlt1=""; 
    $qlt2="";
    while(defined($line1=<INFILE1>) && defined($line2=<INFILE2>)){
	if($line1=~ /^@/ && $line2=~ /^@/){ 
	    if(defined($line1=<INFILE1>) && defined($line2=<INFILE2>)){
		chomp($line1);
		$seq1=$line1;
		chomp($line2);
		$seq2=$line2;
	    }else{
		last;
	    }
	    if(defined($line1=<INFILE1>) && defined($line2=<INFILE2>)){
		if($line1=~ /^+/ && $line2=~ /^+/){
		}else{
		    last;
		}
	    }
	    if(defined($line1=<INFILE1>) && defined($line2=<INFILE2>)){
		chomp($line1);
		$qlt1=$line1;
		chomp($line2);
		$qlt2=$line2;
	    }else{
		last;
	    }
	    if($seq1 !~ /[^ACGTN]/){
		print "\@","$library$readnumber\n$seq1\n+\n$qlt1\n";
	    }
	    $readnumber++;
	    if($seq2 !~ /[^ACGTN]/){
		print "\@","$library$readnumber\n$seq2\n+\n$qlt2\n";
	    }
	    $readnumber++;
	}
    }
}
