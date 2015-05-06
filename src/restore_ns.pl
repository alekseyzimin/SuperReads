#!/usr/bin/env perl
my $contig_file=$ARGV[0];
my $scaff_pos_file=$ARGV[1];

my $ctg=-1;
my $seq="";
#here we note positions of all N sections in the contigs
open(FILE,$contig_file);
while($line=<FILE>){
	chomp($line);
	if($line=~/^>/){
		if($ctg>-1){
		$contigs[$ctg]=$seq;
		}
	$seq="";
	($ctg)=split(/\s+/,substr($line,1));
	}else{
		$seq.=$line;
	}
}
if($ctg>-1){
	$contigs[$ctg]=$seq;
}

#now we read the contig positions file
my $scf="";
my $scfseq="";
my $last_pos=0;
open(FILE,$scaff_pos_file);
while($line=<FILE>){
        chomp($line);  
        if($line=~/^>/){   
		print $scfseq,"\n" if(not($scfseq eq ""));
		print "$line 0.0\n";
		$scfseq="";
		$last_pos=0;
	
	}else{
		@f=split(/\s+/,$line);
		if($f[1]>0){
		my $gap=$f[1]-$last_pos;
		#print "gap=$gap\n";
		$scfseq.="N"x$gap;
		}
		$last_pos=$f[1]+$f[3];	
		if($f[2] eq "+"){
		$scfseq.=$contigs[$f[0]];
		}else{
		my $ctg=$contigs[$f[0]];
		$ctg=~tr/ACGTNacgtn/TGCANtgcan/;
		$scfseq.=reverse($ctg);
		}
		#print $line," ",length($contigs[$f[0]]), " ",length($scfseq),"\n";
		$contigs[$f[0]]="";
	}
}
print $scfseq,"\n" if(not($scfseq eq ""));

my $i=0;
foreach $c(@contigs){
print ">C$i  0.0\n$c\n" if(not($c eq ""));
$i++;
}



