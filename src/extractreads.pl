#!/usr/bin/perl
$infile1 = $ARGV[0]; #file with read names
$infile2 = $ARGV[1]; #fasta file
$fieldnum = $ARGV[2]-1; #field number

open (FILE1, $infile1);
open (FILE2, $infile2);
my %readnames;
while ($line = <FILE1>)
  {
    chomp($line);
    $readnames{$line}=1;
  }
close(FILE1);

my $sequence="";
my $readname="";
open (OUTFILE, ">$outfile");
while ($line = <FILE2>)
{
    if ($line =~ /^>/)
    {
	if(defined $readnames{$readname})
	{
	print ">$readname\n$sequence\n";
	}
	chomp($line);
        @f=split(/\s+/,substr($line,1));
        $readname=$f[$fieldnum];
        $sequence="";
    }
    else
    {
     chomp($line);
     $sequence.=$line;
    }
}
        if(defined $readnames{$readname})
        {
        print ">$readname\n$sequence\n";
        }

close (FILE2);

