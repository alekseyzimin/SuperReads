#!/usr/bin/env perl
my $kunitigNumberOffset=0;
my $maxLen=32767;
my $kmer=$ARGV[0];

while($line=<STDIN>){
  chomp($line);
  if(substr($line,0,1) eq ">"){
    @f=split(/\s+/,$line);
    my $kuNum=substr($f[0],1);
    my $len=substr($f[1],7);
    my $seq;
    $line=<STDIN>;
    chomp($line);
    $seq=$line;
    if($len>$maxLen){
      my $offset=0;
      my $outLen;
      while($len-$offset>=$kmer){
        $outLen=$len-$offset;
        $outLen=$maxLen if($outLen>$maxLen);
        print ">",$kuNum+$kunitigNumberOffset," length:$outLen\n",substr($seq,$offset,$outLen),"\n";
        $kunitigNumberOffset++;
        $offset+=($maxLen-$kmer+1);
      }
    }else{
      print ">",$kuNum+$kunitigNumberOffset," length:$len\n$seq\n";
    }
  }
}


