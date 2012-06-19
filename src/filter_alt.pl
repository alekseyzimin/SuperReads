#!/usr/bin/perl

$last_prefix="";
$last_readnum=-2;
$last_superread="";
$last_pos;

while($line=<STDIN>){
  chomp($line);
  @F=split(" ",$line);
  $prefix=substr($F[0],0,2);
  $readnum=int(substr($F[0],2));
  if($readnum%2==1){
    if($prefix eq $last_prefix && ($readnum-1)==$last_readnum){
      if($F[1] eq $last_superread){
        if(abs($last_pos-$F[2])<700){
          print "$prefix$last_readnum\n$prefix$readnum\n";
        }
      }else{
        @k_u=split(/_/,$last_superread);
        %h=();
        foreach $k(@k_u){
          $h{substr($k,0,length($k)-1)}=1;
        }
        $flag=0;
        @k_u=split(/_/,$F[1]);
        foreach $k(@k_u){
          if(defined($h{substr($k,0,length($k)-1)})){
            $flag=1;
            last;
          }
        }
        if($flag==1){
          print "$prefix$last_readnum\n$prefix$readnum\n";
        }
      }
    }
  }else{
    $last_prefix=$prefix;
    $last_readnum=$readnum;
    $last_superread=$F[1];
    $last_pos=$F[2];
  }
} 

