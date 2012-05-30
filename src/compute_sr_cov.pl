#!/usr/bin/env perl
# SuperRead pipeline
# Copyright (C) 2012  Genome group at University of Maryland.
# 
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#
my $total_rho=0;
my $total_count=0;
my $uidfile=$ARGV[0];
my $countsfile=$ARGV[1];
my $readlen=$ARGV[2];
my $renamed_sr=$ARGV[3];

my %r_sr;
open(FILE,$renamed_sr);
while($line=<FILE>){
  chomp($line);
  @f=split(/\s+/,$line);
  $r_sr{$f[0]}=$f[1];
}
close(FILE);


my @uid;
open(FILE,$uidfile);
while($line=<FILE>){
  chomp($line);
  push(@uid,$line);
}
close(FILE);

my %counts;
open(FILE,$countsfile);
while($line=<FILE>){
  @f=split(/\s+/,$line);
  if(defined($counts{$f[1]})){
    $counts{$f[1]}++;
  }
  else{
    $counts{$f[1]}=1;
  }
}
close(FILE);

#now if we need to find out the count by iid, it is here $counts{$uid[$iid]}

my $total_rho=0;
my $total_count=0;
my $utg=-1;
while($line=<STDIN>){
  if($line=~/^unitig/){
    chomp($line);
    @l=split(/\s+/,$line);
    if($utg ==-1){
      $utg=$l[1];
      $c=0;
      $r=-1;
    }
    else{
      if($c>0){
        $count{$utg}=$c;
        $rho{$utg}=$r-$readlen;
        $rho{$utg}=1 if($rho{$utg}<0);
        if($r>2000){
          $total_rho+=$rho{$utg};
          $total_count+=$c;
        }
      }
      $c=0;
      $r=-1;
      $utg=$l[1];
    }
  }
  elsif($line =~ /^cns/){
    chomp($line);
    @s=split(//,substr($line,4));
    $cg=0;
    $base_count=1;
    foreach $v(@s){
      $cg++ if($v eq "G" || $v eq "C");
      $base_count++;
    } 
    $cg_content{$utg}=$cg/$base_count;
  }
  elsif($line =~ /^FRG/){
    @l=split(/\s+/,$line);
    @f=split(/\:/,$uid[$l[4]]);
    
    if(defined($counts{$f[0]})){
      $r=$l[13] if($l[13]>$r);
      $r=$l[14] if($l[14]>$r);
      
      $c+=$counts{$f[0]};
      $counts{$f[0]}=0;
    }
    elsif(defined($counts{$r_sr{$f[0]}})){
      $r=$l[13] if($l[13]>$r);
      $r=$l[14] if($l[14]>$r);
      
      $c+=$counts{$r_sr{$f[0]}};
      $counts{$r_sr{$f[0]}}=0;
    }
  }
}
if($c>0)
{
  $count{$utg}=$c;
  $rho{$utg}=($r-$readlen);
  $rho{$utg}=1 if($rho{$utg}<0);
  $total_rho+=$rho{$utg};
  $total_count+=$c;
}

my $global_arrival_rate=$total_count/$total_rho;
print STDERR "total_rho= $total_rho total_count= $total_count GAR= $global_arrival_rate\n";
foreach $v(keys %count)
{
  print STDERR "$v $cg_content{$v}\n";
  my $astat=($rho{$v}*$global_arrival_rate)-(0.6931471805599453094*$count{$v});
#if($rho{$v}>2000&&$astat<0)
#{
#print "unitig_coverage_stat $v 10\n";
#}
#else
#{
  print "unitig_coverage_stat $v $astat\n";
#}
}

