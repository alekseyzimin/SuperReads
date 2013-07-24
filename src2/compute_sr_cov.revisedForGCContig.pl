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
my $uidfile = $ARGV[0];
my $countsfile = $ARGV[1];
my $readlen = $ARGV[2];
my $srFRGfile = $ARGV[3];
my $adjstForGCContentfile = "";
if ($#ARGV >= 4) {
    $adjstForGCContentfile = $ARGV[4]; }

my @uid;
open (FILE, $uidfile);
while ($line = <FILE>) {
    chomp ($line);
    push (@uid, $line);
}
close(FILE);

my @counts;
open (FILE, $countsfile);
while ($line = <FILE>) {
    @flds = split (" ", $line);
    $counts[$flds[1]]++; }
close(FILE);

# The following is for file with 3 fields per line
# unitig# ratioAdjustment GCproportion      (the last in [0,1])
if ($adjstForGCContentfile) {
    open (FILE, $adjstForGCContentfile);
    while ($line = <FILE>) {
	chomp ($line);
	@flds = split (" ", $line);
#	$factor[$flds[0]] = $flds[1];
	$factorBasedOnWholeUnitigGC[$flds[0]] = $flds[1];
	$gc_content[$flds[0]] = $flds[2];
    }
    close (FILE); }
my @parts;
open (FILE, $srFRGfile);
while ($line = <FILE>) {
    if ($line =~ /^acc:/) {
   @flds = split (/\:/, $line);
   $parts[$flds[1]]++; }
}
close(FILE);

#now if we need to find out the count by iid, it is here $counts{$uid[$iid]}

my $total_rho = 0;
my $total_count = 0;
my $utg = -1;
while ($line = <STDIN>) {
    if ($line =~ /^unitig/) {
	chomp ($line);
	@flds = split (" ", $line);
	if ($utg == -1) {
	    $utg = $flds[1];
	    $num_superread_reads_in_unitig = 0;
	    $max_offset_for_superread_in_unitig = -1; }
	else {
	    if ($num_superread_reads_in_unitig > 0) {
		if ($adjstForGCContentfile =~ /\S/) {
		    $num_superread_reads_in_unitig *= $factorBasedOnWholeUnitigGC[$utg]; }
		$count{$utg} = $num_superread_reads_in_unitig;
		$rho{$utg} = $max_offset_for_superread_in_unitig - $readlen;
		$rho{$utg} = 1 if ($rho{$utg} < 0);
		if ($max_offset_for_superread_in_unitig > 2000) {
		    $total_rho += $rho{$utg};
		    $total_count += $num_superread_reads_in_unitig; } }
	    $num_superread_reads_in_unitig = 0;
	    $max_offset_for_superread_in_unitig = -1;
	    $utg = $flds[1]; }
    }
#    elsif ($line =~ /^cns/) {
#    }
    elsif ($line =~ /^FRG/) {
	@flds = split (" ", $line);
	@flds2 = split (/\:/, $uid[$flds[4]]);
    
	next if ($flds2[1] !~ /^super-read/);
	if (($counts[$flds2[0]] > 0) && ($parts[$flds2[0]] > 0)) {
	    $max_offset_for_superread_in_unitig = $flds[13] if ($flds[13] > $max_offset_for_superread_in_unitig);
	    $max_offset_for_superread_in_unitig = $flds[14] if ($flds[14] > $max_offset_for_superread_in_unitig);
      
	    $num_superread_reads_in_unitig += $counts[$flds2[0]] / $parts[$flds2[0]]; } }
}

if ($num_superread_reads_in_unitig > 0) {
    if ($adjstForGCContentfile =~ /\S/) {
	$num_superread_reads_in_unitig *= $factorBasedOnWholeUnitigGC[$utg]; }
    $count{$utg} = $num_superread_reads_in_unitig;
    $rho{$utg} = ($max_offset_for_superread_in_unitig - $readlen);
    $rho{$utg} = 1 if ($rho{$utg} < 0);
    $total_rho += $rho{$utg};
    $total_count += $num_superread_reads_in_unitig;
}

my $global_arrival_rate = $total_count / $total_rho;
print STDERR "total_rho= $total_rho total_count= $total_count GAR= $global_arrival_rate\n";
foreach $v (keys %count) {
#    print STDERR "$v $gc_content[$v]\n";
    my $astat = ($rho{$v}*$global_arrival_rate) - (0.6931471805599453094*$count{$v});
    print "unitig_coverage_stat $v $astat\n";
}

