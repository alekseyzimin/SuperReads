#!/usr/bin/env perl
#
# We output "TMN" (too many nodes) immediately when we see 1 record because we
# only need 1 read of the pair to send it to the next pass. (If one read has
# errors then it won't be covered by k-unitigs and thus won't have a super-read.
# To avoid this possibility we only require 1 of the 2 reads to be covered.
$inputReadPositionsFile = $ARGV[0];
$outputJoinedReadsFile = $ARGV[1]; # We will append to this
$outputListOfReadsForNextPassFile = $ARGV[2]; # We will append to this
$minDiffToAllow = 0;
$maxDiffToAllow = 800;
open (FILE, $inputReadPositionsFile);
open (OUTFILE1, ">>$outputJoinedReadsFile");
open (OUTFILE2, ">>$outputListOfReadsForNextPassFile");
while ($line = <FILE>) {
    chomp ($line);
    ($rd, $super, $offset, $ori, $type) = split (" ", $line);
    if ($rd =~ /[02468]$/) {
	$evenRd = $rd;
	$typeHold = $type;
	$oriHold = $ori;
	$offsetHold = $offset;
	$superHold = $super;
	if ($type eq "TMN") {
	    $outStr = &createOutputTMNString ($rd);
	    print OUTFILE2 $outStr; }
	next; }
    # Only output when we hit the odd numbered read of the pair
    ($prefix, $last) = ($rd =~ /^(.+)(.)$/);
    --$last;
    $tempReadName = $prefix . $last;
    if ($tempReadName ne $evenRd) {
	if ($type eq "TMN") {
	    $outStr = &createOutputTMNString ($rd);
	    print OUTFILE2 $outStr; }
	next; }
    # If we get here ($tempReadName == $evenRd)
    next if (($type eq "MS") || ($type eq "TMN")); # TMN already output
    if (($type eq "A") || ($type eq "J")) {
	print OUTFILE1 "$evenRd\n$rd\n";
	next; }
    # Now we only have the case of 'SU' (same unitig)
    if ($super ne $superHold) {
	print OUTFILE2 "$evenRd\n$rd\n";
	next; }
    # Now 'SU' and both super-reads are the same (i.e. one unitig)
    next if ($ori eq $oriHold);
    if ($ori eq "F") {
	$diff = $offsetHold - $offset; }
    else {
	$diff = $offset - $offsetHold; }
    # $diff is the oriented separation; 'R' end - 'F' end
    if (($diff >= $minDiffToAllow) && ($diff <= $maxDiffToAllow)) {
	print OUTFILE1 "$evenRd\n$rd\n"; }
    next;
}
close (FILE);
close (OUTFILE1);
close (OUTFILE2);

sub createOutputTMNString
{
    my ($inputRd) = @_;
    my ($rd1, $rd2, $pref, $last, $outstr);

    ($pref, $last) = ($inputRd =~ /^(.+)(.)$/);
    if ($last % 2 == 0) {
	$rd1 = $inputRd;
	++$last;
	$rd2 = $pref . $last; }
    else {
	$rd2 = $inputRd;
	--$last;
	$rd1 = $pref . $last; }
    $outstr = "$rd1\n$rd2\n";
    return ($outstr);
}

