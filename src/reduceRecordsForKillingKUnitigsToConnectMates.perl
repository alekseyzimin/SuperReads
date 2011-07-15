#!/usr/bin/perl
if ($#ARGV >= 0) {
    $outputDir = $ARGV[0]; }
$kmerLen = 31;
while ($line = <STDIN>) {
    chomp ($line);
    push (@lines, $line);
    $curLine = $#lines;
    @flds = split (" ", $line);
    next unless ($#flds > 0);
    if ($flds[6] < $flds[7]) {
	$begin[$curLine] = $flds[6];
	$end[$curLine] = $flds[7];
	$ori[$curLine] = "F"; }
    else {
	$begin[$curLine] = $flds[7];
	$end[$curLine] = $flds[6];
	$ori[$curLine] = "R"; }
    $begin[$curLine]--;
    ($ahg[$curLine]) = ($flds[3] =~ /^.([^,]+),/);
    $kUnitigNum[$curLine] = $flds[10];
    $readName[$curLine] = $flds[11];
}

for ($i=1; $i<=$#lines; $i++) {
    next unless ($begin[$i] =~ /\d/);
    if (($readName[$i] eq $readName[$i-1]) && ($end[$i-1]-$begin[$i] >= $kmerLen)) {
	$isChimeric{$readName[$i]} = 1; } }

$i=0;
# print $lines[$i];
# if ($isChimeric{$readName[$i]}) { print " CHIMERIC"; }
# print "\n";
for ($i=1; $i<=$#lines; $i++) {
#    print "i = ${i}: ";
#    print $lines[$i];
    if ($lines[$i] !~ /\S/) {
#	print "\n";
	next; }
    if ($isChimeric{$readName[$i]}) {
#	print " CHIMERIC\n";
	next; }
    if ($lines[$i-1] !~ /\S/) {
#	print "\n";
	next; }
    if ($readName[$i] eq $readName[$i-1]) {
	$diff = abs ($ahg[$i]-$ahg[$i-1]);
	if (! $gcd{$kUnitigNum[$i]}) {
	    $gcd{$kUnitigNum[$i]} = $diff; }
	else {
	    $gcd{$kUnitigNum[$i]} = &gcd ($gcd{$kUnitigNum[$i]}, $diff); }
#	print " kUnitigNum = $kUnitigNum[$i] $diff $gcd{$kUnitigNum[$i]}";
    }
#    print "\n";
}

if (! $outputDir) { print "Endgame:\n"; }
@chimericReads = keys %isChimeric;
@chimericReads = sort @chimericReads;
@repetitiveKUnitigs = keys %gcd;
if ($outputDir) {
    $outfile = "$outputDir/chimeric_read.txt";
    open (OUTFILE, ">$outfile"); }
if (! $outputDir) { print "Chimeric reads:\n"; }
for (@chimericReads) {
    if ($outfile) {
	print OUTFILE "$_\n"; }
    else {
	print "$_\n"; } }
if ($outfile) {
    close (OUTFILE); }
if (! $outputDir) {
    print "\nRepetitive k-unitigs and minimum spacing:\n"; }
if ($outputDir) {
    $outfile = "$outputDir/multiCopyKUnitigs.kUnitig.rptLength.txt";
    open (OUTFILE, ">$outfile"); }
for (@repetitiveKUnitigs) {
    $kUnitig = $_;
    if ($outfile) {
	print OUTFILE "$kUnitig $gcd{$kUnitig}\n"; }
    else {
	print "$kUnitig $gcd{$kUnitig}\n"; }
}
if ($outfile) {
    close (OUTFILE); }

sub gcd
{
    my ($a, $b) = @_;
    my ($swap);

    if ($a < $b) {
	$swap = $a;
	$a = $b;
	$b = $swap; }
    while ($b > 0) {
	$a %= $b;
	$swap = $a;
	$a = $b;
	$b = $swap; }
    return ($a);
}

