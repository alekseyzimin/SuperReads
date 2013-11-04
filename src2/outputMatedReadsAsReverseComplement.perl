#! /usr/bin/env perl
# Input a fasta file of reads where mated reads are together
# Output a fasta file of the mates reads with sequence reverse-complemented
# Input via the file in the 1 argument; output to STDOUT
$inputFilename = $ARGV[0];
if ($#ARGV >= 1) {
    if ($ARGV[1] eq "--reverse-complement") {
	$reverseComplement = 1; }
}
open (FILE, $inputFilename);
$fauxNumber = 0;
while ($line = <FILE>) {
    chomp ($line);
    if ($line =~ /^>/) {
	@flds = split (" ", $line);
	($readName2, $readName2suffix) = ($flds[0] =~ /^.(.*)(.)$/);
	if ($outputTheSequence eq "True") {
	    if ($reverseComplement) {
		$readSequence = &reverse_complement ($readSequence); }
	    print "$readSequence\n";
	    $outputTheSequence = "False"; }
	elsif (($readName2 eq $readName1) && ($readName2suffix == $readName1suffix+1) && ($readName1suffix % 2 == 0)) {
	    if ($reverseComplement) {
		$readSequence = &reverse_complement ($readSequence); }
	    $fauxNumberPlus1 = $fauxNumber + 1;
	    print ">cc$fauxNumber ${readName1}$readName1suffix r\n$readSequence\n>cc$fauxNumberPlus1 ${readName2}$readName2suffix r\n";
	    $fauxNumber += 2;
	    $outputTheSequence = "True"; }
	else {
	    $outputTheSequence = "False"; }
	$readSequence = "";
	$readName1 = $readName2;
	$readName1suffix = $readName2suffix;
	next; }
    # If we get here we're collecting sequence
    $readSequence .= $line;
}
close (FILE);

if ($outputTheSequence eq "True") {
    if ($reverseComplement) {
	$readSequence = &reverse_complement ($readSequence); }
    print "$readSequence\n"; }

exit (0);

sub reverse_complement
{
    my ($string) = @_;

    $string =~ tr/ACGTacgt/TGCAtgca/;
    $string = reverse ($string);
    return ($string);
}

