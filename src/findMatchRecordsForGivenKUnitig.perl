#!/usr/bin/perl
$kUni = $ARGV[-1];
$cmd = "grep -10 \" $kUni \" Illumina_data/work1/arrangedCoordsResultsByRead.txt |";
open (FILE, $cmd);
while ($line = <FILE>) {
    chomp ($line);
    if ($line =~ /^readNum/) {
	($readNum) = ($line =~ /^readNum = (\d+)\D/);
	$lines{$readNum} .= "$line\n";
	next; }
    if ($line !~ /\S/) {
	$lines{$readNum} .= "\n";
	next; }
    $lines{$readNum} .= "$line\n";
    ++$lineCount{$readNum};
    @flds = split (" ", $line);
    if ($flds[6] == $kUni) {
	$isNeeded{$readNum} = 1;
	if ($flds[2] < $flds[3]) { $minOffset = $flds[2]; $maxOffset = $flds[3]; }
	else { $minOffset = $flds[3]; $maxOffset = $flds[2]; }
	if (($minOffset == 1) && ($maxOffset == $flds[5])) {
	    $isCovered{$readNum} = 1; } }
}

@readNames = keys %lines;
@readNames = sort byNum @readNames;
for (@readNames) {
    $readName = $_;
    next unless ($isNeeded{$readName});
    print "readName = \"$readName\"; ";
    if ($isNeeded{$readName}) { print "isNeeded "; }
    print "lineCount = $lineCount{$readName}; ";
    print "isCovered = $isCovered{$readName}\n";
    print $lines{$readName};
}

sub byNum
{
    return ($a <=> $b);
}

