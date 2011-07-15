#!/usr/bin/perl
# Cat the superReadGroupsForEachReadWhichHasAGroup.txt file in and this
# outputs the superReadGroups with the read info for them to STDOUT
$isFirstLine = 1;
while ($line = <STDIN>) {
    chomp ($line);
    if ($isFirstLine) {
	@flds = split (" ", $line);
	if ($flds[0] eq "readNum") {
	    $hasReadNum = 1; }
	else {
	    $hasReadNum = 0; }
	$isFirstLine = 0; }
    if ($hasReadNum) {
	($readNum, $readName) = ($line =~ /^readNum = (\d+); readName = (\S+)\s/); }
    else {
	($readName) = ($line =~ /^readName = (\S+)\s/);
	$readNum = $readName; }
    $superReadGroup = <STDIN>;
    if (! $readsInSuperReadGroup{$superReadGroup}) {
	push (@superReadGroups, $superReadGroup); }
    else {
	$readsInSuperReadGroup{$superReadGroup} .= " ; "; }
    $readsInSuperReadGroup{$superReadGroup} .= "$readNum $readName";
}

for (@superReadGroups) {
    $superReadGroup = $_;
    $readsInSuperReadGroup = $readsInSuperReadGroup{$superReadGroup};
    chomp ($superReadGroup);
    print "$superReadGroup : $readsInSuperReadGroup\n"; }
