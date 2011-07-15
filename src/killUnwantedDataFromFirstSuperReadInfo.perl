#!/usr/bin/perl
# Input on STDIN, output on STDOUT
# The working directory is expected as an argument
if ($#ARGV >= 0) {
    $workingDir = $ARGV[0]; }
else {
    $workingDir = "."; }
$readKillFile = "$workingDir/chimeric_read.txt";
$superReadKillFile1 = "$workingDir/createFastaSuperReadSequences.errors.txt";
$superReadKillFile2 = "$workingDir/superReadsWKUnitigOvlGeKmerLen.numReads.superReadName.txt";
$infile5 = "";
open (FILE, $readKillFile);
while ($readName = <FILE>) {
    chomp ($readName);
    $isChimeric{$readName} = 1; }
close (FILE);
open (FILE, $superReadKillFile1);
while ($line = <FILE>) {
    chomp ($line);
    next unless ($line =~ /^\d/);
    @flds = split (" ", $line);
    $superRead = $flds[0];
    $eliminateSuperRead{$superRead} = 1; }
close (FILE);

open (FILE, $superReadKillFile2);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line); 
    $superRead = $flds[1];
    $eliminateSuperRead{$superRead} = 1; }
close (FILE);

while ($line = <STDIN>) {
    chomp ($line);
    @flds = split (" ", $line);
    $readName = $flds[5];
    $line2 = <STDIN>;
    chomp ($line2);
    @flds = split (" ", $line2);
    $superRead = $flds[0] . $flds[1];
    for ($i=2; $i<$#flds; $i+=3) {
	$superRead .= ("_" . $flds[$i] . "_" . $flds[$i+1] . $flds[$i+2]); }
    next if ($eliminateSuperRead{$superRead});
    next if ($isChimeric{$readName});
    print "$line\n$line2\n"; }
