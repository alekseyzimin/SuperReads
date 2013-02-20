#!/usr/bin/env perl
# 1 argument: the name of the directory plus prefix for the Celera assembly (see the default value for an example)
# Output to STDOUT (except for the intermediate file)
use File::Basename;
$exeDir = dirname ($0);
if ($#ARGV >= 0) {
    $prefix = $ARGV[0]; }
else {
    $prefix = "/genome2/raid/tri/rhodobacterAssembly_alekseyFromGenome10/CA/9-terminator/rh"; }

$partPrefix = basename ($prefix);
$intermFile = "${partPrefix}.utgWithTypes.txt";
$asmFile = "${prefix}.asm";
$utglenFile = "${prefix}.posmap.utglen";
$frgutgFile = "${prefix}.posmap.frgutg";
$utgctgFile = "${prefix}.posmap.utgctg";
$frgctgFile = "${prefix}.posmap.frgctg";
$fail = 0;
if (! -e $asmFile) { print STDERR "$asmFile doesn't exist!\n"; $fail = 1; }
if (! -e $utglenFile) { print STDERR "$utglenFile doesn't exist!\n"; $fail = 1; }
if (! -e $frgutgFile) { print STDERR "$frgutgFile doesn't exist!\n"; $fail = 1; }
if (! -e $utgctgFile) { print STDERR "$utgctgFile doesn't exist!\n"; $fail = 1; }
if (! -e $frgctgFile) { print STDERR "$frgctgFile doesn't exist!\n"; $fail = 1; }
if ($fail) { print STDERR "Bye!\n"; exit (1); }

$cmd = "$exeDir/getUnitigTypeFromAsmFile.perl $asmFile > $intermFile";
print STDERR "$cmd\n"; system ($cmd);

$cmd = "$exeDir/addSurrogatesToFrgCtgFile $intermFile $utglenFile $frgutgFile $utgctgFile $frgctgFile";
print STDERR "$cmd\n"; system ($cmd);


