#!/usr/bin/perl
$shootingMethodSuperReadGroupFile = $ARGV[0];
$origSuperReadGroupFile = $ARGV[1];

open (FILE, $shootingMethodSuperReadGroupFile);
while ($line = <FILE>) {
    print $line;
    ($readName) = ($line =~ /^readName = (\d+)/);
    $alreadyOutput{$readName} = 1;
    $line = <FILE>;
    print $line;
}
close (FILE);

open (FILE, $origSuperReadGroupFile);
while ($line = <FILE>) {
    ($readName) = ($line =~ /^readName = (\d+)/);
    print $line unless ($alreadyOutput{$readName});
    $line = <FILE>;
    print $line unless ($alreadyOutput{$readName});
}
close (FILE);

