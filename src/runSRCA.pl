#! /usr/bin/env perl

use FindBin qw($Bin);

print(STDERR "The runSRCA.pl is deprecated. Run 'masurca' instead.\n");
exec($Bin . "/masurca", @ARGV);
