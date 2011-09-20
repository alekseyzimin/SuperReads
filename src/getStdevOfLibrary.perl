#!/usr/bin/perl
# Input a list of values in STDIN
# Use a -stdev # flag to eliminate mates more than # std devs from the mean
#   (iteratively) (default is 6)
# Outputs the mean and std dev to STDOUT
&processArgs;
while ($val = <STDIN>) {
    chomp ($val);
    push (@vals, $val);
}
&calculateFirstMeanAndStdev;
print "Initial: $numVals $mean $stdev\n";

while (1) {
    $changesWereMade = 0;
    &calculateSubsequentMeanAndStdev;
    last unless ($changesWereMade);
    print "$numVals $mean $stdev\n";
}
&calculateFinalStdev;
print "Final: $numVals $mean $stdev\n";

# The following is necessary because it's more stable than the loop calc
sub calculateFinalStdev
{
    my ($sumDiffSqVals, $i, $val, $tval);

    $sumDiffSqVals = 0;
    for ($i=0; $i<=$#vals; $i++) {
	next if ($isExcluded[$i]);
	$val = $vals[$i];
	$tval = $val - $mean;
	$sumDiffSqVals += ($tval * $tval);
    }
    $stdev = sqrt ($sumDiffSqVals / $numVals);
}

sub calculateSubsequentMeanAndStdev
{
    my ($sumVals, $sumSqVals, $i, $val, $isExcluded);

    $sumVals = $sumSqVals = $numVals = 0;
    for ($i=0; $i<=$#vals; $i++) {
	$val = $vals[$i];
	$isExcluded = 0;
	if ($val < $mean - $stdev * $numStdev) { $isExcluded = 1; }
	elsif ($val > $mean + $stdev * $numStdev) { $isExcluded = 1; }
	if ($isExcluded != $isExcluded[$i]) { ++ $changesWereMade; }
	$isExcluded[$i] = $isExcluded;
	next if ($isExcluded);
	$sumVals += $val;
	$sumSqVals += ($val * $val);
	$numVals++;
    }
    $mean = $sumVals / $numVals;
    $stdev = sqrt ($sumSqVals / $numVals - $mean * $mean);
}

sub calculateFirstMeanAndStdev
{
    my ($i, $val, $sumVals, $sumSqVals);

    for ($i=0; $i<=$#vals; $i++) {
	$val = $vals[$i];
	$sumVals += $val;
	$sumSqVals += ($val * $val);
	$numVals++;
    }
    $mean = $sumVals / $numVals;
    $stdev = sqrt ($sumSqVals / $numVals - $mean * $mean);
}

sub processArgs
{
    my ($i, $arg);

    $numStdev = 6;
    for ($i=0; $i<=$#ARGV; $i++) {
	$arg = $ARGV[$i];
	if ($arg eq "-stdev") {
	    ++$i; $arg = $ARGV[$i];
	    $numStdev = $arg;
	}
    }
}
