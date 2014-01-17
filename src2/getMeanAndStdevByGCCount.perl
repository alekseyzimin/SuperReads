#!/usr/bin/env perl
$intervalLen = 200;
&processArgs;
sub processArgs
{
    for ($i=0; $i<=$#ARGV; ++$i) {
	if ($ARGV[$i] eq "--interval-len") {
	    ++$i;
	    $intervalLen = $ARGV[$i];
	    next; }
    }
}

while ($line = <STDIN>) {
    chomp ($line);
    ($GCcount, $coverage, $countForPoint) = split (" ", $line);
    next if ($coverage == 0);
    next if ($countForPoint == 0);
    $GCpct = $GCcount / $intervalLen;
    $GCpercentage = int ($GCpct * 100);
    $diffLow = ($GCpct*100) - $GCpercentage;
    $diffHigh = $GCpercentage + 1 - $GCpct*100;
    if ($diffHigh > 0) { # Work on $GCpercentage
	$curGCpct = $GCpercentage;
	$oldAverage = $avg[$curGCpct];
	$oldWeightSum = $numLocsForGCByPct[$curGCpct];
	$weight = $diffHigh * $countForPoint;
	$numLocsForGCByPct[$curGCpct] += $weight;
	$avg[$curGCpct] = &calculateNewAverage ($oldAverage, $weight, $numLocsForGCByPct[$curGCpct], $coverage);
	$sumOfVariances[$curGCpct] = &calculateNewSumOfVariances ($sumOfVariances[$curGCpct], $weight, $coverage, $oldAverage, $avg[$curGCpct]);
    }
    if ($diffLow > 0) { # Work on $GCpercentage+1
	$curGCpct = $GCpercentage+1;
	$oldAverage = $avg[$curGCpct];
	$oldWeightSum = $numLocsForGCByPct[$curGCpct];
	$weight = $diffLow * $countForPoint;
	$numLocsForGCByPct[$curGCpct] += $weight;
	$avg[$curGCpct] = &calculateNewAverage ($oldAverage, $weight, $numLocsForGCByPct[$curGCpct], $coverage);
	$sumOfVariances[$curGCpct] = &calculateNewSumOfVariances ($sumOfVariances[$curGCpct], $weight, $coverage, $oldAverage, $avg[$curGCpct]);
    }
}

$maxAvg = 0;
for ($GCcount=40; $GCcount<=60; ++$GCcount) {
    if ($numLocsForGCByPct[$GCcount] == 0) {
	next; }
    $avgCoverage = $avg[$GCcount];
#    $avgCoverage = $sumReadsCoverageByPct[$GCcount] * 1.0 / $numLocsForGCByPct[$GCcount];
    if ($maxAvg < $avgCoverage) {
	$maxAvg = $avgCoverage; }
}

for ($GCcount=0; $GCcount<=$#numLocsForGCByPct; ++$GCcount) {
    if ($numLocsForGCByPct[$GCcount] == 0) {
        next; }
    $avgCoverage = $avg[$GCcount];
    $variance = $sumOfVariances[$GCcount] / $numLocsForGCByPct[$GCcount];
    $stdev = sqrt ($variance);
    $stdevToAvg = $stdev / $avgCoverage;
    $adjustmentFactor = $maxAvg / $avgCoverage;
    if ($adjustmentFactor < 1) {
	$adjustmentFactor = 1; }
    $numLocsForGCByPct[$GCcount] = int ($numLocsForGCByPct[$GCcount] + .5);
    if ($numLocsForGCByPct[$GCcount] < 1) {
	$numLocsForGCByPct[$GCcount] = 1; }
    $GCpct = $GCcount/100;
    if (0) {
	print "GCpct = $GCpct ";
	print "adjustmentFactor = $adjustmentFactor ";
	print "avgCoverage = $avgCoverage ";
	print "stdev = $stdev ";
	print "numLocsForGCByPct = $numLocsForGCByPct[$GCcount]\n";
    }
    print "$GCpct ";
    print "$adjustmentFactor ";
    print "$avgCoverage ";
    print "$stdev ";
    print "$numLocsForGCByPct[$GCcount]\n";
#    print "numLocsForGCByPct = $numLocsForGCByPct[$GCcount] avgCoverage = $avgCoverage stdev = $stdev stdev/avgCoverage = $stdevToAvg\n";
}

sub calculateNewAverage
{
    my ($priorAverage, $weight, $totalWeight, $sampleValue) = @_;
    my ($avg);
    $weight *= 1.0;
    $avg = $priorAverage + ($weight / $totalWeight) * ($sampleValue - $priorAverage);
    return ($avg);
}

sub calculateNewSumOfVariances
{
    my ($oldSumOfVariances, $weight, $coverage, $oldAverage, $newAverage) = @_;
    my ($newSumOfVariances);

    $newSumOfVariances = $oldSumOfVariances + $weight * ($coverage - $oldAverage) * ($coverage - $newAverage);
    return ($newSumOfVariances);
}

