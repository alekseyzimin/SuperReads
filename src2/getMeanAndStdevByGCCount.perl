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

$minRatioStdErrToMean = .02;
# $GCcount is actually percent now
for ($GCcount=0; $GCcount<=$#numLocsForGCByPct; ++$GCcount) {
    if ($numLocsForGCByPct[$GCcount] == 0) {
	$isBadValue[$GCcount] = 1;
        next; }
    $avgCoverage = $avg[$GCcount];
    $variance = $sumOfVariances[$GCcount] / $numLocsForGCByPct[$GCcount];
    $stdev = sqrt ($variance);
    $stdev /= sqrt ($numLocsForGCByPct[$GCcount]);
    $stdevToAvg = $stdev / $avgCoverage;
    if (($stdevToAvg >= $minRatioStdErrToMean) || ($numLocsForGCByPct[$GCcount] < 4)) {
	$isBadValue[$GCcount] = 1; }
    $adjustmentFactor = $maxAvg / $avgCoverage;
    if ($adjustmentFactor < 1) {
	$adjustmentFactor = 1; }
    $adjustmentFactor[$GCcount] = $adjustmentFactor;
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
    $outputLine[$GCcount] = "$adjustmentFactor $avgCoverage $stdev $stdevToAvg $numLocsForGCByPct[$GCcount]";
    if (0) {
	print "$GCpct ";
	print "$adjustmentFactor ";
	print "$avgCoverage ";
	print "$stdev ";
	print "$stdevToAvg ";
	print "$numLocsForGCByPct[$GCcount]\n"; }
#    print "numLocsForGCByPct = $numLocsForGCByPct[$GCcount] avgCoverage = $avgCoverage stdev = $stdev stdev/avgCoverage = $stdevToAvg\n";
}
for ($GCcount=50; $GCcount<=100; ++$GCcount) {
    if ($isBadValue[$GCcount]) {
	$lastGoodValue = $GCcount-1; 
	last; } }
if ($lastGoodValue !~ /\d/) {
    $lastGoodValue = 100; }
for ($GCcount=49; $GCcount>=0; --$GCcount) {
    if ($isBadValue[$GCcount]) {
	$firstGoodValue = $GCcount+1; 
	last; } }
if ($firstGoodValue !~ /\d/) {
    $firstGoodValue = 0; }

for ($GCcount=0; $GCcount<=100; ++$GCcount) {
    $GCpct = $GCcount/100;
    if ($GCcount < $firstGoodValue) {
	$adjustmentFactor = $adjustmentFactor[$firstGoodValue]; }
    elsif ($GCcount > $lastGoodValue) {
	$adjustmentFactor = $adjustmentFactor[$lastGoodValue]; }
    else {
	$adjustmentFactor = $adjustmentFactor[$GCcount]; }
    print "$GCpct $adjustmentFactor $outputLine[$GCcount]\n";
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

