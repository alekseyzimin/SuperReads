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

#if ($#ARGV >= 0) {
#    $intervalLen = $ARGV[0]; }
while ($line = <STDIN>) {
    chomp ($line);
    ($GCcount, $coverage, $countForPoint) = split (" ", $line);
    next if ($coverage == 0);
    $numLocsForGC[$GCcount] += $countForPoint;
    $sumReadsCoverage[$GCcount] += ($coverage * $countForPoint);
    $sumSquareCoverage[$GCcount] += ($coverage * $coverage * $countForPoint);
}

for ($GCcount=0; $GCcount<=$#numLocsForGC; ++$GCcount) {
    if ($numLocsForGC[$GCcount] == 0) {
	next; }
    $GCpct = $GCcount / $intervalLen;
    $GCpercentage = int ($GCpct * 100);
    $diffLow = ($GCpct*100) - $GCpercentage;
    $diffHigh = $GCpercentage + 1 - $GCpct*100;
    $sumReadsCoverageByPct[$GCpercentage] += ($diffHigh * $sumReadsCoverage[$GCcount]);
    $sumSquareCoverageByPct[$GCpercentage] += ($diffHigh * $sumSquareCoverage[$GCcount]);
    $sumReadsCoverageByPct[$GCpercentage+1] += ($diffLow * $sumReadsCoverage[$GCcount]);
    $sumSquareCoverageByPct[$GCpercentage+1] += ($diffLow * $sumSquareCoverage[$GCcount]);
    $numLocsForGCByPct[$GCpercentage] += ($diffHigh * $numLocsForGC[$GCcount]);
    $numLocsForGCByPct[$GCpercentage+1] += ($diffLow * $numLocsForGC[$GCcount]);
}

$maxAvg = 0;
for ($GCcount=40; $GCcount<=60; ++$GCcount) {
    if ($numLocsForGCByPct[$GCcount] == 0) {
	next; }
    $avgCoverage = $sumReadsCoverageByPct[$GCcount] * 1.0 / $numLocsForGCByPct[$GCcount];
    if ($maxAvg < $avgCoverage) {
	$maxAvg = $avgCoverage; }
}

for ($GCcount=0; $GCcount<=$#numLocsForGCByPct; ++$GCcount) {
    if ($numLocsForGCByPct[$GCcount] == 0) {
        next; }
    $avgCoverage = $sumReadsCoverageByPct[$GCcount] * 1.0 / $numLocsForGCByPct[$GCcount];
    $avgCoverageSquared = $sumSquareCoverageByPct[$GCcount] / $numLocsForGCByPct[$GCcount];
    $variance = $avgCoverageSquared - $avgCoverage * $avgCoverage;
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

