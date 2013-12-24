#!/usr/bin/perl
$minAStat = 5;
$minFromEnd = 200;
$shortestUnitigUsed = 2000;
$maxSlippagePerPiece = 5;
$prefix = "genome";
$dir = "/genome3/raid/alekseyz/rhodobacter/assembly2.1.0";
$CAdir = "$dir/CA";
$superReadDir = "$dir/work1";
$cmd = "gatekeeper -dumpfragments -tabular $CAdir/${prefix}.gkpStore | grep super-read |";
open (FILE, $cmd);
while ($line = <FILE>) {
    ($uidLocal, $iid, $len) = ($line =~ /^\s*(\S+)\s+(\S+)\s.*\s(\S+)\s*$/);
    $uid{$iid} = $uidLocal;
    $len{$uidLocal} = $len;
#    print "iid = $iid uid = $uidLocal len = $len{$uidLocal}\n";
}
close (FILE);

$fn = "$dir/unitig_cov.txt";
open (FILE, $fn);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $Astat{$flds[1]} = $flds[2];
}
close (FILE);

$cmd = "tigStore -g $CAdir/genome.gkpStore -t $CAdir/${prefix}.tigStore 2 -U -d layout | grep -v -E \"^(cns|qlt) \" |";
open (FILE, $cmd);
while ($line = <FILE>) {
    if ($line =~ /^FRG\s/) {
	($iid) = ($line =~ /^\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s/);
	next unless ($uid{$iid});
    }
    chomp ($line);
    @flds = split (" ", $line);
    if ($line =~ /^unitig/) {
	if (($len >= $shortestUnitigUsed) &&
	    ($Astat{$unitig} >= $minAStat)) {
	    &analyzeResults; }
	&clearResults;
	$unitig = $flds[1]; }
    elsif ($line =~ /^len\s/) {
	$len = $flds[1]; }
    next unless ($line =~ /^FRG/);
    $uid = $uid{$flds[4]};
    $begin = $flds[-2]; $end = $flds[-1];
#    print "$line\n";
    $Astat = $Astat{$unitig};
#    print "$unitig $len $Astat $uid $len{$uid} $begin $end\n";
    if ($uid =~ /\./) {
	($srNum,$srPiece) = ($uid =~ /^(\d+):.+\.(\d+)$/); }
    else {
	($srNum) = ($uid =~ /^(\d+):/); $srPiece = -1; }
    # $srPiece2 needed since must use for offset but srPiece = -1 if 1 piece
    $srPiece2 = $srPiece;
    if ($srPiece2 < 0) {
	$srPiece2 = 0; }
    if ($begin < $end) {
	$ori = "F";
	$impliedBegin = $begin - $srPiece2; }
    else {
	$ori = "R";
	$impliedBegin = $begin + $srPiece2; }
    push (@srNums, $srNum); push (@srPieces, $srPiece);
    push (@srNames, $uid); push (@begins, $begin); push (@ends, $end);
    push (@oris, $ori); push (@impliedBegins, $impliedBegin);
}
if (($len >= $shortestUnitigUsed) &&
    ($Astat{$unitig} >= $minAStat)) {
    &analyzeResults; }
&clearResults;
$fn = "$superReadDir/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt";
open (FILE, $fn);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $super = $flds[1];
    next unless ($unitigFromSuper{$super} =~ /\d/);
    if ($superOri{$super} eq "F") {
	$newOri = $flds[3];
	$beginOffset = $beginFromSuper{$super} + $flds[2]; }
    else {
	if ($flds[3] eq "F") {
	    $newOri = "R"; }
	else {
	    $newOri = "F"; }
	$beginOffset = $beginFromSuper{$super} - $flds[2]; }
    @newFlds = ($flds[0], $unitigFromSuper{$super}, $beginOffset, $newOri);
    print "@newFlds\n";
}
close (FILE);

sub analyzeResults
{

    my ($i, $im1, $index, $index0, $isGood, $super, $uid, $uidNum, $uidNumLast);
    my ($beg, $end);
    my (%hasFirstPiece, %hasLastPiece, %isSuperForThisUnitig);
    my (@hasNoFirstPiece, @hasNoLastPiece, @indices, @supers);
    @indices = (0..$#srNums);
    @indices = sort spcl @indices;
#    print "After sorting...\n";
#    for (@indices) {
#	$index = $_;
#	print "$unitig $len $Astat $srNames[$index] $len{$srNames[$index]} $begins[$index] $ends[$index]\n";
#    }
    undef %isSuperForThisUnitig; undef %hasFirstPiece; undef %hasLastPiece;
    undef %superLen; undef %superBegin; # undef %superOri;
    @hasNoFirstPiece = @hasNoLastPiece = ();
    $isGood = 1;
    for ($i=0; $i<=$#indices; ++$i) {
	$index = $indices[$i];
	$uid = $srNames[$index];
	$uidNum = $srNums[$index];
	$isSuperForThisUnitig{$uidNum} = 1;
	if ($srPieces[$index] <= 0) {
	    $superBegin{$uidNum} = $impliedBegins[$index];
	    $superOri{$uidNum} = $oris[$index];
	    $hasFirstPiece{$uidNum} = 1; }
	if ($len{$uid} < 2047) {
	    $superLen{$uidNum} = $srPieces[$index] + $len{$uid};
	    $hasLastPiece{$uidNum} = 1; }
	if ($srPieces[$index] < 0) {
	    $superLen{$uidNum} = $len{$uid};
	    $hasLastPiece{$uidNum} = 1; }
	if ($i > 0) {
	    $im1 = $i-1;
	    $index0 = $indices[$im1];
	    $uidNumLast = $srNums[$index0];
	    if ($uidNumLast ne $uidNum) {
		next; }
	    if ($oris[$index] ne $oris[$index0]) {
		push (@unitigFailsChangingOri, $uid); }
	    if (abs ($begins[$index] - $begins[$index0]) > 2047) {
		push (@unitigFailsForGap, $uid); }
	    if (abs ($impliedBegins[$index] - $impliedBegins[$index0]) > $maxSlippagePerPiece) {
		push (@unitigFailsForSlippage, $uid); }
	}
	if (($#unitigFailsForGap >= 0) ||
	    ($#unitigFailsForSlippage >= 0) ||
	    ($#unitigFailsChangingOri >= 0)) {
	    $isGood = 0;
	    last; }
    }
    if ($isGood) {
	@supers = keys %isSuperForThisUnitig;
	for (@supers) {
	    $super = $_;
	    if (! $hasFirstPiece{$super}) {
		$isGood = 0;
		push (@hasNoFirstPiece, $super); }
	    if (! $hasLastPiece{$super}) {
		$isGood = 0;
		push (@hasNoLastPiece, $super); }
	}
    }
#    print "Analyzing unitig $unitig...\n";
#    print "Super-reads with gaps: @unitigFailsForGap\n";
#    print "Super-reads with match slippage: @unitigFailsForSlippage\n";
#    print "Super-reads fail for changing ori: @unitigFailsChangingOri\n";
#    print "Super-read fails (missing first piece): @hasNoFirstPiece\n";
#    print "Super-read fails (missing last piece): @hasNoLastPiece\n";
    if (! $isGood) {
	return; }
    for (@supers) {
	$super = $_;
	$beg = $superBegin{$super};
	if ($superOri{$super} eq "F") {
	    $beg = $superBegin{$super};
	    $end = $beg + $superLen{$super}; }
	else {
	    $end = $superBegin{$super};
	    $beg = $end - $superLen{$super}; }
	
	print "$unitig $super $beg $end ";
	if ($superOri{$super} eq "F") {
	    $beginFromSuper{$super} = $beg;
	    $endFromSuper{$super} = $end;
	    print "0 $superLen{$super}\n"; }
	else {
	    $beginFromSuper{$super} = $end;
	    $endFromSuper{$super} = $beg;
	    print "$superLen{$super} 0\n"; }
	$unitigFromSuper{$super} = $unitig;
    }
}


sub clearResults
{
    @srNums = @srPieces = @srNames = @begins = @ends = @oris = @impliedBegins = 
	@unitigFailsChangingOri = @unitigFailsForGap = @unitigFailsForSlippage = ();
}

sub spcl
{
    if ($srNums[$a] != $srNums[$b]) {
	return ($srNums[$a] <=> $srNums[$b]); }
    return ($srPieces[$a] <=> $srPieces[$b]);
}

