#$1 path to the stores
#$2 prefix of the stores
#$3 file with ALL jump library UIDs
#$4 threshold for the mate pair sizes -- I recommend setting it to the 1/2 of the library mean

echo "Filtering libraries for redundancy and chimerism...";

#dump all read uids
gatekeeper -dumpfragments -tabular $1/$2.gkpStore |awk '{print $1" "$3}' > $2.uidMuid

#create POSMAP - style unitig file containing all reads -- it may be useful later
tigStore -g $1/$2.gkpStore -t $1/$2.tigStore 1 -U -d layout | awk '{if($1 ~ /^u/){unitig=$2}else if($1 ~ /^F/){if($(NF-1)<$NF){print $5" "unitig" "$(NF-1)" "$NF" f";}else{print $5" "unitig" "$NF" "$(NF-1)" r";}}}' | perl -e '{open(FILE, $ARGV[0]);while($line=<FILE>){@l=split(/\s+/,$line);push(@uid,$l[0]);}while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);print "$uid[$l[0]] $l[1] $l[2] $l[3] $l[4]\n";}}' $2.uidMuid > $2.posmap.frgutg

#here we extract the mated reads from the relevant library and convert them into compatible read names for redundancy computation
awk '{if($2 !="0") print $0}' $2.uidMuid |  perl -e '{open(FILE,$ARGV[0]);while($line=<FILE>){chomp($line);@l=split(/\s+/,$line);$h{$l[0]}=1;};while($line=<STDIN>){@l=split(/\s+/,$line);print $line if(defined $h{$l[0]})}}' $3 | perl -e '{$rn=1;while($line=<STDIN>){chomp($line); @l=split(/\s+/,$line); next if(defined $used{$l[0]}); printf("%s il%012da\n%s il%012db\n",$l[0],$rn,$l[1],$rn);$used{$l[0]}=1;$used{$l[1]}=1;$rn++;}}' > $2.conversion

#convert the posmap file to the compatible names
perl -e '{open(FILE,$ARGV[0]);while($line=<FILE>){chomp($line);@l=split(/\s+/,$line);$newname{$l[0]}=$l[1]};while($line=<STDIN>){chomp($line);@flds=split(/\s+/,$line);next if(not(defined($newname{$flds[0]}))); print "$newname{$flds[0]} $flds[1]  $flds[2] $flds[3] $flds[4]\n";}}'  $2.conversion < $2.posmap.frgutg > $2.posmap.frgutg.compatible

#compute the redundancy
awk '{if($5=="f"){print $2"_"int($4/2)" "substr($1,1,14)}else{print $2"_"int($3/2)" "substr($1,1,14)}}' $2.posmap.frgutg.compatible | perl -e '{while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);@flds=split(/_/,$l[0]);$flag=$flds[0]+$flds[1];if(defined($mps{$l[1]})){@flds1=split(/_/,$fwd{$l[1]});if($flag<$flds1[0]+$flds1[1]){$rev{$l[1]}=$fwd{$l[1]};$fwd{$l[1]}=$l[0];}else{$rev{$l[1]}=$l[0];}}else{$fwd{$l[1]}=$l[0]}$mps{$l[1]}=1;}foreach $v(keys %mps){if(defined($fwd{$v}) && defined($rev{$v})){$se{"$fwd{$v} $rev{$v}"}.="$v "}elsif(defined($fwd{$v})){$se{"$fwd{$v} xxx"}.="$v "}elsif(defined($rev{$v})){$se{"xxx $rev{$v}"}.="$v "}} foreach $v(keys %se){@l=split(/\s+/,$se{$v});print $#l+1," $v $se{$v}\n";}}' |perl -e '{while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);for($i=4;$i<=$#l;$i++){print $l[$i],"\n";}}}'  > $2.inserts.redundant

awk '{if($5=="f"){print $2"_"int($3/2)" "substr($1,1,14)}else{print $2"_"int($4/2)" "substr($1,1,14)}}' $2.posmap.frgutg.compatible | perl -e '{while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);@flds=split(/_/,$l[0]);$flag=$flds[0]+$flds[1];if(defined($mps{$l[1]})){@flds1=split(/_/,$fwd{$l[1]});if($flag<$flds1[0]+$flds1[1]){$rev{$l[1]}=$fwd{$l[1]};$fwd{$l[1]}=$l[0];}else{$rev{$l[1]}=$l[0];}}else{$fwd{$l[1]}=$l[0]}$mps{$l[1]}=1;}foreach $v(keys %mps){if(defined($fwd{$v}) && defined($rev{$v})){$se{"$fwd{$v} $rev{$v}"}.="$v "}elsif(defined($fwd{$v})){$se{"$fwd{$v} xxx"}.="$v "}elsif(defined($rev{$v})){$se{"xxx $rev{$v}"}.="$v "}} foreach $v(keys %se){@l=split(/\s+/,$se{$v});print $#l+1," $v $se{$v}\n";}}' |perl -e '{while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);for($i=4;$i<=$#l;$i++){print $l[$i],"\n";}}}'  >> $2.inserts.redundant


#convert to UID's
perl -e '{open(FILE,"$ARGV[0]");while($line=<FILE>){chomp($line);$h{substr($line,0,14)}=1}while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);print $l[0],"\n" if(defined($h{substr($l[1],0,14)}));}}' $2.inserts.redundant < $2.conversion  > $2.redundant.uid

#report results
echo -n "Total paired reads in the assembly "
wc -l $2.conversion
echo -n "Redundant reads "
wc -l $2.redundant.uid

#now we figure out which inserts are too short/chimeric

awk '{print substr($1,1,14)" "$2;}'  $2.posmap.frgutg.compatible | sort -k1,1 -S 5% |uniq -d |awk '{print $1}' > $2.sameunitiginserts
perl -e '{open(FILE,$ARGV[0]); while($line=<FILE>){chomp($line);$h{$line}=1;}while($line=<STDIN>){@l=split(/\s+/,$line);next if(not defined($h{substr($l[0],0,length($l[0])-1)})); if($l[4] eq "f"){$h{substr($l[0],0,length($l[0])-1)}="$l[2] $l[4] $h{substr($l[0],0,length($l[0])-1)}";}else{$h{substr($l[0],0,length($l[0])-1)}="$h{substr($l[0],0,length($l[0])-1)} $l[3] $l[4]";}}foreach $v(keys %h){@l=split(/\s+/,$h{$v});print "$v $h{$v} ",$l[3]-$l[0],"\n" if($l[1] eq "f" && $l[4] eq "r");}}' $2.sameunitiginserts < $2.posmap.frgutg.compatible > $2.posmap.frgutg.sizeinfo

awk '{if($NF <'$4'){print $1}}' $2.posmap.frgutg.sizeinfo > $2.inserts.chimeric

#convert to UID's
perl -e '{open(FILE,"$ARGV[0]");while($line=<FILE>){chomp($line);$h{substr($line,0,14)}=1}while($line=<STDIN>){chomp($line);@l=split(/\s+/,$line);print $l[0],"\n" if(defined($h{substr($l[1],0,14)}));}}' $2.inserts.chimeric < $2.conversion  > $2.chimeric.uid

#report results
echo -n "Chimeric reads "
wc -l $2.chimeric.uid

