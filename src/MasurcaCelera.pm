package MasurcaCelera;

# Celera commands. This function depends on many global variables from the masurca script!

sub runCA {
  my ($rerun_pe, $rerun_sj,$out, $frg_files, $tmplist, %config) = @_;
  my ($ovlMerSize, $other_parameters) = (30, "doFragmentCorrection=0 doOverlapBasedTrimming=0 doExtendClearRanges=0 ovlMerSize=30 ovlMemory=4GB");
  if(@{$config{OTHER_INFO}} || @{$config{MOLECULO_INFO}}) {
    ($ovlMerSize, $other_parameters) = (22, "doFragmentCorrection=1 doOverlapBasedTrimming=1 doExtendClearRanges=1 ovlMerSize=22 ovlMemory=4GB");
  }

  print $out <<"EOS";
rm -f CA/0-overlaptrim-overlap/overlap.sh CA/1-overlapper/overlap.sh
echo "computeInsertSize=0
gkpFixInsertSizes=0
unitigger=bog
cnsConcurrency=$config{NUM_THREADS}" >runCA.spec
EOS

  if(not(-d "CA/7-0-CGW")|| $rerun_pe || $rerun_sj){
    if(not(-e "CA/5-consensus/consensus.success")|| $rerun_pe || $rerun_sj){
      if(not(-d "CA/genome.ovlStore")|| $rerun_pe || $rerun_sj){
        print $out <<"EOS";
runCA -s runCA.spec stopAfter=initialStoreBuilding -p genome -d CA $config{CA_PARAMETERS} $other_parameters superReadSequences_shr.frg $frg_files   1> runCA0.out 2>&1

TOTAL_READS=`gatekeeper -dumpinfo -lastfragiid CA/genome.gkpStore | awk '{print \$NF}'`
save TOTAL_READS
ovlRefBlockSize=`perl -e '\$s=int('\$TOTAL_READS'/10); if(\$s>100000){print \$s}else{print \"100000\"}'`
save ovlRefBlockSize
ovlHashBlockSize=`perl -e '\$s=int('\$TOTAL_READS'/100); if(\$s>10000){print \$s}else{print \"10000\"}'`
save ovlHashBlockSize
ovlCorrBatchSize=`perl -e '\$s=int('\$TOTAL_READS'/100); if(\$s>10000){print \$s}else{print \"10000\"}'`
save ovlCorrBatchSize
ovlMerThreshold=`jellyfish histo -t $config{NUM_THREADS} k_u_hash_0 | awk '{thresh=75;if(\$1>1) {dist+=\$2;if(dist>int(\"'\$ESTIMATED_GENOME_SIZE'\")*0.98&&flag==0){if(\$1>thresh) thresh=\$1;flag=1}}}END{print thresh}'`
log ovlMerThreshold=\$ovlMerThreshold\n
save ovlMerThreshold

echo "gkpFixInsertSizes=0
jellyfishHashSize=\$JF_SIZE
ovlRefBlockSize=\$ovlRefBlockSize 
ovlHashBlockSize=\$ovlHashBlockSize
ovlCorrBatchSize=\$ovlCorrBatchSize
ovlMerThreshold=\$ovlMerThreshold
merylThreads=$config{NUM_THREADS}
merylMemory=16384
ovlStoreMemory=16384
ovlThreads=1
frgCorrThreads=1
frgCorrConcurrency=$config{NUM_THREADS} 
ovlCorrConcurrency=$config{NUM_THREADS} 
ovlConcurrency=$config{NUM_THREADS}" >> runCA.spec
EOS
      }

print $out <<"EOS";
rm -f CA/5-consensus/consensus.sh 
runCA -s runCA.spec stopAfter=consensusAfterUnitigger -p genome -d CA $config{CA_PARAMETERS} $other_parameters superReadSequences_shr.frg $frg_files   1> runCA1.out 2>&1

if [[ -e \"CA/4-unitigger/unitigger.err\" ]];then
  log \"Overlap/unitig success\"
else
  fail Overlap/unitig failed, check output under CA/ and runCA1.out
fi

recompute_astat_superreads.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt

if [ ! -e CA/overlapFilter.success ];then
NUM_SUPER_READS=`cat superReadSequences_shr.frg $tmplist | grep -c --text '^{FRG' `
save NUM_SUPER_READS
( cd CA && tigStore -g genome.gkpStore -t genome.tigStore 2 -d layout -U | tr -d '-' | awk 'BEGIN{print \">unique unitigs\"}{if(\$1 == \"cns\"){seq=\$2}else if(\$1 == \"data.unitig_coverage_stat\" && \$2>=5){print seq\"N\"}}' | jellyfish count -L 2 -C -m $ovlMerSize -s \$ESTIMATED_GENOME_SIZE -t $config{NUM_THREADS} -o unitig_mers /dev/fd/0 &&  cat <(overlapStore -b 1 -e \$NUM_SUPER_READS -d genome.ovlStore  | awk '{if(\$1<'\$NUM_SUPER_READS' && \$2<'\$NUM_SUPER_READS') print \$0}'|filter_overlap_file -t $config{NUM_THREADS} <(gatekeeper -dumpfastaseq genome.gkpStore ) unitig_mers /dev/fd/0) <(overlapStore -d genome.ovlStore | awk '{if(\$1>='\$NUM_SUPER_READS' || \$2>='\$NUM_SUPER_READS') print \$1\" \"\$2\" \"\$3\" \"\$4\" \"\$5\" \"\$6\" \"\$7}')  |convertOverlap -b -ovl > overlaps.ovb &&  rm -rf 4-unitigger 5-consensus genome.tigStore && mkdir ovlStoreBackup && mv genome.ovlStore ovlStoreBackup && overlapStore -c genome.ovlStore -M 4096 -t $config{NUM_THREADS} -g genome.gkpStore overlaps.ovb 1>overlapstore.err 2>&1 && rm overlaps.ovb && touch overlapFilter.success
)

runCA -s runCA.spec stopAfter=consensusAfterUnitigger -p genome -d CA $config{CA_PARAMETERS} $other_parameters superReadSequences_shr.frg $frg_files   1> runCA2.out 2>&1

recompute_astat_superreads.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt
fi

EOS
    }
  }


#here we do the scaffolding, but first we have to check and fix unitig consensus

  print $out <<"EOS";
if [[ -e \"CA/5-consensus/consensus.success\" ]];then
  echo \"Unitig consensus success\"
else
  echo \"Fixing unitig consensus...\"
  mkdir CA/fix_unitig_consensus
  ( cd CA/fix_unitig_consensus
    cp `which fix_unitigs.sh` .
    ./fix_unitigs.sh genome 
   )
#last resort if fixing failed -- we simply delete the ones that were not fixed
echo \"Fixing unitig consensus... last resort\"
(cd CA && tigStore -g genome.gkpStore -t genome.tigStore 2 -U -d layout |grep -P '^unitig|^cns' |awk 'BEGIN{flag=0}{if(flag==0){unitig=\$2}else{if(length(\$2)==0) print \"tigStore -g genome.gkpStore -t genome.tigStore 2 -u \"unitig\" -d layout > layout.tmp && tigStore -g genome.gkpStore -t genome.tigStore 2 -R <(head -n 12 layout.tmp |awk \\47{if(\$1 ~ \\\"data.num_frags\\\") print \\\"data.num_frags 1\\\";else print \$0}\\47) && tail -n +13 layout.tmp |awk \\47{print \\\"frg iid \\\"\$5\\\" mateiid 0\\\"}\\47 > gkp.edits.msg && gatekeeper --edit gkp.edits.msg genome.gkpStore 1>gkp.out 2>&1 && awk \\47{if(\$5>0){print \\\"frg iid \\\"\$5\\\" mateiid 0\\\"}}\\47 gkp.out > gkp.edits1.msg && gatekeeper --edit gkp.edits1.msg genome.gkpStore 1>gkp1.out 2>&1\\n"; }flag=1-flag;}' > dump_delete_unitigs.sh && bash dump_delete_unitigs.sh && touch 5-consensus/consensus.success
)
fi

runCA -s runCA.spec $config{CA_PARAMETERS} -p genome -d CA $other_parameters 1>runCA3.out 2>&1
if [[ -e \"CA/9-terminator/genome.qc\" ]];then
  log \"CA success\"
else
  fail CA failed, check output under CA/ and runCA3.out
fi
EOS
}

1;
