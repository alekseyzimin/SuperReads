package MasurcaCelera;

# Celera commands. This function depends on many global variables from the masurca script!

sub runCA {
  my ($rerun_pe, $rerun_sj,$out, $frg_files, $tmplist, %config) = @_;
  my ($ovlMerSize, $other_parameters) = (30, "doFragmentCorrection=0 doOverlapBasedTrimming=0 doExtendClearRanges=0 ovlMerSize=30");
  if(@{$config{OTHER_INFO}} || @{$config{MOLECULO_INFO}}) {
    ($ovlMerSize, $other_parameters) = (22, "doFragmentCorrection=1 doOverlapBasedTrimming=1 doExtendClearRanges=1 ovlMerSize=22");
  }

  print $out <<"EOS";
rm -f CA/0-overlaptrim-overlap/overlap.sh CA/1-overlapper/overlap.sh
echo "computeInsertSize=0
cgwDemoteRBP=0
cnsReuseUnitigs=1
unitigger=bogart
utgGraphErrorLimit=1000
utgGraphErrorRate=0.015
utgMergeErrorLimit=1000
utgMergeErrorRate=0.025
utgErrorLimit=0
cnsConcurrency=$config{NUM_THREADS}" >runCA.spec
EOS

  if(not(-d "CA/7-0-CGW")|| $rerun_pe || $rerun_sj){
    if(not(-e "CA/5-consensus/consensus.success")|| $rerun_pe || $rerun_sj){
      if(not(-d "CA/genome.ovlStore")|| $rerun_pe || $rerun_sj){
        print $out <<"EOS";
TOTAL_READS=`cat  *.frg | grep -c --text '^{FRG' `
save TOTAL_READS
ovlRefBlockSize=`perl -e '\$s=int('\$TOTAL_READS'/8); if(\$s>100000){print \$s}else{print \"100000\"}'`
save ovlRefBlockSize
ovlCorrBatchSize=`perl -e '\$s=int('\$TOTAL_READS'/80); if(\$s>10000){print \$s}else{print \"10000\"}'`
save ovlCorrBatchSize
cnsMinFrags=`perl -e 'print int('\$TOTAL_READS'/256)'`
save cnsMinFrags
ovlMerThreshold=`jellyfish histo -t $config{NUM_THREADS} k_u_hash_0 | awk '{thresh=75;if(\$1>1) {dist+=\$2;if(dist>int(\"'\$ESTIMATED_GENOME_SIZE'\")*0.98&&flag==0){if(\$1>thresh) thresh=\$1;flag=1}}}END{print thresh}'`
log ovlMerThreshold=\$ovlMerThreshold\n
save ovlMerThreshold

echo "gkpFixInsertSizes=0
cnsMinFrags=\$cnsMinFrags
ovlRefBlockSize=\$ovlRefBlockSize 
ovlCorrBatchSize=\$ovlCorrBatchSize
ovlMerThreshold=\$ovlMerThreshold
merylThreads=$config{NUM_THREADS}
merylMemory=16384
ovlStoreMemory=16384
ovlThreads=2
frgCorrThreads=2
frgCorrConcurrency=$config{NUM_THREADS} 
ovlCorrConcurrency=$config{NUM_CNS_THREADS} 
ovlConcurrency=$config{NUM_THREADS}" >> runCA.spec
EOS
      }

print $out <<"EOS";
rm -rf CA/5-consensus 5-consensus-coverage-stat 5-consensus-insert-sizes
runCA -s runCA.spec stopAfter=consensusAfterUnitigger -p genome -d CA $config{CA_PARAMETERS} $other_parameters superReadSequences_shr.frg $frg_files   1> runCA1.out 2>&1

if [[ -e \"CA/4-unitigger/unitigger.err\" ]];then
  log \"Overlap/unitig success\"
else
  fail Overlap/unitig failed, check output under CA/ and runCA1.out
fi

recompute_astat_superreads.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt

NUM_SUPER_READS=`cat superReadSequences_shr.frg $tmplist | grep -c --text '^{FRG' `
save NUM_SUPER_READS
if [ ! -e CA/overlapFilter.success ];then
( cd CA &&  tigStore -g genome.gkpStore -t genome.tigStore 5 -d layout -U | tr -d '-' | awk 'BEGIN{print \">unique unitigs\"}{if(\$1 == \"cns\"){seq=\$2}else if(\$1 == \"data.unitig_coverage_stat\" && \$2>=5){print seq\"N\"}}' | jellyfish count -L 2 -C -m $ovlMerSize -s \$ESTIMATED_GENOME_SIZE -t $config{NUM_THREADS} -o unitig_mers /dev/fd/0 &&  cat <(overlapStore -b 1 -e \$NUM_SUPER_READS -d genome.ovlStore  | awk '{if(\$1<\$2 && \$1<'\$NUM_SUPER_READS' && \$2<'\$NUM_SUPER_READS') print \$0}'|filter_overlap_file -t $config{NUM_THREADS} <(gatekeeper  -dumpfragments -withsequence genome.gkpStore| grep -P '^fragmentIdent|^fragmentSequence' |perl -ane 'BEGIN{\$flag=1}{if(\$flag){print ">";}print "\$F[2]\\n";\$flag=1-\$flag;}') unitig_mers /dev/fd/0) <(overlapStore -d genome.ovlStore | awk '{if(\$1<\$2 && (\$1>='\$NUM_SUPER_READS' || \$2>='\$NUM_SUPER_READS')) print \$1\" \"\$2\" \"\$3\" \"\$4\" \"\$5\" \"\$6\" \"\$7}')  |convertOverlap -ovl |gzip -c > overlaps.ovb.gz &&  mkdir -p ovlStoreBackup && mv 4-unitigger 5-consensus 5-consensus-coverage-stat 5-consensus-insert-sizes genome.tigStore genome.ovlStore ovlStoreBackup && overlapStoreBuild -o genome.ovlStore -M 16384 -g genome.gkpStore overlaps.ovb.gz 1>overlapstore.err 2>&1 && rm overlaps.ovb.gz && touch overlapFilter.success)
fi

runCA -s runCA.spec stopAfter=consensusAfterUnitigger -p genome -d CA $config{CA_PARAMETERS} $other_parameters superReadSequences_shr.frg $frg_files   1> runCA1.out 2>&1

recompute_astat_superreads.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt
EOS
    }
  }
#scaffolding
  print $out <<"EOS";
runCA -s runCA.spec $config{CA_PARAMETERS} -p genome -d CA $other_parameters 1>runCA3.out 2>&1
if [[ -e \"CA/9-terminator/genome.qc\" ]];then
  log \"CA success\"
else
  fail CA failed, check output under CA/ and runCA3.out
fi
EOS
}

1;
