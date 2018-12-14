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
echo "gkpFixInsertSizes=0
merylThreads=$config{NUM_THREADS}
merylMemory=32768
ovlStoreMemory=32768
ovlThreads=2
frgCorrThreads=2
frgCorrConcurrency=$config{NUM_THREADS} 
ovlCorrConcurrency=$config{NUM_THREADS} 
ovlConcurrency=$config{NUM_THREADS}
useGrid=$config{USE_GRID}
gridEngine=$config{GRID_ENGINE}
unitigger=bogart
utgGraphErrorLimit=1000
utgMergeErrorLimit=1000
utgGraphErrorRate=0.015
utgMergeErrorRate=0.025
ovlCorrBatchSize=100000
doUnitigSplitting=0
cgwDemoteRBP=0
doChimeraDetection=normal
merylThreads=$config{NUM_THREADS}
computeInsertSize=0
cnsOnGrid=0
cnsConcurrency=$config{NUM_THREADS}
cnsMinFrags=10000
cnsMaxCoverage=7
cnsReuseUnitigs=1
cgwErrorRate=0.1" > runCA.spec
EOS

  if(not(-d "CA/7-0-CGW")|| $rerun_pe || $rerun_sj){
    if(not(-e "CA/5-consensus/consensus.success")|| $rerun_pe || $rerun_sj){
      if(not(-d "CA/genome.ovlStore")|| $rerun_pe || $rerun_sj){
        print $out <<"EOS";
runCA -s runCA.spec stopAfter=initialStoreBuilding -p genome -d CA $config{CA_PARAMETERS} $other_parameters superReadSequences_shr.frg $frg_files   1> runCA0.out 2>&1

TOTAL_READS=`gatekeeper -dumpinfo -lastfragiid CA/genome.gkpStore | awk '{print \$NF}'`
save TOTAL_READS
ovlRefBlockSize=`perl -e '\$s=int('\$TOTAL_READS'/5); if(\$s>100000){print \$s}else{print \"100000\"}'`
save ovlRefBlockSize
ovlHashBlockLength=10000000
save ovlHashBlockLength
ovlCorrBatchSize=`perl -e '\$s=int('\$TOTAL_READS'/100); if(\$s>10000){print \$s}else{print \"10000\"}'`
save ovlCorrBatchSize

echo "ovlRefBlockSize=\$ovlRefBlockSize 
ovlHashBlockLength=\$ovlHashBlockLength
ovlCorrBatchSize=\$ovlCorrBatchSize
" >> runCA.spec
EOS
      }

print $out <<"EOS"; 
runCA -s runCA.spec stopAfter=consensusAfterUnitigger -p genome -d CA $config{CA_PARAMETERS} $other_parameters superReadSequences_shr.frg $frg_files   1> runCA1.out 2>&1
rm -f CA/5-consensus/consensus.sh
if [[ -e \"CA/4-unitigger/unitigger.err\" ]];then
  log \"Overlap/unitig success\"
else
  fail Overlap/unitig failed, check output under CA/ and runCA1.out
fi
log \"Recomputing A-stat for super-reads\"
recompute_astat_superreads_CA8.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt superReadSequences_shr.frg

if [ ! -e CA/overlapFilter.success ];then
NUM_SUPER_READS=`cat superReadSequences_shr.frg $tmplist | grep -c --text '^{FRG' `
save NUM_SUPER_READS
log \"Filtering overlaps\"
( cd CA && \
tigStore -g genome.gkpStore -t genome.tigStore 5 -U -d consensus | \
awk -F "=" 'BEGIN{print ">unique unitigs";flag=0}{if(\$1 ~ /^>/){if(\$6>=5){flag=1}}else{if(flag){print \$1"N"}flag=0}}' | \
jellyfish count -L 2 -C -m $ovlMerSize -s \$ESTIMATED_GENOME_SIZE -t $config{NUM_THREADS} -o unitig_mers /dev/fd/0 && \
cat <(overlapStore -b 1 -e \$NUM_SUPER_READS -d genome.ovlStore  | awk '{if(\$1<\$2 && (\$1<'\$NUM_SUPER_READS' && \$2<'\$NUM_SUPER_READS')) print \$0}'|filter_overlap_file -t $config{NUM_THREADS} <(gatekeeper  -dumpfragments -withsequence genome.gkpStore| grep -P '^fragmentIdent|^fragmentSequence' | awk 'BEGIN{flag=1}{if(flag){print ">"\$3}else{ print \$3;} flag=1-flag; }') unitig_mers /dev/fd/0) <(overlapStore -d genome.ovlStore | awk '{if(\$1<\$2 && (\$1>='\$NUM_SUPER_READS' || \$2>='\$NUM_SUPER_READS')) print \$1\" \"\$2\" \"\$3\" \"\$4\" \"\$5\" \"\$6\" \"\$7}')  |convertOverlap -ovl |gzip > overlaps_dedup.ovb.gz &&  \
overlapStoreBuild -o genome.ovlStore.BUILDING -M 32768 -g genome.gkpStore overlaps_dedup.ovb.gz 1>overlapStore.rebuild.err 2>&1 && \ 
rm -rf ovlStoreBackup && \
mkdir ovlStoreBackup && \
mv 4-unitigger 5-consensus 5-consensus-coverage-stat 5-consensus-insert-sizes genome.tigStore genome.ovlStore ovlStoreBackup && \
mv genome.ovlStore.BUILDING genome.ovlStore && \
touch overlapFilter.success
)

runCA -s runCA.spec stopAfter=consensusAfterUnitigger -p genome -d CA $config{CA_PARAMETERS} $other_parameters superReadSequences_shr.frg $frg_files   1> runCA2.out 2>&1
log \"Recomputing A-stat for super-reads\"
recompute_astat_superreads_CA8.sh genome CA \$PE_AVG_READ_LENGTH work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt superReadSequences_shr.frg
fi

EOS
    }
  }
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
