package MasurcaSoap;

sub runSOAP {
  my ($out, $reads_file, %config) = @_;
  
  my $cmdline="splitFileByPrefix.pl";
  foreach my $v (@{$config{JUMP_INFO}}) {
    my @f = split(" ", $v);
    $cmdline .= " $f[0]";
  }

  print $out <<"EOS";
$cmdline < sj.cor.clean2.fa
log 'SOAPdenovo'
mkdir -p SOAP_assembly
( cd SOAP_assembly
  if [ \$KMER -le 63 ];then
    SOAPdenovo-63mer all -u -w -p $config{NUM_THREADS} -D 0 -d 0 -K \$KMER -k 33 -R -o asm -s ../soap_config 1>../SOAPdenovo.err 2>\&1
  else
    SOAPdenovo-127mer all -u -w -p $config{NUM_THREADS} -D 0 -d 0 -K \$KMER -k 33 -R -o asm -s ../soap_config 1>../SOAPdenovo.err 2>\&1
  fi
)
EOS
  
  if($config{CLOSE_GAPS}){
    my $reads_argument= join(" ", map { "--reads-file '$_'" } @$reads_file);
    print $out <<"EOS";
log 'Gap closing'
if [[ -e \"SOAP_assembly/asm.scafSeq\" ]];then
  closeGapsInScaffFastaFile.perl  --max-reads-in-memory 1000000000 -s $config{JF_SIZE} --scaffold-fasta-file  SOAP_assembly/asm.scafSeq $reads_argument --output-directory SOAP_gapclose --min-kmer-len 19 --max-kmer-len \$((\$PE_AVG_READ_LENGTH-5)) --num-threads $config{NUM_THREADS} --contig-length-for-joining \$((\$PE_AVG_READ_LENGTH-1)) --contig-length-for-fishing 200 --reduce-read-set-kmer-size 25 1>gapClose.err 2>&1
  if [[ -e \"SOAP_gapclose/genome.ctg.fasta\" ]];then
    echo \"Gap close success. Output sequence is in SOAP_gapclose/genome.\{ctg,scf\}.fasta\"
  else
    fail Gap close failed, you can still use pre-gap close scaffold file asm.scafSeq. Check gapClose.err for problems.
  fi
else
  fail SOAPdenovo failed, Check SOAPdenovo.err for problems.
fi
EOS
  }
}

1;
