#!/usr/bin/env perl
#
#assumes numeric super read names
$range=$ARGV[0];
$code=0;
$pr=-1;
$pp="";
$ps="";
%h=();
while($line=<STDIN>){
    chomp($line);
    @f=split(/\s+/,$line);
    $rn=int(substr($f[0],2));
    $ln=substr($f[0],0,2);
    if($rn%2==1 && $pr+1==$rn && $ln eq $ol){
	$sr_index_prev=$ps*25000+$po;
	$sr_index_curr=$f[1]*25000+$f[2];
	if($sr_index_prev>$sr_index_curr){
	    $tag=$ps." ".$f[1];
	    $coord1=$po;
	    $coord2=$f[2];
	}else{
	    $tag=$f[1]." ".$ps;
	    $coord2=$po;
	    $coord1=$f[2];
	}
	$code=0;
	for($i=-$range;$i<=$range;$i++){
	    for($j=-$range;$j<=$range;$j++){
		$code++ if(exists($h{$tag." ".($coord1+$i)." ".($coord2+$j)}));
	    }
	}
	if($code==0){
	    $h{$tag." ".($coord1)." ".($coord2)}=();
	}else{
	    print "$ln$rn\n$ln$pr\n";
	}
    }

    $ol=$ln;
    $pr=$rn;
    $ps=$f[1];
    $po=$f[2];
}


#awk '{rn=int(substr(\$1,3));if(rn%2==1 && int(substr(pr,3))+1==rn){print ps\" \"po\" \"pr\"\\n\"\$2\" \"\$3\" \"pr}else{pr=\$1;ps=\$2;po=\$3}}' |awk 'BEGIN{flag=0}{if(flag==1){index1=int(substr(c1_1,1,length(c1_1)-1))*20000+c1_2;index2=int(substr(\$1,1,length(\$1)-1))*20000+\$2;if(index1>index2){print c1_1\" \"\$1\" \"c1_2\" \"\$2\" \"c}else{print \$1\" \"c1_1\" \"\$2\" \"c1_2\" \"c}}c=\$3;c1_1=\$1;c1_2=\$2;flag=1-flag;}'|perl -ane '{chomp;\$range=2;\$code=0;for(\$i=-\$range;\$i<=\$range;\$i++){for(\$j=-\$range;\$j<=\$range;\$j++){\$code++ if(defined(\$h{\"\$F[0] \$F[1] \".(\$F[2]+\$i).\" \".(\$F[3]+\$j)}))}}if(\$code==0){\$h{\"\$F[0] \$F[1] \$F[2] \$F[3]\"}=1}else{print \"\$F[4]\\n\",substr(\$F[4],0,2),int(substr(\$F[4],2))+1,\"\\n\"}}' > redundant_sj.txt\n"
