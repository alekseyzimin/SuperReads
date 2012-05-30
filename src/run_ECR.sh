BEGIN_SCF=$1;
END_SCF=$2;
LAST_CKP=$3;
PREFIX=$4;

echo "extendClearRanges  -g ../$PREFIX.gkpStore  -t ../$PREFIX.tigStore  -n $LAST_CKP  -c $PREFIX  -b $BEGIN_SCF -e $END_SCF  -i 1  > extendClearRanges-scaffold.initial.err 2>&1"
extendClearRanges  -g ../$PREFIX.gkpStore  -t ../$PREFIX.tigStore  -n $LAST_CKP  -c $PREFIX  -b $BEGIN_SCF -e $END_SCF  -i 1  > extendClearRanges-scaffold.initial.err 2>&1

let NEXT_CKP=$LAST_CKP+1;
if [[ -e $PREFIX.ckp.$NEXT_CKP ]];then
echo "success"
touch extendClearRanges.success
exit;
else
echo "initial failure"
cp  extendClearRanges-scaffold.initial.err  extendClearRanges-scaffold.continue.err
while [ 1 ];do 
FAILED_SCF=`grep 'examining scaffold' extendClearRanges-scaffold.continue.err|tail -n 1|awk '{print substr($3,1,length($3)-1)}'`
let FAILED_SCF=$FAILED_SCF-1;
echo "extendClearRanges  -g ../$PREFIX.gkpStore -t ../$PREFIX.tigStore  -n $LAST_CKP  -c $PREFIX  -b $BEGIN_SCF -e $FAILED_SCF  -i 1  > extendClearRanges-scaffold.$FAILED_SCF.err 2>&1"
extendClearRanges  -g ../$PREFIX.gkpStore -t ../$PREFIX.tigStore  -n $LAST_CKP  -c $PREFIX  -b $BEGIN_SCF -e $FAILED_SCF  -i 1  > extendClearRanges-scaffold.$FAILED_SCF.err 2>&1
if [[ -e $PREFIX.ckp.$NEXT_CKP ]];then
echo "completed up to  $FAILED_SCF ok, continue";
else
echo "failed again before $FAILED_SCF";
exit;
fi
let BEGIN_SCF=$FAILED_SCF+2;
let LAST_CKP=$NEXT_CKP;
echo "extendClearRanges  -g ../$PREFIX.gkpStore  -t ../$PREFIX.tigStore  -n $LAST_CKP  -c $PREFIX  -b $BEGIN_SCF -e $END_SCF  -i 1  > extendClearRanges-scaffold.continue.err  2>&1"
extendClearRanges  -g ../$PREFIX.gkpStore  -t ../$PREFIX.tigStore  -n $LAST_CKP  -c $PREFIX  -b $BEGIN_SCF -e $END_SCF  -i 1  > extendClearRanges-scaffold.continue.err  2>&1
let NEXT_CKP=$LAST_CKP+1;
if [[ -e $PREFIX.ckp.$NEXT_CKP ]];then
touch extendClearRanges.success
exit
fi
#sleep 10
done
fi

