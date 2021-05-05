#!/bin/bash -f
# n runs
let a=0
let b=0
let i=0
let j=0
let k=1
let m=0
let n=0
# number of concurrent jobs to run on a server
# should ideally be either 16 or 32 for salviati's 32-32.q
# 16 cores but 32 threads per server, YMMV depending on your job - so try both!
NWAY=16
# where the jobs live, change as you like
MYNAME="/scratch/mff/hxyquench_jobs/L16/T070_200/16runs1"
# base name for job
MYNAME2="v2run"
# input and binary
MYINP="initial.in"
MYBIN="hxyghost16beta.exe"
# n is number of runs
# set to 50 as per Lyon
# should ideally be a multiple of either 16 or 32 for salviati's 32-32.q
let n=16
echo $n" runs"
# output with these is useful down the road
echo "Got ${NSLOTS} slots."
IPWD=`pwd`
echo "in ${IPWD}"
# work out how many times to go around spawning batches of jobs
let m=`echo ${NWAY}`
# there's be 'a' n-way runs
let a=$n/$m
# and one 'b'-way run at the end
let b=$n%$m
# loop over runs
for i in `seq 1 $a`
do
# spawn the concurrent runs
  for j in `seq 1 $m`
  do
    MYDIR="${MYNAME}/${MYNAME2}${k}"
    $( mkdir -p ${MYDIR} ; cp ${MYBIN} ${MYINP} ${MYDIR} && cd ${MYDIR} ; export MYRANDOM=`python -c "import string; import random; random.seed(); print random.randint(100000000,199999999)"` ; sed -i.bak -e "s#^123456789#${MYRANDOM}#g" ${MYINP} ; ./${MYBIN} < ./${MYINP} ) &
    let k+=1
  done
# wait here until all concurrent runs are done
  wait
done
if [[ $b != 0 ]]
then
# spawn the remaining concurrent runs
  for j in `seq 1 $b`
  do
    MYDIR="${MYNAME}/${MYNAME2}${k}"
    $( mkdir -p ${MYDIR} ; cp ${MYBIN} ${MYINP} ${MYDIR} && cd ${MYDIR} ; export MYRANDOM=`python -c "import string; import random; random.seed(); print random.randint(100000000,199999999)"` ; sed -i.bak -e "s#^123456789#${MYRANDOM}#g" ${MYINP} ; ./${MYBIN} < ./${MYINP} ) &
    let k+=1
  done
# wait here until all runs are done
  wait
fi
exit 0
