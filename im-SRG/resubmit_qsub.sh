#/bin/bash

export LANG=C 
#export BLCR_CHECKFILE="/mnt/scratch/parzuch6/CPOINTS/checkfile.blcr" 
#export BLCR_OUTFILE="output.txt"
#export BLCR_WAIT_TIME=10

if [ ! -f ${BLCR_CHECKFILE} ] 
then
    echo "Running for first time" 
    cr_run $* 1> ${BLCR_OUTFILE} 2>&1 &
    ### NOW CR_RUN IS RUNNING IN THE BACKGROUND
    export PID=$!   ## PID IS NOW THE PROCESS ID OF cr_run
else
    echo "Restarting ${BLCR_CHECKFILE}"
    cr_restart --no-restore-pid ${BLCR_CHECKFILE} >> ${BLCR_OUTFILE} 2>&1 &
    export PID=$!   ## PID IS NOW THE PROCESS ID OF cr_restart    
fi

(sleep ${BLCR_WAIT_TIME}; echo "done sleeping, checkpointing now."; cr_checkpoint -v -f ${BLCR_CHECKFILE} --term ${PID};) & 
timeout=$!
###  WRITE CHECKPOINT, Terminate PROCESS. 
echo "Running..."
wait ${PID} 
RET=$? #exit status of cr_run supposedly...
echo "Main process is done. Checking status" 

if [ ${RET} = "1" ] # if return status is "1" the job was terminated
then 
    echo "Job terminated mid-stride, waiting for checkpoint"
    wait ${timeout} 
    echo "Done, resubmitting..."
    cd ${PATH_TO_PBS_FILE}
    pwd
    qsub ${BLCR_PBS_FILE} 
else
    echo " Job finished before timer, killing it now" 
    kill ${timeout} #otherwise we just kill the timer
fi

echo "THE END STATE OF THE MAIN PROGRAM WAS: ${RET}" 

exit 0
