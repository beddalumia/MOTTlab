#!/usr/bin/env bash
#
#SBATCH --job-name=TMIT             #(custom jobname)
#SBATCH --mail-user=user@hpc.it     #(e-mail address)
#SBATCH --mail-type=FAIL            #(e-mail request)        
#SBATCH --nodes=1                   #(shared-memory enforcement)
#SBATCH --ntasks=1                  #(MPI serial job)
#SBATCH --cpus-per-task=32          #(multithreading rank)         
#SBATCH --mem=0                     #(autodetermined memory)
#SBATCH --array=01-60%1             #(multiple jobs with parametric input)
#SBATCH --partition=queue.name      #(Multiple partitions are possible)
#SBATCH --time=03:05:07             #(Time limit < hrs:min:sec)
#SBATCH --output=sLOG_%x_out%A_%a   #(STDOUT-log)
#SBATCH --error=sLOG_%x_err%A_%a    #(STDERR-log)
#
# > MODULES AND LAUNCH DIRECTORY <
    module load matlab
    cd $SLURM_SUBMIT_DIR
#
# > JOB sLOG INFO PRINTS <
    echo 'JOBINFO:'
    echo 'job identifier is         '$SLURM_JOBID
    echo 'job name is               '$SLURM_JOB_NAME
    echo ""
    echo '------------------------------------------------------'
    echo 'This job is allocated on '$SLURM_JOB_CPUS_PER_NODE' cpu(s)'
    echo 'Job is running on node(s): '
    echo  $SLURM_JOB_NODELIST
    echo '------------------------------------------------------'
    START_TIME=`date +%H:%M-%a-%d/%b/%Y`
    echo ""
    echo 'WORKINFO:'
    echo 'job starting at           '$START_TIME
    echo 'sbatch is running on      '$SLURM_SUBMIT_HOST
    echo 'executing on cluster      '$SLURM_CLUSTER_NAME
    echo 'executing on partition    '$SLURM_JOB_PARTITION
    echo 'working directory is      '$SLURM_SUBMIT_DIR
    echo 'current home directory is '$(getent passwd $SLURM_JOB_ACCOUNT | cut -d: -f6)
    echo ""
    echo 'NODEINFO:'
    echo 'number of nodes is        '$SLURM_JOB_NUM_NODES
    echo 'number of cpus/node is    '$SLURM_JOB_CPUS_PER_NODE
    echo 'number of gpus/node is    '$SLURM_GPUS_PER_NODE
    echo '------------------------------------------------------'
#
# > JOB COMMAND SCRIPT <
    matlab -batch main -logfile "log_${SLURM_JOB_NAME}_U${SLURM_ARRAY_TASK_ID}e-1"
    cat "log_${SLURM_JOB_NAME}_U${SLURM_ARRAY_TASK_ID}e-1" >> "LOG_${SLURM_JOB_NAME}"
#
# Wait for processes, if any.
echo "Waiting for all the processes to finish..."
wait








