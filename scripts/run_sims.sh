#!/bin/bash 

Variable=(Run N_R alpha_R  Beta_R  NSW_R NB_R)
Name=( Run N_R alpha_R  Beta_R  NSW_R NB_R )
#      0    1    2     3       4       5    
{
read -a Variable
echo ${Variable[@]}
echo ${Variable[0]}
while [ ! ${Variable[0]} = "stop" ];   do
    if [ ${Variable[0]} = "Y" ]; then
        export B_R_dir=`echo ${Variable[3]} | sed s/"\.0"//`
        
        export Dir="N"${Variable[1]}"_alpha"${Variable[2]}"_beta"${Variable[3]}        
        echo $Dir
        if [ ! -e $Dir ]; then
            mkdir $Dir
            cd  $Dir
            cp  ../Start/* . 
            cd ..
        else
            cd $Dir
            cp ../Start/parameters .
            cp ../Start/job.sh .
            cd ..
        fi
        cd $Dir
        
        let i=1
        while [  $i -lt 7 ]; do
            sed s/${Name[$i]}/${Variable[$i]}/    parameters  > tmp
            mv tmp parameters
            let i=i+1 
        done
        
        $KITAEV_DIR/src/main $KITAEV_DIR

        # sed s/Dir_R/$Dir/   job.sh  >  tmp
        # mv  tmp  job.sh
        # sed s/NCPU_R/${Variable[13]}/   job.sh  >  tmp
        # mv  tmp  job.sh
  
        # job_id=$(sbatch job.sh | sed -e s/[^0-9]//g)
        # echo "run $job_id with name $Dir"
        # for i in {1..1}; do
        #   job_id=$(sbatch --dependency=afterok:$job_id job.sh | sed -e s/[^0-9]//g)
        #   echo "run $job_id of dependency $i on $Dir"
        # done

        cd ..
    fi
    read -a Variable
    echo ${Variable[@]}
    echo ${Variable[0]}
done
}<sims

