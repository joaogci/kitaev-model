#!/bin/bash 

Variable=(Run Lx_R Ly_R alpha_R  Beta_R  NSW_R NB_R)
Name=( Run Lx_R Ly_R alpha_R  Beta_R  NSW_R NB_R )
#      0    1    2     3       4       5      6  
{
read -a Variable
echo ${Variable[@]}
echo ${Variable[0]}
while [ ! ${Variable[0]} = "stop" ];   do
    if [ ${Variable[0]} = "Y" ]; then
        export B_R_dir=`echo ${Variable[4]} | sed s/"\.0"//`
        
        export Dir="Lx"${Variable[1]}"_Ly"${Variable[2]}"_alpha"${Variable[3]}"_beta"${Variable[4]}        
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
        
        $KITAEV_DIR/src/main 

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

