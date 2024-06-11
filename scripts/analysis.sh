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
        cd $Dir
        
        $KITAEV_DIR/src/ana *

        cd ..
    fi
    read -a Variable
    echo ${Variable[@]}
    echo ${Variable[0]}
done
}<sims

