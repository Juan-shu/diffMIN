#!/bin/bash

# compute background data
nohup Rscript _1background_step1_step3.R &
  wait

nohup Rscript _2background_step4.R &
  wait

# compute foreground data
for i in "1" "2" "3" "4" "5" "6" "7"
do
	for step in _3foreground_step2_step3.R _4foreground_perturb_matrix_NMI.R _5foreground_denominator_matrix_NMI.R _6Compute_foreground_network.R _7foreground_step4.R
	do
		echo "Running ${step}"
		if [ "${step}" = "_3foreground_step2_step3.R" ];then
			sed -i "1s/.*/class <- \"g$i\"/" "${step}"
		fi
		
		if [ "${step}" = "_7foreground_step4.R" ];then
			sed -i "2s/.*/class <- \"g$i\"/" "${step}"
		fi
		nohup Rscript ${step} &
		wait
	done
done
