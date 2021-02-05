#!/bin/bash
#Script by Maria Beatriz Walter Costa and adapted by Leticia Costa
#This script downloads Mg-rast files, from the IDs of array @mgm_list

#Usage: $bash SCRIPT
#To adapt the script for your needs, remove the old IDs of the array, and insert your own

#Array of MG-RAST file IDs
list=(mgm4440319.3
mgm4445755.3
mgm4445756.3
mgm4480739.3
mgm4480740.3
mgm4480741.3
mgm4480742.3
mgm4480743.3
mgm4484839.3
mgm4486661.3
mgm4486662.3
mgm4486663.3
mgm4486664.3
mgm4486665.3
mgm4486666.3
mgm4486667.3
mgm4486668.3
mgm4486669.3
mgm4487909.3
mgm4487910.3
mgm4487911.3
mgm4516541.3
mgm4516694.3
mgm4694758.3
mgm4694759.3
mgm4694760.3
)

#This following block of code is a loop to download each of the MG-RAST file IDs
for i in "${list[@]}"; do

	wget https://api.mg-rast.org/download/"${i}"?file=299.1
	#curl https://api.mg-rast.org/download/"${i}"?file=299.1
	#echo "$i has been processed"
done
