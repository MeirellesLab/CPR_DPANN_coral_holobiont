#!/bin/bash
#Script by Maria Beatriz Walter Costa and adapted by Leticia Costa
#This script downloads SRA files, from the IDs of array @sra_list

#Usage: $bash SCRIPT
#To adapt the script for your needs, remove the old IDs of the array, and insert your own

#Array of NCBI file IDs
list=(SRR1275409
SRR1275449
SRR1822516
SRR2937345
SRR2937346
SRR2937347
SRR2937348
SRR2937349
SRR2937350
SRR2937351
SRR2937352
SRR2937353
SRR2937354
SRR2937355
SRR2937356
SRR3499156
SRR3569370
SRR3694369
SRR3694370
SRR3694371
SRR3694372
SRR5215424
SRR5215455
SRR5215457
SRR5215462
SRR5605611
SRR6784973
SRR6784974
SRR6784975
SRR6784976
SRR6784977
SRR6784978
SRR6784979
SRR6784980
SRR6784981
SRR6784982
SRR6784983
SRR6784984
SRR6784985
SRR6784986
SRR6784987
SRR6784988
SRR6784989
SRR6784990
SRR6784991
SRR6784992
SRR6784993
SRR6784994
SRR6784995
SRR6784996
SRR6784997
SRR6784998
SRR6784999
SRR6785000
SRR6785005
SRR6785006
SRR6785009
SRR6785010
SRR6785011
SRR6785012
SRR6785013
SRR6785014
SRR6785015
SRR6785017
SRR6785018
SRR6785019
SRR6785020
SRR6785021
SRR6785022
SRR6785023
SRR6785024
SRR6785026
SRR6785027
SRR6785028
SRR6785029
SRR6785030
SRR6785031
SRR6785032
SRR6785033
SRR6785034
SRR6785035
SRR6785036
SRR6785037
SRR6785038
SRR6785039
SRR6785040
SRR6785041
SRR6785042
SRR6785043
SRR6785044
SRR6785045
SRR6785046
SRR6785047
SRR6785048
SRR6785049
SRR6785050
SRR6785051
SRR6785052
SRR6785053
SRR6785054
SRR6785055
SRR6785056
SRR6785057
SRR6785058
)

#This following block of code is a loop to download each of the NCBI file IDs
for i in "${list[@]}"; do
	
	#Variable accession is the NCBI ID
	accession=$i
	#Variable sixletters is the first six letters of variable accession
	sixletters=${accession:0:6}

	#Testing the variables
	#echo "File ${i}.sra"
	#echo "Accession $accession"
	#echo "Sixletters $sixletters"

	#After getting the variables, the command line to download the data goes below 
	#NCBI help-page: https://www.ncbi.nlm.nih.gov/books/NBK158899/ - /sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/<first 6 characters of accession>/<accession>/<accession>.sra
	#curl ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$sixletters/$accession/${accession}.sra
	wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$sixletters/$accession/${accession}.sra"

	#Old command line below
	#curl http://api.metagenomics.anl.gov/download/"${i}"?file=299.1 > ${i}.299.1.gz
	#echo "$i has been processed"
done
