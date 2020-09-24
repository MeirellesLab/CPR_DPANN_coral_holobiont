#!/bin/bash
##Version: Jan 4, 2019 - 16:09

initdir=$PWD
exec 2> /dev/null

while getopts s:d:t:k:o:f:c: option
do
case "${option}"
in
s) SAMPLE=${OPTARG};;
d) DIR=(${OPTARG});;
t) SRAPATH=${OPTARG};;
k) DBPATH=$OPTARG;;
o) OUTDIR=${OPTARG};;
f) FILE=${OPTARG};;
c) c=${OPTARG};;

esac
done

exec 2> /dev/tty

if [ "$1" == "--help" ] || [ "$1" == "-h" ]
then

echo "Taxon Profiling Script
This script was designed to process samples from NCBI (in the SRA format) and from the MGRAST database
(fastq or fasta format, beginning with mgm + id). The script runs all necessary steps for quality and
uniformity filtering of samples (a shorter pipeline for MGRAST samples, but longer for SRA files), and
executes taxonomic profiling using kraken2. For SRA files, the script runs vdb-validate, fastq-dump,
prinseq and kraken2. For MGRAST files, the script checks the file format, and then runs prinseq and
kraken. The script requires the file to be already downloaded. Since it requires no internet connection,
it can be readily submitted to a cluster partition for analysis. This script requires you to have prinseq
-lite (perl script without the .pl extension) and kraken2 (script created when installing kraken2) on your
bin folder, as executable. It also requires that SRA Toolkit is downloaded and decompressed in your computer.

Usage:
./taxonprofiling.sh [OPTIONS]

Another option is to remove the .sh extension and move the script to your bin folder, so it can be
immediately called in the command line, without having to explicit the location of the script.

OPTIONS

-s	Path to the  sample to be analyzed (when analysing only one sample, optional if -d is used)
-d	Path to the directory containing all samples to be analyzed. (Optional if -s is used)
-f	File format (either MGRAST or SRA). Mandatory.
-k	Path to the directory where the Kraken Reference Database to be used is located. Mandatory.
-t	Path to the directory where SRA Toolkit is located (without expliciting the 'bin' folder)
	mandatory if setting -f to SRA.
-o	Path to the directory where to place final output files. Optional (default is current directory)
	with one or more of the following: p (for Phyla), c (for Class), o (for Order), f (for Family),
	g (for Genus) or s (for Species). When using more than one option, they MUST be separated by a comma.
-c	Number of threads to be used in Kraken (default is 1)
--help	Shows this help information."


else
if [ -z "${SAMPLE}" ] && [ -z "${DIR}" ]
then
	echo 'You must provide path to sample (s). Run script with --help for more info.'
	exit 1
fi

if [ -z "${DBPATH}" ]
then
	echo 'You must provide path to kraken DB. Run script with --help for more info.'
	exit 1
fi

if [ -z "${OUTDIR}" ]
then
	OUTDIR=./
fi

if [ -z "${threads}" ]
then
	c=1
fi


##Checking type of file (MGRAST or SRA)
if [ -z "${FILE}" ]
then
	echo "You must provide file type (either MGRAST or SRA). Run script with --help for more info."
	exit 1
elif [ "${FILE}" == "SRA" ]
then

	##All lines below are for SRA files
	if [ -z "${SRAPATH}" ]
	then
		echo 'You must provide path to SRA Toolkit. Run script with --help for more info.'
		exit 1
	else
		checkvdb=$(find ${SRAPATH}/bin/ -name vdb-validate)
		checkfqdump=$(find ${SRAPATH}/bin/ -name fastq-dump)
		if [ -z "${checkvdb}" ] && [ -z "${checkfqdump}" ]
		then
			echo "The provided SRA path does not contain SRA Toolkit. Please provide a valid path!"
			exit 1
		elif [ -z "${checkvdb}" ] || [ -z "${checkfqdump}" ]
		then
			echo "The provided SRA path does not contain required files. Please provide a valid path!"
			exit 1
		fi
	fi

	if [ -z "${SAMPLE}" ]
	then
		cd ${DIR}
		samples=($( ls *.sra))
		cd ${initdir}
		for i in ${samples[@]}
		do
			#Retrieving SRA ID
			if [[ ${i} =~ SRR(.+).sra ]]
			then
				id=$(echo SRR${BASH_REMATCH[1]})
			fi
			if [[ ${i} =~ ERR(.+).sra ]]
			then
				id=$(echo ERR${BASH_REMATCH[1]})
			fi
			if [[ ${i} =~ DRR(.+).sra ]]
			then
				id=$(echo DRR${BASH_REMATCH[1]})
			fi

			##Fastq-dump
			if [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"*.fastq") == "" ]]
                        then
				echo "Running vdb validation..."
				${SRAPATH}/bin/vdb-validate ${DIR}/${i} > vdbout${id}.tmp
				if grep -q 'err' vdbout${id}.tmp
				then
					echo 'Sample '${id}' is corrupted. Please download it again used prefetch tool!'
					echo ${id} >> vdb_fail
				else
					echo "Sample "${id}" is okay."
					echo "Running fastq-dump..."
					${SRAPATH}/bin/fastq-dump --outdir ${OUTDIR}/fastqdump_output/ --skip-technical  --readids -R pass -B --split-files --clip ${DIR}/${i}
					echo 'Sample '${id}' processed by fastqdump, saved to '${OUTDIR}'/fastqdump_output/'
				fi
				rm vdbout${id}.tmp

				echo 'Modifying fastq file for prinseq...'
				N=$(ls ${OUTDIR}/fastqdump_output/${id}*.fastq | wc -l)
				if [ $N -eq 2 ]
				then
					sed -i 's/.\{2\}\s/ /' ${OUTDIR}/fastqdump_output/${id}_pass_1.fastq
					sed -i 's/.\{2\}\s/ /' ${OUTDIR}/fastqdump_output/${id}_pass_2.fastq
				else
					sed -i 's/.\{2\}\s/ /' ${OUTDIR}/fastqdump_output/${id}_pass_1.fastq
				fi
			fi
		done

		for i in ${samples[@]}
		do
			#Retrieving SRA ID
			if [[ ${i} =~ SRR(.+).sra ]]
			then
				id=$(echo SRR${BASH_REMATCH[1]})
			fi
			if [[ ${i} =~ ERR(.+).sra ]]
			then
				id=$(echo ERR${BASH_REMATCH[1]})
			fi
			if [[ ${i} =~ DRR(.+).sra ]]
			then
				id=$(echo DRR${BASH_REMATCH[1]})
			fi
			if [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"*.fastq") != "" ]]
			then
				echo "Checking if single or paired..."
				N=$(ls ${OUTDIR}/fastqdump_output/${id}*.fastq | wc -l)
				if [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"_prinseq_good*.fasta") == "" ]]
				then
					if [ $N -eq 2 ]
					then
						if  [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"*prinseq*.fasta") == "" ]]
						then
							echo "Running prinseq..."
							##Prinseq line for paired
							prinseq-lite -verbose -fastq ${OUTDIR}/fastqdump_output/${id}_pass_1.fastq -fastq2 ${OUTDIR}/fastqdump_output/${id}_pass_2.fastq -min_len 50 -min_qual_score 13 -trim_qual_left 13 -trim_qual_right 13 -ns_max_n 2 -out_format 1
							echo "Prinseq succesfully run on sample "${id}
						fi

						##Moving singletons to other directory and renaming prinseq output files
						mkdir ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq
						mkdir ${OUTDIR}/fastqdump_output/${id}_bad_prinseq
						mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq*singletons* ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq/
						mv ${OUTDIR}/fastqdump_output/${id}_pass_2_prinseq*singletons* ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq/
						mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq*bad* ${OUTDIR}/fastqdump_output/${id}_bad_prinseq/
						mv ${OUTDIR}/fastqdump_output/${id}_pass_2_prinseq*bad* ${OUTDIR}/fastqdump_output/${id}_bad_prinseq/
						mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq_good*.fasta ${OUTDIR}/fastqdump_output/${id}_prinseq_good_1.fasta
						mv ${OUTDIR}/fastqdump_output/${id}_pass_2_prinseq_good*.fasta ${OUTDIR}/fastqdump_output/${id}_prinseq_good_2.fasta
					else
						if  [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"*prinseq*.fasta") == "" ]]
						then
							echo "Running prinseq..."
							##Prinseq line for single
							prinseq-lite -verbose -fastq ${OUTDIR}/fastqdump_output/${id}_pass_1.fastq -min_len 50 -min_qual_score 13 -trim_qual_left 13 -trim_qual_right 13 -ns_max_n 2 -out_format 1
							echo "Prinseq succesfully run on sample "${id}
						fi
						##Moving singletons to other directory and renaming prinseq output files
						mkdir ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq
						mkdir ${OUTDIR}/fastqdump_output/${id}_bad_prinseq
						mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq*singletons* ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq/
						mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq*bad* ${OUTDIR}/fastqdump_output/${id}_bad_prinseq/
						mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq_good*.fasta ${OUTDIR}/fastqdump_output/${id}_prinseq_good_1.fasta
					fi
				fi
				if [[ $(find ${OUTDIR}/ -name ${id}"*_report") != "" ]]
				then
					if [ $N -eq 2 ]
					then
						##Kraken line for paired and good
						echo "Running Kraken..."
						kraken2 --db ${DBPATH} --paired ${OUTDIR}/fastqdump_output/${id}_prinseq_good_1.fasta ${OUTDIR}/fastqdump_output/${id}_prinseq_good_2.fasta --classified-out ${OUTDIR}/${id}_kraken_class# --unclassified-out ${OUTDIR}/${id}_kraken_unclass#  -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
						echo "Kraken successfully run for sample "${id}
					else
						##Kraken line for single
						echo "Running Kraken..."
						kraken2 --db ${DBPATH} ${OUTDIR}/fastqdump_output/${id}_prinseq_good_1.fasta --classified-out ${OUTDIR}/${id}_kraken_class --unclassified-out ${OUTDIR}/${id}_kraken_unclass  -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
						echo "Kraken successfully run for sample "${id}
					fi
				fi
			else
				echo "Fastq-dump for sample "${id}" failed (saved in file fastqdump_failed_samples. Please run the script again for this sample"
				echo ${id} >> fastqdump_failed_samples
			fi
		done
	else
		i=${SAMPLE}
		#Retrieving SRA ID
		if [[ ${i} =~ SRR(.+).sra ]]
		then
			id=$(echo SRR${BASH_REMATCH[1]})
		fi
		if [[ ${i} =~ ERR(.+).sra ]]
		then
			id=$(echo ERR${BASH_REMATCH[1]})
		fi
		if [[ ${i} =~ DRR(.+).sra ]]
		then
			id=$(echo DRR${BASH_REMATCH[1]})
		fi
		##Fastq-dump
		if [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"*.fastq") == "" ]]
		then
			echo "Running vdb validation..."
			${SRAPATH}/bin/vdb-validate ${DIR}/${i} > vdbout${id}.tmp
			if grep -q 'err' vdbout${id}.tmp
			then
				echo 'Sample '${id}' is corrupted. Please download it again used prefetch tool!'
				exit 1
			else
				echo "Sample "${id}" is okay."
				echo "Running fastq-dump..."
				${SRAPATH}/bin/fastq-dump --outdir ${OUTDIR}/fastqdump_output/ --skip-technical  --readids -R pass -B --split-files --clip ${i}
				echo 'Sample '${id}' processed by fastqdump, saved to '${OUTDIR}'/fastqdump_output/'

				echo 'Modifying fastq file for prinseq...'
				N=$(ls ${OUTDIR}/fastqdump_output/${id}*.fastq | wc -l)
				if [ $N -eq 2 ]
				then
					sed -i 's/.\{2\}\s/ /' ${OUTDIR}/fastqdump_output/${id}_pass_1.fastq
					sed -i 's/.\{2\}\s/ /' ${OUTDIR}/fastqdump_output/${id}_pass_2.fastq
				else
					sed -i 's/.\{2\}\s/ /' ${OUTDIR}/fastqdump_output/${id}_pass_1.fastq
				fi
			fi
			rm vdbout${id}.tmp
		fi

		if [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"*.fastq") != "" ]]
		then
			echo "Checking if single or paired..."
			N=$(ls ${OUTDIR}/fastqdump_output/${id}*.fastq | wc -l)
			if [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"_prinseq_good*.fasta") == "" ]]
			then
				if [ $N -eq 2 ]
				then
					if  [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"*prinseq*.fasta") == "" ]]
					then
						echo "Running prinseq..."
						##Prinseq line for paired
						prinseq-lite -verbose -fastq ${OUTDIR}/fastqdump_output/${id}_pass_1.fastq -fastq2 ${OUTDIR}/fastqdump_output/${id}_pass_2.fastq -min_len 50 -min_qual_score 13 -trim_qual_left 13 -trim_qual_right 13 -ns_max_n 2 -out_format 1
						echo "Prinseq succesfully run on sample "${id}
					fi

					##Moving singletons to other directory and renaming prinseq output files
					mkdir ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq
					mkdir ${OUTDIR}/fastqdump_output/${id}_bad_prinseq
					mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq*singletons* ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq/
					mv ${OUTDIR}/fastqdump_output/${id}_pass_2_prinseq*singletons* ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq/
					mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq*bad* ${OUTDIR}/fastqdump_output/${id}_bad_prinseq/
					mv ${OUTDIR}/fastqdump_output/${id}_pass_2_prinseq*bad* ${OUTDIR}/fastqdump_output/${id}_bad_prinseq/
					mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq_good*.fasta ${OUTDIR}/fastqdump_output/${id}_prinseq_good_1.fasta
					mv ${OUTDIR}/fastqdump_output/${id}_pass_2_prinseq_good*.fasta ${OUTDIR}/fastqdump_output/${id}_prinseq_good_2.fasta
				else
					if  [[ $(find ${OUTDIR}/fastqdump_output/ -name ${id}"*prinseq*.fasta") == "" ]]
					then
						echo "Running prinseq..."
						##Prinseq line for single
						prinseq-lite -verbose -fastq ${OUTDIR}/fastqdump_output/${id}_pass_1.fastq -min_len 50 -min_qual_score 13 -trim_qual_left 13 -trim_qual_right 13 -ns_max_n 2 -out_format 1
						echo "Prinseq succesfully run on sample "${id}
					fi
					##Moving singletons to other directory and renaming prinseq output files
					mkdir ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq
					mkdir ${OUTDIR}/fastqdump_output/${id}_bad_prinseq
					mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq*singletons* ${OUTDIR}/fastqdump_output/${id}_singleton_prinseq/
					mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq*bad* ${OUTDIR}/fastqdump_output/${id}_bad_prinseq/
					mv ${OUTDIR}/fastqdump_output/${id}_pass_1_prinseq_good*.fasta ${OUTDIR}/fastqdump_output/${id}_prinseq_good_1.fasta
				fi
			fi
			if [ $N -eq 2 ]
			then
				##Kraken line for paired and good
				echo "Running Kraken..."
				kraken2 --db ${DBPATH} --paired ${OUTDIR}/fastqdump_output/${id}_prinseq_good_1.fasta ${OUTDIR}/fastqdump_output/${id}_prinseq_good_2.fasta --classified-out ${OUTDIR}/${id}_kraken_class# --unclassified-out ${OUTDIR}/${id}_kraken_unclass#  -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
				echo "Kraken successfully run for sample "${id}
			else
				##Kraken line for single
				echo "Running Kraken..."
				kraken2 --db ${DBPATH} ${OUTDIR}/fastqdump_output/${id}_prinseq_good_1.fasta --classified-out ${OUTDIR}/${id}_kraken_class --unclassified-out ${OUTDIR}/${id}_kraken_unclass  -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
				echo "Kraken successfully run for sample "${id}
			fi
		else
			echo "Fastq-dump for sample "${id}" failed (saved in file fastqdump_failed_samples. Please run the script again for this sample"
			echo ${id} >> fastqdump_failed_samples
		fi
	fi
else
	##All lines below are for MGRAST files
	if [ -z "${SAMPLE}" ]
	then
		cd ${DIR}
		samples=($( ls mgm*))
		cd ${initdir}

		for i in ${samples[@]}
		do
			#Retrieving MGRAST ID
			if [[ $i =~ mgm(.+)_prinseq_good.fasta ]]
			then
				id=$(echo mgm${BASH_REMATCH[1]})
			fi

			if [[ $(find ${OUTDIR}/ -name ${id}"_kraken_report") == "" ]]
			then

			if [[ $(find ${OUTDIR}/ -name ${id}"*prinseq*") == "" ]]
			then
				echo "Checking if MGRAST file is fastq or fasta format..."
				if [ $(head -c 1 ${DIR}/${i}) == '>' ]
				then
					echo "File is in fasta format."
					if [[ $i =~ mgm(.+).fasta ]]
        	then
						echo "File already has fasta extension."
						mod=0
          else
						echo "Modifying extension..."
						mod=1
						mv ${DIR}/${i} ${DIR}/${i}.fasta
					fi
					echo "Running prinseq..."
					##Prinseq line for fasta
					if [ ${mod} -eq 1 ]
					then
						prinseq-lite -verbose -fasta ${DIR}/${i}.fasta -min_len 50 -ns_max_n 2 -out_format 1
						echo "Prinseq succesfully run on sample "${id}
					else
						prinseq-lite -verbose -fasta ${DIR}/${i} -min_len 50 -ns_max_n 2 -out_format 1
          	echo "Prinseq succesfully run on sample "${id}
					fi
			    mv ${DIR}/${id}*_prinseq_good*.fasta ${DIR}/${id}_prinseq_good.fasta

					##Kraken line
					echo "Running Kraken..."
					kraken2 --db ${DBPATH} ${DIR}/${id}_prinseq_good.fasta --classified-out ${OUTDIR}/${id}_kraken_class --unclassified-out ${OUTDIR}/${id}_kraken_unclass  -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
					echo "Kraken successfully run for sample "${id}

				elif [ $(head -c 1 ${DIR}/${i}) == "@" ]
				then
					echo "File is in fastq format."
					if [[ $i =~ mgm(.+).fastq ]]
					then
						echo "File already has fastq extension."
						mod=0
					else
						echo "Modifying extension..."
						mod=1
						mv ${DIR}/${i} ${DIR}/${i}.fastq
					fi
					echo "Running prinseq..."

					##Prinseq line for fasta
					if [ ${mod} -eq 1 ]
					then
						prinseq-lite -verbose -fastq ${DIR}/${i}.fastq -min_len 50 -trim_qual_right 13 -trim_qual_left 13 -min_qual_score 13 -ns_max_n 2 -out_format 1
						echo "Prinseq succesfully run on sample "${id}
					else
						prinseq-lite -verbose -fastq ${DIR}/${i} -min_len 50 -trim_qual_right 13 -trim_qual_left 13 -min_qual_score 13 -ns_max_n 2 -out_format 1
						echo "Prinseq succesfully run on sample "${id}
					fi
					mv ${DIR}/${id}*_prinseq_good*.fasta ${DIR}/${id}_prinseq_good.fasta

					##Kraken line
					echo "Running Kraken..."
					kraken2 --db ${DBPATH} ${DIR}/${id}_prinseq_good.fasta --classified-out ${OUTDIR}/${id}_kraken_class --unclassified-out ${OUTDIR}/${id}_kraken_unclass  -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
					echo "Kraken successfully run for sample "${id}

				else
					echo "Sample "${id}" is not in a recognized format. Supported formats are FASTA and FASTQ. Please, check file content!"
					echo "Sample "${id}" was not processed."
				fi
			else
				##Kraken line
        echo "Running Kraken..."
        kraken2 --db ${DBPATH} ${DIR}/${id}_prinseq_good.fasta --classified-out ${OUTDIR}/${id}_kraken_class --unclassified-out ${OUTDIR}/${id}_kraken_unclass -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
        echo "Kraken successfully run for sample "${id}
			fi
			fi
		done
	else
		i=${SAMPLE}
		if [[ $i =~ (.+)mgm[0-9]* ]]
		then
			path=${BASH_REMATCH[1]}
		fi
		#Retrieving MGRAST ID
		if [[ $i =~ mgm(.+).3 ]]
		then
			id=$(echo mgm${BASH_REMATCH[1]})
		fi
		if [[ $(find ${OUTDIR}/ -name ${id}"*prinseq*") == "" ]]
    then
			echo "Checking if MGRAST file is fastq or fasta..."
			if [ $(head -c 1 ${i}) == '>' ]
			then
				echo "File is in fasta format."
				if [[ $i =~ mgm(.+).fasta ]]
        then
					echo "File already has fasta extension."
					mod=0
        else
					echo "Modifying extension..."
					mod=1
					mv ${i} ${i}.fasta
				fi
				echo "Running prinseq..."
				##Prinseq line for fasta
				if [ ${mod} -eq 1 ]
				then
					prinseq-lite -verbose -fasta ${i}.fasta -min_len 50 -ns_max_n 2 -out_format 1
					echo "Prinseq succesfully run on sample "${id}
				else
					prinseq-lite -verbose -fasta ${i} -min_len 50 -ns_max_n 2 -out_format 1
         	echo "Prinseq succesfully run on sample "${id}
				fi
		    mv ${path}${id}*_prinseq_good*.fasta ${path}${id}_prinseq_good.fasta

				##Kraken line
				echo "Running Kraken..."
				kraken2 --db ${DBPATH} ${path}${id}_prinseq_good.fasta --classified-out ${OUTDIR}/${id}_kraken_class --unclassified-out ${OUTDIR}/${id}_kraken_unclass  -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
				echo "Kraken successfully run for sample "${id}

		elif [ $(head -c 1 ${i}) == "@" ]
		then
			echo "File is in fastq format."
			if [[ $i =~ mgm(.+).fastq ]]
			then
				echo "File already has fastq extension."
				mod=0
			else
				echo "Modifying extension..."
				mod=1
				mv ${i} ${i}.fastq
			fi
			echo "Running prinseq..."

			##Prinseq line for fastq
			if [ ${mod} -eq 1 ]
			then
				prinseq-lite -verbose -fastq ${i}.fastq -min_len 50 -trim_qual_left 13 -trim_qual_right 13 -min_qual_score 13 -ns_max_n 2 -out_format 1
				echo "Prinseq succesfully run on sample "${id}
			else
				prinseq-lite -verbose -fastq ${i} -min_len 50 -trim_qual_left 13 -trim_qual_right 13 -min_qual_score 13 -ns_max_n 2 -out_format 1
				echo "Prinseq succesfully run on sample "${id}
			fi
			mv ${path}${id}*_prinseq_good*.fasta ${path}${id}_prinseq_good.fasta

			##Kraken line
			echo "Running Kraken..."
			kraken2 --db ${DBPATH} ${path}${id}_prinseq_good.fasta --classified-out ${OUTDIR}/${id}_kraken_class --unclassified-out ${OUTDIR}/${id}_kraken_unclass  -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
			echo "Kraken successfully run for sample "${id}
		else
			echo "Sample "${id}" is not in a recognized format. Supported formats are FASTA and FASTQ. Please, check file content!"
			exit 1
		fi
		else
    ##Kraken line
    echo "Running Kraken..."
    kraken2 --db ${DBPATH} ${DIR}/${id}_prinseq_good.fasta --classified-out ${OUTDIR}/${id}_kraken_class --unclassified-out ${OUTDIR}/${id}_kraken_unclass  -output ${OUTDIR}/${id}_kraken_output --report ${OUTDIR}/${id}_kraken_report --threads ${c} --memory-mapping
    echo "Kraken successfully run for sample "${id}
	fi ##This is the fi to test if prinseq was already run
fi ##This is the fi to test whether one or more samples were supplied to mgrast
fi ##This is the fi to test whether it was MGRAST or SRA
fi ##This is the fi for the help information
