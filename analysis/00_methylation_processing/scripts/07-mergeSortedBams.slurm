#! /usr/bin/env bash

set -euo pipefail

# figure out basic resources
THREADS=1
TMP="/tmp/${$}"
if [ -n "${SLURM_JOB_ID:-}" ] # basically, if this is managed by slurm vs being run locally
then
	# figure out threads
	THREADS="${SLURM_NTASKS:-0}"
	if [ ${THREADS} -le 1 ]
	then
		THREADS="${SLURM_CPUS_PER_TASK:-1}"
	fi

	# figure out lscratch/tmp space
	TMP="/lscratch/${SLURM_JOB_ID}"
fi

# handle signals
unfinished() { printf 'ERROR: Failure while merging sorted BAM(s) from %s-ing %s'"'"'s %s piece(s).\n' "${MAPPER:-UNKNOWN_MAPPER}" "${SAMPLE:-UNKNOWN_SAMPLE}" "${DTYPE:-UNKNOWN_DATATYPE}" 1>&2; }
report_sig() { printf 'INFO: Exiting after receiving signal %u.\n' "${1}" 1>&2; }
hup_sgnl()   { unfinished; report_sig 1;  exit 1;  } # SIGHUP
control_c()  { unfinished; report_sig 2;  exit 2;  } # SIGINT
control_d()  { unfinished; report_sig 3;  exit 3;  } # SIGQUIT
term_sgnl()  { unfinished; report_sig 15; exit 15; } # SIGTERM
trap hup_sgnl  SIGHUP
trap control_c SIGINT
trap control_d SIGQUIT
trap term_sgnl SIGTERM

# load modules
module purge
module load samtools/1.21  # samtools merge

# snag input from command-line args
if [ ${#} -ne 3 ]
then
	printf 'ERROR: incorrect args. Expected: 1= sample id (e.g., HG02055). 2= data type ("hifi" or "ont"). 3= aligner ("minimap2" or "winnowmap").\n' 1>&2
	exit 1
fi

SAMPLE="${1}"
DTYPE="${2,,}"
if [ "${DTYPE}" == "nanopore" ]
then
	DTYPE="ont"
fi
if [[ "${DTYPE}" != "hifi" ]] && [[ "${DTYPE}" != "ont" ]]
then
	printf 'ERROR: You provided "%s" instead of hifi or ont as the data type (2nd command-line argument).\n' "${DTYPE}" 1>&2
	exit 1
fi
MAPPER="${3,,}"
if [ "${MAPPER}" == "minimap" ]
then
	MAPPER="minimap2"
fi
if [[ "${MAPPER}" != "minimap2" ]] && [[ "${MAPPER}" != "winnowmap" ]]
then
	printf 'ERROR: You provided "%s" instead of minimap2 or winnowmap as the aligner (3rd command-line argument).\n' "${MAPPER}" 1>&2
	exit 1
fi

# set some variables
REF_BASE="${SAMPLE}-verkkoV2.2.1-asmV4"
ALNS_DIR="alns"

SAMTOOLS_EXTRA_THREADS=$((${THREADS}-1))
if [ ${SAMTOOLS_EXTRA_THREADS} -le 0 ]
then
	SAMTOOLS_EXTRA_THREADS=0
fi

printf 'Moving into alignments directory...\n' 1>&2
printf 'CMD: %s\n' "cd ${ALNS_DIR}" 1>&2
cd "${ALNS_DIR}"

MERGED_BAM="${MAPPER%%2}_${DTYPE}-x-${REF_BASE}.F0x904.bam"
SORTED_SPLIT_BAMS=(`find . -maxdepth 1 -type f -name "${MAPPER%%2}_${DTYPE}"'*'"-x-${REF_BASE}.F0x904.sorted.bam" -printf '%P\n' | sort -V`)

if [[ ! -n "${SORTED_SPLIT_BAMS:-}" ]] || [[ ${#SORTED_SPLIT_BAMS[@]} -eq 0 ]]
then
	printf 'ERROR: It looks like the requisite input BAM files could not be found. Searched for files matching glob %s in %s/.\n' \
			"${MAPPER%%2}_${DTYPE}"'*'"-x-${REF_BASE}.F0x904.sorted.bam" \
			"${ALNS_DIR}" \
			1>&2
	exit 1
fi

# merge the sorted split bams into a single sorted bam
printf 'Merging all sorted split BAMs...\n' 1>&2
printf '\t%s\n' "${SORTED_SPLIT_BAMS[@]}" 1>&2


if [[ -e "${MERGED_BAM}.ok" ]] && [[ -s "${MERGED_BAM}" ]]
then
	printf 'INFO: Quitting (with exit code 0) without doing anything because the expected output file (%s) exists, as does the .ok file.\n' "${MERGED_BAM}" 1>&2
	exit 0
else
	if [ ${#SORTED_SPLIT_BAMS[@]} -gt 1 ]
	then
		printf 'CMD: %s\n' "samtools merge --write-index -@${SAMTOOLS_EXTRA_THREADS} -c -p -o ${MERGED_BAM}##idx##${MERGED_BAM}.bai ${SORTED_SPLIT_BAMS[*]}" 1>&2
		time samtools merge \
				--write-index \
				-@${SAMTOOLS_EXTRA_THREADS} \
				-c -p \
				-o "${MERGED_BAM}##idx##${MERGED_BAM}.bai" \
				"${SORTED_SPLIT_BAMS[@]}"

		touch "${MERGED_BAM}.ok"

		printf 'Cleaning up the sorted split BAMs we no longer need...\n' 1>&2
		printf 'CMD: rm -fv %s %s\n' "${SORTED_SPLIT_BAMS[*]}" "${SORTED_SPLIT_BAMS[*]/%/.bai}" 1>&2
		rm -fv "${SORTED_SPLIT_BAMS[@]}" "${SORTED_SPLIT_BAMS[@]/%/.bai}" &> /dev/null

	else
		printf 'Only one input file exists. Merging does not make sense. Will simply rename the file (bam its index).\n' 1>&2
		SPLIT_BAM="${SORTED_SPLIT_BAMS[*]}"
		mv -v "${SPLIT_BAM}"     "${MERGED_BAM}"
		mv -v "${SPLIT_BAM}.bai" "${MERGED_BAM}.bai"
		touch "${MERGED_BAM}.ok"
		unset SPLIT_BAM
	fi
fi

#printf 'Exiting alignments directory...\n' 1>&2
#printf 'CMD: %s\n' "cd .." 1>&2
cd ..

