#! /usr/bin/env bash

set -euo pipefail

# figure out basic resources
THREADS=1
MEM_CPU_IN_MB=0
TMP="/tmp/${$}"
declare -a SPLIT_NUMS
RANGETOOL="rangeTool.py"
if [ -n "${SLURM_JOB_ID:-}" ] # basically, if this is managed by slurm vs being run locally
then
	# figure out threads
	THREADS="${SLURM_NTASKS:-0}"
	if [ ${THREADS} -le 1 ]
	then
		THREADS="${SLURM_CPUS_PER_TASK:-1}"
	fi

	# figure out mem
	if [ -n "${SLURM_JOB_NUM_NODES:-}" ] && [ ${SLURM_JOB_NUM_NODES} -ne 1 ]
	then
		printf "%s\n" "ERROR: This job is meant to be run with a single node" 1>&2
		exit 1
	elif [ -n "${SLURM_MEM_PER_CPU:-}" ]
	then
		MEM_CPU_IN_MB=${SLURM_MEM_PER_CPU}
	elif [ -n "${SLURM_MEM_PER_NODE:-}" ]
	then
		MEM_CPU_IN_MB=$(bc <<< "${SLURM_MEM_PER_NODE}*${SLURM_JOB_NUM_NODES:-1}/${THREADS}")
	else
		printf "%s\n" 'ERROR: $SLURM_MEM_PER_NODE and $SLURM_MEM_PER_CPU not specificed.' 1>&2
		exit 1
	fi

	# figure out lscratch/tmp space
	TMP="/lscratch/${SLURM_JOB_ID}"
	
	# figure out the split number
	if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]
	then
		SPLIT_NUMS+=("${SLURM_ARRAY_TASK_ID}")
	fi

	# figure out the path to rangeTool.py
	if [ -n "${SLURM_SUBMIT_DIR:-}" ]
	then
		RANGETOOL="${SLURM_SUBMIT_DIR}/${RANGETOOL}"
	fi
fi

# figure out rangeTool.py location
if [ ! -s "${RANGETOOL}" ]
then
	RANGETOOL=`basename "${RANGETOOL}"`
	if [ ! -s "${RANGETOOL}" ]
	then
		if [ -s "../${RANGETOOL}" ]
		then
			RANGETOOL="../${RANGETOOL}"
		else
			set +e
			RANGETOOL=`which "${RANGETOOL}" 2> /dev/null`
			set -e
			if [ ! -n "${RANGETOOL:-}" ]
			then
				printf 'ERROR: Could not find a copy of rangeTool.py.\n' 1>&2
				exit 1
			fi
		fi
	fi
fi

# handle signals
unfinished() { printf 'ERROR: Failure while sorting (and indexing) BAM(s) from %s-ing %s'"'"'s %s, piece(s) %s.\n' "${MAPPER:-UNKNOWN_MAPPER}" "${SAMPLE:-UNKNOWN_SAMPLE}" "${DTYPE:-UNKNOWN_DATATYPE}" "${SPLIT_NUMS[*]:-UNKNOWN_PIECES}" 1>&2; }
report_sig() { printf 'INFO: Exiting after receiving signal %u.\n' "${1}" 1>&2; }
hup_sgnl()   { unfinished; report_sig 1;  exit 1;  } # SIGHUP
control_c()  { unfinished; report_sig 2;  exit 2;  } # SIGINT
control_d()  { unfinished; report_sig 3;  exit 3;  } # SIGQUIT
term_sgnl()  { unfinished; report_sig 15; exit 15; } # SIGTERM
trap hup_sgnl  SIGHUP
trap control_c SIGINT
trap control_d SIGQUIT
trap term_sgnl SIGTERM

# snag input from command-line args
if [ ${#} -lt 3 ]
then
	printf 'ERROR: incorrect args. Expected: 1= sample id (e.g., HG02055). 2= data type ("hifi" or "ont"). 3= aligner ("minimap2" or "winnowmap"). OPTIONAL: 4(+)= split numbers (e.g., 3 or 5-8 or 4,6-9 or 2,5-6 8 10-14).\n' 1>&2
	exit 1
fi

SAMPLE="${1}"
DTYPE="${2,,}"
if [ "${DTYPE}" == "nanopore" ]
then
	DTYPE="ont"
fi
MAPPER="${3,,}"
if [ "${MAPPER}" == "minimap" ]
then
	MAPPER="minimap2"
fi
shift 3
if [ ${#} -gt 0 ]
then
	if [[ -n "${SPLIT_NUMS:-}" ]] && [[ ${#SPLIT_NUMS[@]} -gt 0 ]]
	then
		printf 'INFO: overriding split(s) to process ("%s"; presumably originally from $SLURM_ARRAY_TASK_ID) with those provided on the command line: "%s".\n' "${SPLIT_NUMS[*]}" "${*}" 1>&2
	fi
	SPLIT_NUMS=("${@}")
fi
SPLIT_NUMS=(`python3 "${RANGETOOL}" -On -d, -D ' ' -N "${SPLIT_NUMS[@]}"`)
SPLIT_RANGES_STR=`python3 "${RANGETOOL}" -Or -d ' '  -D, -N "${SPLIT_NUMS[@]}"`
SPLIT_NUMS=(`printf '%s ' "${SPLIT_NUMS[@]}" | awk 'BEGIN{FS=" "; OFS=FS}{for(i=1;i<=NF;i++){$i=sprintf("%04u", $i)} print $0}'`)

# load modules and set some variables
REF_BASE="${SAMPLE}-verkkoV2.2.1-asmV4"

module purge
module load samtools/1.21  # samtools sort

if [[ "${MAPPER}" != "minimap2" ]] && [[ "${MAPPER}" != "winnowmap" ]]
then
	printf 'ERROR: You provided "%s" instead of minimap2 or winnowmap as the aligner (3rd command-line argument).\n' "${MAPPER}" 1>&2
	exit 1
fi

ALNS_DIR="alns"

SAMTOOLS_EXTRA_THREADS=$((${THREADS}-1))
if [ ${SAMTOOLS_EXTRA_THREADS} -le 0 ]
then
	SAMTOOLS_EXTRA_THREADS=0
fi

printf 'Moving into alignments directory...\n' 1>&2
printf 'CMD: %s\n' "cd ${ALNS_DIR}" 1>&2
cd "${ALNS_DIR}"

# sort and index the bams
printf 'Sequentially sorting (and indexing) all requested split BAMs (%s)...\n' "${SPLIT_RANGES_STR}" 1>&2
time for SPLIT_NUM in "${SPLIT_NUMS[@]}"
do
	IBAM="${MAPPER%%2}_${DTYPE}${SPLIT_NUM}-x-${REF_BASE}.F0x904.bam"
	OBAM="${IBAM%.bam}.sorted.bam"

	# set memory option
	MEM=${MEM_CPU_IN_MB}
	if [ ${MEM_CPU_IN_MB} -eq 0 ]
	then
		# guess the required memory (the same way I would when assigning mem to
		# the job) and hope for the best that we assigned the correct amount of
		# RAM to the job from the outset
		FSIZE=$(stat -L --printf '%s' "${IBAM}")
		MEM=$(bc <<< "(${FSIZE}*3.5/2^20+4096)/${THREADS}")
		unset FSIZE
	fi

	# do the sorting
	printf 'CMD: %s\n' "samtools sort --write-index -@${SAMTOOLS_EXTRA_THREADS} -m ${MEM}M -O bam -o ${OBAM}##idx##${OBAM}.bai ${IBAM}" 1>&2
	time samtools sort \
			--write-index \
			-@${SAMTOOLS_EXTRA_THREADS} \
			-m ${MEM}M \
			-O "bam" \
			-o "${OBAM}##idx##${OBAM}.bai" \
			"${IBAM}"
	
	touch "${OBAM}.ok"

	printf 'CMD: %s\n' "rm -f ${IBAM}" 1>&2
	rm -f "${IBAM}" &> /dev/null
	
done

#printf 'Exiting alignments directory...\n' 1>&2
#printf 'CMD: %s\n' "cd .." 1>&2
cd ..

