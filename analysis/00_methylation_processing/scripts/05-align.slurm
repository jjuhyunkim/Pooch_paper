#! /usr/bin/env bash

set -euo pipefail

# figure out basic resources
THREADS=8
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

# adjust threads if needed
if [ ${THREADS} -lt 2 ]
then
	THREADS=2
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
unfinished() { printf 'ERROR: Failure while %s-ing %s'"'"'s %s, piece(s) %s.\n' "${MAPPER:-UNKNOWN_MAPPER}" "${SAMPLE:-UNKNOWN_SAMPLE}" "${DTYPE:-UNKNOWN_DATATYPE}" "${SPLIT_NUMS[*]:-UNKNOWN_PIECES}" 1>&2; }
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
REF_FA="${REF_BASE}.fa.gz"
W_OPT=""

module purge
module unuse /usr/local/lmod/modulefiles
module load samtools/1.21  # samtools view

if [ "${MAPPER}" == "minimap2" ]
then
	module load minimap2/2.28  # minimap2

elif [ "${MAPPER}" == "winnowmap" ]
then
	module load winnowmap/2.03 # winnowmap
	W_OPT="-W ${REF_BASE}.k15.repetitive.tsv"
else
	printf 'ERROR: You provided "%s" instead of minimap2 or winnowmap as the aligner (3rd command-line argument).\n' "${MAPPER}" 1>&2
	exit 1
fi

X_OPT=""
if [ "${DTYPE}" == "hifi" ]
then
	X_OPT="map-pb"
	if [ "${MAPPER}" == "minimap2" ] # we know the only two possibilities at this point are "minimap2" and "winnowmap"
	then
		X_OPT="map-hifi"
	fi
elif [ "${DTYPE}" == "ont" ]
then
	X_OPT="map-ont"
else
	printf 'ERROR: You provided "%s" instead of hifi or ont as the data type (2nd command-line argument).\n' "${DTYPE}" 1>&2
	exit 1
fi

READS_DIR="reads/${DTYPE}"
ALNS_DIR="alns"

# create the output directory
printf 'CMD: %s\n' "mkdir -p ${ALNS_DIR} &> /dev/null" 1>&2
mkdir -p "${ALNS_DIR}" &> /dev/null

# set some more variables
ALN_THREADS=$((${THREADS}-2))
SAMTOOLS_EXTRA_THREADS=1
if [ ${THREADS} -le 8 ]
then
	SAMTOOLS_EXTRA_THREADS=0
	ALN_THREADS=${THREADS}
fi
if [ "${MAPPER}" == "winnowmap" ]
then
	ALN_THREADS=$((${ALN_THREADS}/2))
fi

# align the reads to the ref
printf 'Sequentially aligning all requested split read sets (%s) to the ref...\n' "${SPLIT_RANGES_STR}" 1>&2
time for SPLIT_NUM in "${SPLIT_NUMS[@]}"
do
	QRY_FQ="${READS_DIR}/${SPLIT_NUM}.fq.gz"
	OUT_BAM="${ALNS_DIR}/${MAPPER%%2}_${DTYPE}${SPLIT_NUM}-x-${REF_BASE}.F0x904.bam"

	printf 'CMD: %s\n' "${MAPPER} -x ${X_OPT} -I 8G --MD --cs -y -a -Y -L --eqx ${W_OPT} -t ${ALN_THREADS} ${REF_FA} ${QRY_FQ} | samtools view -F 0x904 -h -b -@${SAMTOOLS_EXTRA_THREADS} > ${OUT_BAM}" 1>&2
	time ${MAPPER} \
			-x "${X_OPT}" \
			-I 8G --MD --cs -y \
			-a -Y -L --eqx \
			${W_OPT} \
			-t ${ALN_THREADS} \
			"${REF_FA}" \
			"${QRY_FQ}" \
		| samtools view \
			-F 0x904 \
			-h -b \
			-@${SAMTOOLS_EXTRA_THREADS} \
		> "${OUT_BAM}"

	touch "${OUT_BAM}.ok"
done

