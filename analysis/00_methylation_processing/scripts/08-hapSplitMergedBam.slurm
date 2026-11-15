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

# create list of files to clean up, even upon failure. NOTE: not actually being
# used anymore now that we are just using the bed in the parent directory.
declare -a TO_CLEAN=()

# handle signals
cleanup()    { rm -f "${TO_CLEAN[@]}" &> /dev/null; }
unfinished() { printf 'ERROR: Failure while splitting by haplotype the merged sorted BAM from %s-ing %s'"'"'s %s piece(s).\n' "${MAPPER:-UNKNOWN_MAPPER}" "${SAMPLE:-UNKNOWN_SAMPLE}" "${DTYPE:-UNKNOWN_DATATYPE}" 1>&2; }
report_sig() { printf 'INFO: Exiting after receiving signal %u.\n' "${1}" 1>&2; }
hup_sgnl()   { unfinished; report_sig 1;  exit 1;  } # SIGHUP
control_c()  { unfinished; report_sig 2;  exit 2;  } # SIGINT
control_d()  { unfinished; report_sig 3;  exit 3;  } # SIGQUIT
term_sgnl()  { unfinished; report_sig 15; exit 15; } # SIGTERM
exiting()    { cleanup; exit ${1}; }                 # EXIT
trap  hup_sgnl      SIGHUP
trap  control_c     SIGINT
trap  control_d     SIGQUIT
trap  term_sgnl     SIGTERM
trap 'exiting ${?}' EXIT

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

# set the name of the merged BAM file
MERGED_BAM="${MAPPER%%2}_${DTYPE}-x-${REF_BASE}.F0x904.bam"
OK8="${MERGED_BAM%.bam}.split-haps.ok"

for HAP in {m,p}at hap{1,2}
do
	# set the name of the hap-specific bam
	HAP_BAM="${MERGED_BAM%.bam}.${HAP}.bam"

	# determine the file with the bed of ref seqIDs for this hap
	HAP_SEQID_BED="../${REF_BASE}.${HAP}.bed"
	if [ ! -s "${HAP_SEQID_BED}" ]
	then
		printf 'Due to (likely expected) missing bed of reference seqIDs (%s), skipping subsetting merged BAM (%s) to include only reads mapped to %s (as defined by the (likely expected) missing bed of seqIDs). This is _not_ necessarily a problem (it likely is not) because we expect only mat/pat or hap1/hap2, so seeing this message for both hap1/hap2 or both mat/pat is expected. If you see it for a different combination, you should investigate.\n' "${HAP_SEQID_BED}" "${HAP}" "${MERGED_BAM}" 1>&2
		continue
	fi

	# only subset the merged BAM if we haven't done it already (based on .ok file)
	if [ -e "${HAP_BAM}.ok" ]
	then
		# report that we already did it (based on .ok file)
		printf 'Due to .ok file, skipping subsetting merged BAM (%s) to include only reads mapped to %s based on the reference seqIDs found in %s.\n' "${MERGED_BAM}" "${HAP}" "${HAP_SEQID_BED}" 1>&2
	else
		# we haven't already subset the merged BAM (based on .ok file),
		# so let's go ahead and subset it
		printf 'Subsetting merged BAM (%s) to include only reads mapped to %s based on the reference seqIDs found in %s...\n' "${MERGED_BAM}" "${HAP}" "${HAP_SEQID_BED}" 1>&2
		printf 'CMD: %s\n' "samtools view -O BAM --write-index -@${SAMTOOLS_EXTRA_THREADS} --regions-file ${HAP_SEQID_BED} -o ${HAP_BAM}##idx##${HAP_BAM}.bai ${MERGED_BAM}" 1>&2
		time samtools view \
				-O BAM \
				--write-index \
				-@${SAMTOOLS_EXTRA_THREADS} \
				--regions-file "${HAP_SEQID_BED}" \
				-o "${HAP_BAM}##idx##${HAP_BAM}.bai" \
				"${MERGED_BAM}"

		# we could do this instead to also set the @RG header and the RG tag for
		# each read to ${HAP}. If we did, we'd need to adjust our extra threads calc.
		#samtools view \
		#		-O BAM -u \
		#		-@${SAMTOOLS_EXTRA_THREADS} \
		#		--regions-file "${HAP_SEQID_BED}" \
		#		"${MERGED_BAM}" \
		#	| samtools addreplacerg \
		#		-O BAM \
		#		--write-index \
		#		-@${SAMTOOLS_EXTRA_THREADS} \
		#		-m overwrite_all \
		#		-w \
		#		-r "${HAP}" \
		#		-o "${HAP_BAM}##idx##${HAP_BAM}.bai"

		touch "${HAP_BAM}.ok"
	fi

done

# finally do a sanity check. We should have mat/pat or hap1/hap2.
#if ! ( ( [[ -s "${MERGED_BAM%.bam}.mat.bam" ]] && [[ -s "${MERGED_BAM%.bam}.pat.bam" ]] ) || ( [[ -s "${MERGED_BAM%.bam}.hap1.bam" ]] && [[ -s "${MERGED_BAM%.bam}.hap2.bam" ]] ) ) # <-- .bam files only
#if ! ( ( [[ -e "${MERGED_BAM%.bam}.mat.bam.ok" ]] && [[ -e "${MERGED_BAM%.bam}.pat.bam.ok" ]] ) || ( [[ -e "${MERGED_BAM%.bam}.hap1.bam.ok" ]] && [[ -e "${MERGED_BAM%.bam}.hap2.bam.ok" ]] ) ) # <-- .ok files only
if ! ( ( [[ -e "${MERGED_BAM%.bam}.mat.bam.ok"  ]] && [[ -e "${MERGED_BAM%.bam}.pat.bam.ok"  ]]   &&
		 [[ -s "${MERGED_BAM%.bam}.mat.bam"     ]] && [[ -s "${MERGED_BAM%.bam}.pat.bam"     ]] ) ||
	   ( [[ -e "${MERGED_BAM%.bam}.hap1.bam.ok" ]] && [[ -e "${MERGED_BAM%.bam}.hap2.bam.ok" ]]   &&
		 [[ -s "${MERGED_BAM%.bam}.hap1.bam"    ]] && [[ -s "${MERGED_BAM%.bam}.hap2.bam"    ]] ) )
then
	printf 'ERROR: Expected %s.{m,p}at.bam to both exist and be non-empty or %s.hap{1,2}.bam to both exist and be non-empty (plus existing .ok files in either case), but zero or only one existed in one or both pairs.\n' "${MERGED_BAM%.bam}" "${MERGED_BAM%.bam}" 1>&2
	exit 1
fi

touch "${OK8}"

#printf 'Exiting alignments directory...\n' 1>&2
#printf 'CMD: %s\n' "cd .." 1>&2
cd ..

