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
unfinished() { printf 'ERROR: Failure while modkit pileup-ing merged sorted BAM from %s-ing %s'"'"'s %s piece(s).\n' "${MAPPER:-UNKNOWN_MAPPER}" "${SAMPLE:-UNKNOWN_SAMPLE}" "${DTYPE:-UNKNOWN_DATATYPE}" 1>&2; }
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
module load htslib/1.21     # tabix, bgzip
module load modkit/0.4.1    # modkit
module load bedtools/2.31.1 # bedtools

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

EXTRA_THREADS=$((${THREADS}-1))
if [ ${EXTRA_THREADS} -le 0 ]
then
	EXTRA_THREADS=0
fi

printf 'Moving into alignments directory...\n' 1>&2
printf 'CMD: %s\n' "cd ${ALNS_DIR}" 1>&2
cd "${ALNS_DIR}"

# set some more variables
REF_FA="../${REF_BASE}.fa"
#REF_FA="../${REF_BASE}.fa.gz" # modkit requires it to be uncompressed
MERGED_BAM="${MAPPER%%2}_${DTYPE}-x-${REF_BASE}.F0x904.bam"
BED_METHYL="${MERGED_BAM%.bam}.methyl.bed"
OK9="${MERGED_BAM%.bam}.modkit.ok"

# first run modkit pileup
if [ -e "${BED_METHYL}.ok" ]
then
	printf 'Skipping modkit pileup-ing the (sorted) merged BAM file (%s) because of .ok file.\n' "${MERGED_BAM}" 1>&2
else
	printf 'CMD: %s\n' "modkit pileup -t ${THREADS:-1} -n 51200 --preset traditional --ref ${REF_FA} --log-filepath ${BED_METHYL%.bed}.log ${MERGED_BAM} ${BED_METHYL}" 1>&2
	time modkit pileup \
		-t ${THREADS:-1} \
		-n 51200 \
		--preset traditional \
		--ref "${REF_FA}" \
		--log-filepath "${BED_METHYL%.bed}.log" \
		"${MERGED_BAM}" \
		"${BED_METHYL}"

	touch "${BED_METHYL}.ok"
fi

# then bgzip the resulting bedMethyl file
if [ -e "${BED_METHYL}.gz.ok" ]
then
	printf 'Skipping bgzipping the bedMethyl file (%s) because of .ok file.\n' "${BED_METHYL}" 1>&2
else
	printf 'CMD: %s\n' "bgzip -@${THREADS:-1} ${BED_METHYL}" 1>&2
	time bgzip -@${THREADS:-1} "${BED_METHYL}"
	touch "${BED_METHYL}.gz.ok"
fi


# finally, index the bedMethyl file w/ tabix
if [ -e "${BED_METHYL}.gz.tbi.ok" ]
then
	printf 'Skipping tabix indexing the bedMethyl file (%s) because of .ok file.\n' "${BED_METHYL}.gz" 1>&2
else
	printf 'CMD: %s\n' "tabix -@${EXTRA_THREADS:-0} -p bed ${BED_METHYL}.gz" 1>&2
	time tabix \
		-@${EXTRA_THREADS:-0} \
		-p bed \
		"${BED_METHYL}.gz"
	touch "${BED_METHYL}.gz.tbi.ok"
fi


# Now we'll split the bedMethyl file into two haplotypes (i.e., mat/pat or hap1/hap2)

for HAP in {m,p}at hap{1,2}
do
	# set the name of the hap-specific bedMethyl
	HAP_BED="${BED_METHYL%.bed}.${HAP}.bed"

	# determine the file with the bed of ref seqIDs for this hap
	HAP_SEQID_BED="../${REF_BASE}.${HAP}.bed"
	if [ ! -s "${HAP_SEQID_BED}" ]
	then
		printf 'Due to (likely expected) missing bed of reference seqIDs (%s), skipping subsetting bedMethyl (%s) to include only reads mapped to "%s" (as defined by the (likely expected) missing bed of seqIDs). This is _not_ necessarily a problem (it likely is not) because we expect only mat/pat or hap1/hap2, so seeing this message for both hap1/hap2 or both mat/pat is expected. If you see it for a different combination, you should investigate.\n' "${HAP_SEQID_BED}" "${BED_METHYL}" "${HAP}" 1>&2
		continue
	fi

	# only subset the bedMethyl if we haven't done it already (based on .ok file)
	if [ -e "${HAP_BED}.ok" ]
	then
		# report that we already did it (based on .ok file)
		printf 'Due to .ok file, skipping subsetting bedMethyl (%s) to include only reads mapped to %s based on the reference seqIDs found in %s.\n' "${BED_METHYL}" "${HAP}" "${HAP_SEQID_BED}" 1>&2
	else
		# we haven't already subset the bedMethyl (based on .ok file),
		# so let's go ahead and subset it

		# ALTERNATE METHOD 0
		# we're not going to use this method because then we needn't maintain the
		# regex in two places (i.e., this script and step 7)
		## build the search regex
		#PARTIAL_REGEX=""
		#if [[ "${HAP}" =~ ^[mp]at$ ]]
		#then
		#	PARTIAL_REGEX='(chr[0-9A-Za-z]+_)?'"${HAP}"'(ernal)?'
		#elif [[ "${HAP}" =~ ^hap[12]$ ]]
		#then
		#	PARTIAL_REGEX='hap(lotype)?'"${HAP##hap}"
		#else
		#	printf 'ERROR: unrecognized value for $HAP (%s). You should not see this error.\n' "${HAP}" 1>&2
		#	exit 1
		#fi

		#printf '[ALTERNATE 0] Subsetting bedMethyl (%s) to include only reads mapped to %s based on the reference seqIDs found in %s...\n' "${BED_METHYL}" "${HAP}" "${HAP_SEQID_BED}" 1>&2
		#printf 'CMD: %s\n' "bgzip -d -@${EXTRA_THREADS} < ${BED_METHYL}.gz | awk -v \"p=${PARTIAL_REGEX}\" 'BEGIN{FS=\"\\t\"; OFS=FS; r=\"^\" p \"-[0-9]+\$\"} \$1 ~ r {print \$0} /^(#|track|browser)/{print \$0}' > ${HAP_BED}" 1>&2
		#bgzip -d -@${EXTRA_THREADS} \
		#	< "${BED_METHYL}.gz" \
		#	| awk -v "p=${PARTIAL_REGEX}" \
		#		'BEGIN{FS="\t"; OFS=FS; r="^" p "-[0-9]+$"} $1 ~ r {print $0} /^(#|track|browser)/{print $0}' \
		#	> "${HAP_BED}"

		# ALTERNATE METHOD 1
		# We _are_ going to use this method. It used little RAM (similar to
		# method 0) and has the benefit of us not needing to maintin the regex
		# here in this script in addition to the other script for step 7.
		printf 'Subsetting bedMethyl (%s) to include only reads mapped to %s based on the reference seqIDs found in %s...\n' "${BED_METHYL}" "${HAP}" "${HAP_SEQID_BED}" 1>&2
		printf 'CMD: %s\n' "bedtools intersect -wa -u -header -nobuf -a ${BED_METHYL}.gz -b ${HAP_SEQID_BED} > ${HAP_BED}" 1>&2
		time bedtools intersect \
				-wa -u \
				-header \
				-nobuf \
				-a "${BED_METHYL}.gz" \
				-b "${HAP_SEQID_BED}" \
			> "${HAP_BED}"

		touch "${HAP_BED}.ok"

		# ALTERNATE METHOD 2
		# This version used a lot of RAM (for the sort, obviously), so we're
		# not going to use it because the input is already sorted.
		#printf '[ALTERNATE 2] Subsetting bedMethyl (%s) to include only reads mapped to %s based on the reference seqIDs found in %s...\n' "${BED_METHYL}" "${HAP}" "${HAP_SEQID_BED}" 1>&2
		#printf 'CMD: %s\n' "bedtools intersect -wa -u -header -nobuf -a ${BED_METHYL}.gz -b ${HAP_SEQID_BED} | bedtools sort -g <(cut -d $'\\t' -f 1-2 ${HAP_SEQID_BED} | sort -t $'\\t' -V -k1,1) -i - > ${HAP_BED%.bed}.alt2.bed" 1>&2
		#time bedtools intersect \
		#		-wa -u \
		#		-header \
		#		-nobuf \
		#		-a "${BED_METHYL}.gz" \
		#		-b "${HAP_SEQID_BED}" \
		#	| bedtools sort \
		#		-g <(cut -d $'\t' -f 1-2 "${HAP_SEQID_BED}" | sort -t $'\t' -V -k1,1) \
		#		-i - \
		#	> "${HAP_BED%.bed}.alt2.bed"

	fi

	# bgzip the hap-specific bedMethyl
	if [ -e "${HAP_BED}.gz.ok" ]
	then
		printf 'Skipping bgzipping the hap-specific bedMethyl file (%s) because of .ok file.\n' "${HAP_BED}" 1>&2
	else
		printf 'CMD: %s\n' "bgzip -@${THREADS:-1} ${HAP_BED}" 1>&2
		time bgzip -@${THREADS:-1} "${HAP_BED}"
		touch "${HAP_BED}.gz.ok"
	fi

	#  Index the bedMethyl file w/ tabix
	if [ -e "${HAP_BED}.gz.tbi.ok" ]
	then
		printf 'Skipping tabix indexing the hap-specific bedMethyl file (%s) because of .ok file.\n' "${HAP_BED}.gz" 1>&2
	else
		printf 'CMD: %s\n' "tabix -@${EXTRA_THREADS:-0} -p bed ${HAP_BED}.gz" 1>&2
		time tabix \
			-@${EXTRA_THREADS:-0} \
			-p bed \
			"${HAP_BED}.gz"
		touch "${HAP_BED}.gz.tbi.ok"
	fi

done

# finally do a sanity check. We should have mat/pat or hap1/hap2.
if ! ( ( [[ -e "${BED_METHYL%.bed}.mat.bed.ok"  ]] && [[ -e "${BED_METHYL%.bed}.pat.bed.ok"  ]]   &&
		 [[ -s "${BED_METHYL%.bed}.mat.bed.gz"  ]] && [[ -s "${BED_METHYL%.bed}.pat.bed.gz"  ]] ) ||
	   ( [[ -e "${BED_METHYL%.bed}.hap1.bed.ok" ]] && [[ -e "${BED_METHYL%.bed}.hap2.bed.ok" ]]   &&
		 [[ -s "${BED_METHYL%.bed}.hap1.bed.gz" ]] && [[ -s "${BED_METHYL%.bed}.hap2.bed.gz" ]] ) )
then
	printf 'ERROR: Expected %s.{m,p}at.bed.gz to both exist and be non-empty or %s.hap{1,2}.bed.gz to both exist and be non-empty (plus existing .ok files in either case), but zero or only one existed in one or both pairs.\n' "${BED_METHYL%.bed}" "${BED_METHYL%.bed}" 1>&2
	exit 1
fi

touch "${OK9}"

#printf 'Exiting alignments directory...\n' 1>&2
#printf 'CMD: %s\n' "cd .." 1>&2
cd ..

