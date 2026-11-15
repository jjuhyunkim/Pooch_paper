suppressPackageStartupMessages(library(DSS))
suppressPackageStartupMessages(library(data.table))
library(bsseq)

cat("working dir:", getwd(), "\n")

args <- commandArgs(trailingOnly = TRUE)
ncores <- 24

# Access individual arguments
chromosome <- as.character(args[1])

print(paste("chromosome:", chromosome))


skip_dict <- list(
    chr12 = c("HG03942"),
    chr19 = c("HG02602")
)

skip_sample = c()

fofn = fread("dss.fofn", header = FALSE)

if ( chromosome == "chrX" ) {
  sampleList = fread("chrX.samples", header = FALSE)
  sampleList = sampleList$V1
  print(paste("num of samples:", length(sampleList)))

}else {
    sampleList= unique(fofn$V1)
    print(paste("num of samples:", length(sampleList)))
}

# sampleList = sampleList[1:10] # debugging test
fofn = fofn[fofn$V1 %in% sampleList,]
print(paste("num of samples:", length(sampleList)))

swap_samples <- c("HG03579", "HG02717", "HG02451")

file_list <- list()
file_name <- character()

skip_chr_samples <- if (chromosome %in% names(skip_dict)) skip_dict[[chromosome]] else character(0)

for (sample in sampleList) {  # Process only the first 10 samples for testing

    if (sample %in% skip_sample) {
        print(paste("Skipping sample:", sample))
        next
    }

    if (sample %in% skip_chr_samples) {
        print(paste("Skipping sample:", sample, "for chromosome:", chromosome))
        next
    }

    file_pat <- paste0("/path/to/pooch/pooch_data/04.methyl_chm13/winnowmap_ont-x-",sample,"-verkkoV2.2.1-asmV4.F0x904.methyl.filt.hap1.chm13.bed")
    file_mat <- paste0("/path/to/pooch/pooch_data/04.methyl_chm13/winnowmap_ont-x-",sample,"-verkkoV2.2.1-asmV4.F0x904.methyl.filt.hap2.chm13.bed")
    
    # file_pat <- paste0("//path/to/pooch/pooch_data/04.methyl_chm13/winnowmap_ont-x-",sample,"-verkkoV2.2.1-asmV4.F0x904.methyl.filt.hap1.chm13.bed")
    # file_mat <- paste0("//path/to/pooch/pooch_data/04.methyl_chm13/winnowmap_ont-x-",sample,"-verkkoV2.2.1-asmV4.F0x904.methyl.filt.hap2.chm13.bed")

    if (!file.exists(file_pat) || !file.exists(file_mat)) {
        print(paste("Missing file(s), skipping sample:", sample))
        next
    }
     
    pat <- fread(cmd = paste0("grep '^", chromosome, "\\b' ", file_pat),header = FALSE,col.names = c("chr", "pos", "end", "N", "X"))[, .(chr, pos, N, X)]
    mat <- fread(cmd = paste0("grep '^", chromosome, "\\b' ", file_mat),header = FALSE,col.names = c("chr", "pos", "end", "N", "X"))[, .(chr, pos, N, X)]
    
    # pat <- fread(cmd = paste0("grep '^", chromosome, "\\b' ", file_pat, "| head -1000000"),header = FALSE,col.names = c("chr", "pos", "end", "N", "X"))[, .(chr, pos, N, X)] # debugging test
    # mat <- fread(cmd = paste0("grep '^", chromosome, "\\b' ", file_mat, "| head -1000000"),header = FALSE,col.names = c("chr", "pos", "end", "N", "X"))[, .(chr, pos, N, X)] # debugging test
    

    if (sample %in% swap_samples) {
        print(paste("Swapping sample:", sample))
        tmp <- pat
        pat <- mat
        mat <- tmp
    }
    print(paste("Processing sample:", sample, "with", nrow(pat), "CpG sites in pat and", nrow(mat), "CpG sites in mat"))
    file_list <- c(file_list, list(mat), list(pat))
    file_name <- c(file_name, paste0(sample, ".mat"), paste0(sample, ".pat"))
}

BSobj_unfilt = makeBSseqData(file_list, file_name)
BSobj_unfilt

# Get coverage matrix (rows = loci, cols = samples)
cov_mat <- getCoverage(BSobj_unfilt, type = "Cov")

# Compute proportion of samples with zero coverage per locus
zero_prop <- rowMeans(cov_mat == 0)

# Keep loci where <= 5% of samples have zero coverage
keep <- zero_prop <= 0.05

n_keep = sum(zero_prop <= 0.05)
n_lost = sum(zero_prop > 0.05)
write.table(data.frame(chromosome, n_keep, n_lost), paste0("dmls.N05_filt.",chromosome,".nCpG.csv"), row.names = FALSE, quote = FALSE, sep = "\t")

# Subset BSseq object
BSobj <- BSobj_unfilt[keep, ] # debugging test

n_row =  dim(BSobj@colData)[1]/2
haplotype = factor(rep(c("mat","pat"), n_row))
pair = factor( rep(1:n_row, each=2) )
design = data.frame(haplotype, pair)
design

BSobj

DMLfit = DMLfit.multiFactor(BSobj, design, formula = ~ haplotype + pair)
dmlTest = DMLtest.multiFactor(DMLfit, term="haplotype")

# dmlTest = DMLtest(BSobj, group1=file_name[grep("mat", file_name)], group2=file_name[grep("pat", file_name)], smoothing=TRUE, ncores=ncores-1,equal.disp=FALSE,smoothing.span=500)
# dmlTest.sm <- DMLtest(BSobj, group1=matv, group2=patv, smoothing=TRUE, equal.disp=FALSE, smoothing.span=500, ncores=ncores-1)

# dmls = callDML(dmlTest, delta=0, p.threshold=1)
write.table(dmlTest, paste0("dmls.N05_filt.",chromosome,".all.csv"), row.names = FALSE, quote = FALSE, sep = "\t")

# identify DMRs
## delta = if non-zero, posterior prob that diff in means is >delta
## p.threshold = p-val threshold for including loci in the DMR
## minlen = the minimum length of a DMR
## minCG = the minimum number of CpG sites needed for a DMR
## dis.merge = if two DMRs are within this distance of each other, they will be merged
## pct.sig = the minimum percentage of CpG sites that must be significant in a region

dmrs_5e2 = callDMR(dmlTest, p.threshold=0.05, delta=0, minlen=50, minCG=3, dis.merge=100, pct.sig=0.5)
dmrs_1e3 = callDMR(dmlTest, p.threshold=0.001, delta=0, minlen=50, minCG=3, dis.merge=100, pct.sig=0.5)
dmrs_1e5 = callDMR(dmlTest, p.threshold=1e-5, delta=0, minlen=50, minCG=3, dis.merge=100, pct.sig=0.5)

head(dmrs_5e2)
head(dmrs_1e3)
head(dmrs_1e5)

write.table(dmrs_5e2, paste0("dmrs.N05_filt.",chromosome,".5e2.csv"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(dmrs_1e3, paste0("dmrs.N05_filt.",chromosome,".1e3.csv"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(dmrs_1e5, paste0("dmrs.N05_filt.",chromosome,".1e5.csv"), row.names = FALSE, quote = FALSE, sep = "\t")

sessionInfo()