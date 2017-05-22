library('TCGA2STAT')

rnaseq.ov <- getTCGA(disease="OV", data.type="RNASeq", type="RPKM")
str(rnaseq.ov)
head(rnaseq.ov$dat[,1:3])

# Get the RPKM of genes along with all clinical data, and combine the RPKM with overall survival (OS)
rnaseq_os.ov <- getTCGA(disease="OV", data.type="RNASeq", type="RPKM", clinical=TRUE)
dim(rnaseq_os.ov$merged.dat)
head(rnaseq_os.ov$merged.dat[,1:5])

# Get the RPKM gene expression profiles along with all clinical data, and combine the expression with age
rnaseq_age.ov <- getTCGA(disease="OV", data.type="RNASeq", type="RPKM", clinical=TRUE, cvars="yearstobirth")
head(rnaseq_age.ov$merged.dat[,1:5])

# Get somatic non-silent mutations for ovarian cancer patients
mut.ov <- getTCGA(disease="OV", data.type="Mutation")
# or equivalently
mut.ov <- getTCGA(disease="OV", data.type="Mutation", type="somatic")
str(mut.ov)
head(mut.ov)

# Get SNP array CNA for ovarian cancer patients, no Y chr
cnasnp <- getTCGA(disease="OV", data.type="CNA_SNP")
# or
# Get SNP array CNV for ovarian cancer patients, no Y chr
cnvsnp <- getTCGA(disease="OV", data.type="CNV_SNP")
head(cnvsnp$dat[,1:3])
