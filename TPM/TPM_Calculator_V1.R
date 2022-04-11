setwd("C:/BIOINF545/FinalProject/v1/DATA")
## V1 
data_files <- list.files("C:/BIOINF545/FinalProject/v1/DATA")
# PolyA
SRR4421939_polyA <- read.table('SRR4421939.count.txt',header = T)
SRR4421940_polyA <- read.table('SRR4421940.count.txt',header = T)
# riboErase
SRR4422432_riboErase <- read.table('SRR4422432.count.txt',header = T)
SRR4422433_riboErase <- read.table('SRR4422433.count.txt',header = T)
#SRR4421939_polyA
idx_1 <- SRR4421939_polyA$Length <50
SRR4421939_polyA <- SRR4421939_polyA[idx_1==FALSE,]
#SRR4421940_polyA
idx_2 <- SRR4421940_polyA$Length <50
SRR4421939_polyA <- SRR4421939_polyA[idx_1==FALSE,]
#SRR4422432_riboErase
idx_3 <- SRR4422432_riboErase$Length <50
SRR4422432_riboErase <- SRR4422432_riboErase[idx_1==FALSE,]
#SRR4422433_riboErase
idx_4 <- SRR4422433_riboErase$Length <50
SRR4422433_riboErase <- SRR4422433_riboErase[idx_1==FALSE,]

# Create a TPM matrix by dividing each column of the counts matrix by estimate of gene length
# SRR4421939_polyA
SRR4421939_polyA_Counts<- SRR4421939_polyA$......alignment.SRR4421939.sam
SRR4421939_polyA_Length <- SRR4421939_polyA$Length
x_SRR4421939_polyA <- SRR4421939_polyA_Counts/SRR4421939_polyA_Length
x_SRR4421939_polyA<-as.matrix(x_SRR4421939_polyA)
tpm_SRR4421939_polyA <- t (t(x_SRR4421939_polyA)*1e6 / colSums(x_SRR4421939_polyA))
SRR4421939_polyA$TPM <- tpm_SRR4421939_polyA
#SRR4421940_polyA
SRR4421940_polyA_Counts<- SRR4421940_polyA$...alignment.SRR4421940.sam
SRR4421940_polyA_Length <- SRR4421940_polyA$Length
x_SRR4421940_polyA <- SRR4421940_polyA_Counts/SRR4421940_polyA_Length
x_SRR4421940_polyA<-as.matrix(x_SRR4421940_polyA)
tpm_SRR4421940_polyA <- t (t(x_SRR4421940_polyA)*1e6 / colSums(x_SRR4421940_polyA))
SRR4421940_polyA$TPM <- tpm_SRR4421940_polyA
#SRR4422432_riboErase
SRR4422432_riboErase_Counts <- SRR4422432_riboErase$...alignment.SRR4422432.sam
SRR4422432_riboErase_Length <- SRR4422432_riboErase$Length
x_SRR4422432_riboErase <- SRR4422432_riboErase_Counts/SRR4422432_riboErase_Length
x_SRR4422432_riboErase <- as.matrix(x_SRR4422432_riboErase)
tpm_SRR4422432_riboErase <- t (t(x_SRR4422432_riboErase)*1e6 / colSums(x_SRR4422432_riboErase))
SRR4422432_riboErase$TPM <- tpm_SRR4422432_riboErase
#SRR4422433_riboErase
SRR4422433_riboErase_Counts <- SRR4422433_riboErase$...alignment.SRR4422433.sam
SRR4422433_riboErase_Length <- SRR4422433_riboErase$Length 
x_SRR4422433_riboErase <- SRR4422433_riboErase_Counts/SRR4422433_riboErase_Length
x_SRR4422433_riboErase <- as.matrix(x_SRR4422433_riboErase)
tpm_SRR4422433_riboErase <- t (t(x_SRR4422433_riboErase)*1e6 / colSums(x_SRR4422433_riboErase))
SRR4422433_riboErase$TPM <- tpm_SRR4422433_riboErase
