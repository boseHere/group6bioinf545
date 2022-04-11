### Processed Data ######## V4 
setwd("C:/BIOINF545/FinalProject/counts/Data/")
## Load the data 
# PolyA
SRR4421939_polyA_V4 <- read.table('SRR4421939.count.txt',header = T)
SRR4421940_polyA_v4 <- read.table('SRR4421940.count.txt',header = T)
# riboErase
SRR4422432_riboErase_v4 <- read.table('SRR4422432.count.txt',header = T)
SRR4422433_riboErase_v4 <- read.table('SRR4422433.count.txt',header = T)
# Throw out anything less than 50 bp 
#SRR4421939_polyA
idx_1 <- SRR4421939_polyA_V4$Length <50
SRR4421939_polyA_V4 <- SRR4421939_polyA_V4[idx_1==FALSE,]
#SRR4421940_polyA
idx_2 <- SRR4421940_polyA_v4$Length <50
SRR4421940_polyA_v4 <- SRR4421940_polyA_v4[idx_2==FALSE,]
#SRR4422432_riboErase
idx_3 <- SRR4422432_riboErase_v4$Length <50
SRR4422432_riboErase_v4 <- SRR4422432_riboErase_v4[idx_3==FALSE,]
#SRR4422433_riboErase
idx_4 <- SRR4422433_riboErase_v4$Length <50
SRR4422433_riboErase_v4 <- SRR4422433_riboErase_v4[idx_4==FALSE,]

# Create a TPM matrix by dividing each column of the counts matrix by estimate of gene length
# SRR4421939_polyA
SRR4421939_polyA_v4_Counts<- SRR4421939_polyA_V4$...alignment.SRR4421939.sam
SRR4421939_polyA_v4_Length <- SRR4421939_polyA_V4$Length
x_SRR4421939_polyA <- SRR4421939_polyA_v4_Counts/SRR4421939_polyA_v4_Length
x_SRR4421939_polyA<-as.matrix(x_SRR4421939_polyA)
tpm_SRR4421939_polyA <- t (t(x_SRR4421939_polyA)*1e6 / colSums(x_SRR4421939_polyA))
SRR4421939_polyA_V4$TPM <- tpm_SRR4421939_polyA

#SRR4421940_polyA
SRR4421940_polyA_v4_Counts<- SRR4421940_polyA_v4$...alignment.SRR4421940.sam
SRR4421940_polyA_v4_Length <- SRR4421940_polyA_v4$Length
x_SRR4421940_polyA <- SRR4421940_polyA_v4_Counts/SRR4421940_polyA_v4_Length
x_SRR4421940_polyA<-as.matrix(x_SRR4421940_polyA)
tpm_SRR4421940_polyA <- t (t(x_SRR4421940_polyA)*1e6 / colSums(x_SRR4421940_polyA))
SRR4421940_polyA_v4$TPM <- tpm_SRR4421940_polyA

#SRR4422432_riboErase
SRR4422432_riboErase_v4_Counts <- SRR4422432_riboErase_v4$...alignment.SRR4422432.sam
SRR4422432_riboErase_v4_Length <- SRR4422432_riboErase_v4$Length
x_SRR4422432_riboErase <- SRR4422432_riboErase_v4_Counts/SRR4422432_riboErase_v4_Length
x_SRR4422432_riboErase <- as.matrix(x_SRR4422432_riboErase)
tpm_SRR4422432_riboErase <- t (t(x_SRR4422432_riboErase)*1e6 / colSums(x_SRR4422432_riboErase))
SRR4422432_riboErase_v4$TPM <- tpm_SRR4422432_riboErase

#SRR4422433_riboErase
SRR4422433_riboErase_v4_Counts <- SRR4422433_riboErase_v4$...alignment.SRR4422433.sam
SRR4422433_riboErase_v4_Length <- SRR4422433_riboErase_v4$Length 
x_SRR4422433_riboErase <- SRR4422433_riboErase_v4_Counts/SRR4422433_riboErase_v4_Length
x_SRR4422433_riboErase <- as.matrix(x_SRR4422433_riboErase)
tpm_SRR4422433_riboErase <- t (t(x_SRR4422433_riboErase)*1e6 / colSums(x_SRR4422433_riboErase))
SRR4422433_riboErase_v4$TPM <- tpm_SRR4422433_riboErase


###### Write out the outputs ######
write.csv(SRR4421939_polyA_V4, 'SRR4421939_polyA_V4.csv')
write.csv(SRR4421940_polyA_v4,'SRR4421940_polyA_v4.csv')
write.csv(SRR4422432_riboErase_v4,'SRR4422432_riboErase_v4.csv')
write.csv(SRR4422433_riboErase_v4,'SRR4422433_riboErase_v4.csv')