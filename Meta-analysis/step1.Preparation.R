mvp_edge_gwas_bmi_eur <- read.delim("mvp_edge_gwas_bmi_alpha_eur.txt")
ukb_edge_gwas_bmi_eur <- read.delim("ukb_edge_gwas_bmi_alpha_eur.txt")

mvp_edge_gwas_bmi_eur <- tidyr::separate(mvp_edge_gwas_bmi_eur, Variant.ID, into = c("chr", "pos", "ref", "alt"), sep = "[:]")
mvp_edge_gwas_bmi_eur$pos <- as.numeric(mvp_edge_gwas_bmi_eur$pos)
mvp_edge_gwas_bmi_eur$chr <- as.numeric(gsub("\\D", "", mvp_edge_gwas_bmi_eur$chr))

ukb_edge_gwas_bmi_eur <- tidyr::separate(ukb_edge_gwas_bmi_eur, Variant.ID, into = c("chr", "pos", "ref", "alt"), sep = "[:]")
ukb_edge_gwas_bmi_eur$pos <- as.numeric(ukb_edge_gwas_bmi_eur$pos)
ukb_edge_gwas_bmi_eur$chr <- as.numeric(ukb_edge_gwas_bmi_eur$chr)

mvp_edge_gwas_bmi_eur$SNP <- paste0(mvp_edge_gwas_bmi_eur$chr,":",mvp_edge_gwas_bmi_eur$pos,":",mvp_edge_gwas_bmi_eur$ref,":",mvp_edge_gwas_bmi_eur$alt)
mvp_edge_gwas_bmi_eur <- mvp_edge_gwas_bmi_eur[,c(15,9,10,12,1,2,7,6)]
colnames(mvp_edge_gwas_bmi_eur) <- c("SNP","BETA","SE","P","CHR","BP","A1","A2")
#mvp_edge_gwas_bmi_eur$OR <- exp(mvp_edge_gwas_bmi_eur$OR)

ukb_edge_gwas_bmi_eur$SNP <- paste0(ukb_edge_gwas_bmi_eur$chr,":",ukb_edge_gwas_bmi_eur$pos,":",ukb_edge_gwas_bmi_eur$ref,":",ukb_edge_gwas_bmi_eur$alt)
ukb_edge_gwas_bmi_eur <- ukb_edge_gwas_bmi_eur[,c(16,10,11,13,1,2,7,6)]
colnames(ukb_edge_gwas_bmi_eur) <- c("SNP","BETA","SE","P","CHR","BP","A1","A2")
#ukb_edge_gwas_bmi_eur$OR <- exp(ukb_edge_gwas_bmi_eur$OR)

write.table(mvp_edge_gwas_bmi_eur,"~/plink/mvp_edge_gwas_bmi_eur_meta.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ukb_edge_gwas_bmi_eur,"~/plink/ukb_edge_gwas_bmi_eur_meta.txt", quote = FALSE, row.names = FALSE, sep = "\t")

###
mvp_edge_gwas_cad_eur <- read.delim("mvp_edge_gwas_cad_alpha_eur.txt")
ukb_edge_gwas_cad_eur <- read.delim("ukb_edge_gwas_cad_alpha_eur.txt")

mvp_edge_gwas_cad_eur <- tidyr::separate(mvp_edge_gwas_cad_eur, Variant.ID, into = c("chr", "pos", "ref", "alt"), sep = "[:]")
mvp_edge_gwas_cad_eur$pos <- as.numeric(mvp_edge_gwas_cad_eur$pos)
mvp_edge_gwas_cad_eur$chr <- as.numeric(gsub("\\D", "", mvp_edge_gwas_cad_eur$chr))

ukb_edge_gwas_cad_eur <- tidyr::separate(ukb_edge_gwas_cad_eur, Variant.ID, into = c("chr", "pos", "ref", "alt"), sep = "[:]")
ukb_edge_gwas_cad_eur$pos <- as.numeric(ukb_edge_gwas_cad_eur$pos)
ukb_edge_gwas_cad_eur$chr <- as.numeric(ukb_edge_gwas_cad_eur$chr)

mvp_edge_gwas_cad_eur$SNP <- paste0(mvp_edge_gwas_cad_eur$chr,":",mvp_edge_gwas_cad_eur$pos,":",mvp_edge_gwas_cad_eur$ref,":",mvp_edge_gwas_cad_eur$alt)
mvp_edge_gwas_cad_eur <- mvp_edge_gwas_cad_eur[,c(15,9,10,12,1,2,7,6)]
colnames(mvp_edge_gwas_cad_eur) <- c("SNP","OR","SE","P","CHR","BP","A1","A2")
mvp_edge_gwas_cad_eur$OR <- exp(mvp_edge_gwas_cad_eur$OR)

ukb_edge_gwas_cad_eur$SNP <- paste0(ukb_edge_gwas_cad_eur$chr,":",ukb_edge_gwas_cad_eur$pos,":",ukb_edge_gwas_cad_eur$ref,":",ukb_edge_gwas_cad_eur$alt)
ukb_edge_gwas_cad_eur <- ukb_edge_gwas_cad_eur[,c(16,10,11,13,1,2,7,6)]
colnames(ukb_edge_gwas_cad_eur) <- c("SNP","OR","SE","P","CHR","BP","A1","A2")
ukb_edge_gwas_cad_eur$OR <- exp(ukb_edge_gwas_cad_eur$OR)

write.table(mvp_edge_gwas_cad_eur,"~/plink/mvp_edge_gwas_cad_eur_meta.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ukb_edge_gwas_cad_eur,"~/plink/ukb_edge_gwas_cad_eur_meta.txt", quote = FALSE, row.names = FALSE, sep = "\t")


