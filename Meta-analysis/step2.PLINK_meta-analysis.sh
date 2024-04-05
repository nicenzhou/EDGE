./plink --meta-analysis ukb_edge_gwas_cad_eur_meta.txt mvp_edge_gwas_cad_eur_meta.txt --out meta_analysis_edge_cad_eur;
./plink --meta-analysis ukb_edge_gwas_bmi_eur_meta.txt mvp_edge_gwas_bmi_eur_meta.txt + qt --out meta_analysis_edge_bmi_eur;

./plink --meta-analysis ukb_plink_gwas_cad_eur_meta.txt mvp_plink_gwas_cad_eur_meta.txt --out meta_analysis_plink_cad_eur;
./plink --meta-analysis ukb_plink_gwas_bmi_eur_meta.txt mvp_plink_gwas_bmi_eur_meta.txt + qt --out meta_analysis_plink_bmi_eur;
