library(VariantAnnotation)

vcf_w0 <- readVcf("W0.vcf", genome = "hg38")
rs <- geno(i)$GT %>% rownames()
 inter_df <- data.frame(
  Chr = as.character(seqnames(rowRanges(i))),
  Pos = start(rowRanges(i)),
  REF = as.character(ref(i)),
  ALT = as.character(unlist(alt(i))),
  ID  = as.character(rs),
stringsAsFactors = FALSE
)

ann_raw <- info(vcf_orig_s1)$ANN
df_s1$Ann <- sapply(ann_raw, function(x) paste(x, collapse = ","))

bgzip -c OSV_001_P9.norm_with_ann.vcf > OSV_001_P9.norm_with_ann.vcf.gz
bgzip -c OSV_011_P6.norm_with_ann.vcf > OSV_011_P6.norm_with_ann.vcf.gz

tabix -p vcf OSV_001_P9.norm_with_ann.vcf.gz
tabix -p vcf OSV_011_P6.norm_with_ann.vcf.gz

 bcftools annotate -x INFO --force -o OSV_001_P9_trim.vcf OSV_001_P9.final.vcf
bcftools norm -f GRCh38_renamed.fa -Ov -o OSV_001_P9.norm.vcf OSV_001_P9.noInfo.vcf

bcftools isec -p result OSV_001_P9.norm.vcf.gz OSV_011_P6.norm.vcf.gz


dups <- tmp[duplicated(tmp$ID) | duplicated(tmp$ID, fromLast = TRUE), ]

java -jar snpEff.jar download GRCh38.p14
#ann 
java -Xmx8g -jar snpEff.jar GRCh38.p14 -noStats -noLog result/0000.vcf.gz > result/0000.ann.vcf

# 예시 - dbSNP 주석 붙이기
java -Xmx8g -jar /home/fff0996/WES/snpEff/SnpSift.jar annotate -a -noId -name dbSNP138_ID /home/fff0996/WES/references/GCF_000001405.40.gz 0000.ann.vcf > 0000.ann.dbsnp.vcf
#1000Genome
java -Xmx8g -jar /home/fff0996/WES/snpEff/SnpSift.jar annotate -a -noAlt -noId -name p3_1000G_ /home/fff0996/WES/references/1000Genome.vcf.gz 0000.ann.dbsnp.vcf > 0000.ann.1kg.vcf
 bcftools annotate   -a ../references/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz   -c ID,INFO/AF   -h <(bcftools view -h 0000.ann.dbsnp.vcf | grep "^##")   -O z   -o 0000.ann.dbsnp.1kg.vcf.gz   0000.ann.dbsnp.vcf.gz
                    
# ClinVar
java -Xmx8g -jar /home/fff0996/WES/snpEff/SnpSift.jar annotate -a -noId -name CLINVAR_ /home/fff0996/WES/references/clinvar_20250330.vcf.gz 0000.ann.dbsnp.1kg.vcf > 0000.ann.dbsnp.1kg.clinvar.vcf
#ESP
java -Xmx8g -jar /home/fff0996/WES/snpEff/SnpSift.jar annotate -a -noId -name ESP6500_ /home/fff0996/WES/references/ESP/ESP6500.merged.vcf.gz 0000.ann.dbsnp.1kg.clinvar.vcf > 0000.ann.dbsnp.1kg.clinvar.ESP.vcf
#dbNSFP

# 반복해서 ESP6500, 1000G 등 주석도 붙여주기
##DB 여기서 우선 찾아보기 웬만한건 있어보임 
#https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&invt=AbuMiw
#bgzip result/0000.ann.vcf
#tabix -p vcf result/0000.ann.vcf.gz

#어노테이션 엑셀 만들기
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/ANN\t%INFO/LOF\t%INFO/NMD\t%INFO/p3_1000G_AF\t%INFO/p3_1000G_EUR_AF\t%INFO/p3_1000G_EAS_AF\t%INFO/ESP6500_MAF\t%INFO/CLINVAR_CLNSIG\t%INFO/CLINVAR_CLNDN\t%INFO/CLINVAR_CLNREVSTAT\n' \
-o output.tsv 0000.ann.dbsnp.full.1kg.clinvar.ESP.vcf.gz

echo -e 'CHROM\tPOS\tID\tREF\tALT\tFILTER\tANN\tLOF\tNMD\t1KG_AF\t1KG_EUR_AF\t1KG_EAS_AF\tESP_MAF\tCLNSIG\tCLNDN\tCLNREVSTAT' > header.txt
cat header.txt output.tsv > final_output_with_header.tsv

java -Xmx8g -jar ../snpEff/SnpSift.jar extractFields 0000.ann.dbsnp.1kg.clinvar.ESP.vcf.gz CHROM POS ID REF ALT FILTER CLINVAR_CLNSIG CLINVAR_CLNDISDB CLINVAR_CLNDN CLINVAR_CLNREVSTAT LOF NMD ANN[0].ALLELE ANN[0].EFFECT ANN[0].IMPACT ANN[0].GENE ANN[0].GENEID ANN[0].FEATURE ANN[0].FEATUREID ANN[0].BIOTYPE ANN[0].RANK ANN[0].HGVS_C ANN[0].HGVS_P ANN[0].CDNA_POS ANN[0].CDS_POS ANN[0].AA_POS ANN[0].DISTANCE ANN[0].ERRORS > annotated_variants.tsv
