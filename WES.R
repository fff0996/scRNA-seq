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
java -Xmx8g -jar snpEff.jar GRCh38.p14 result/0000.vcf.gz > result/0000.ann.vcf

# 예시 - dbSNP 주석 붙이기
java -Xmx8g -jar SnpSift.jar annotate -noId -info dbSNP138_ID /path/to/GCF_000001405.40.gz result/0000.ann.vcf > result/0000.ann.dbsnp.vcf
#1000Genome
java -Xmx8g -jar ~/WES/snpEff/SnpSift.jar annotate -noId -name p3_1000G_ /path/to/1000Genome.vcf.gz result/0000.ann.dbsnp.vcf > result/0000.ann.1kg.vcf
# ClinVar
java -Xmx8g -jar SnpSift.jar annotate -noId -name CLINVAR_ /path/to/clinvar_20240716.vcf.gz result/0000.ann.dbsnp.vcf > result/0000.ann.clinvar.vcf

#ESP
java -Xmx8g -jar SnpSift.jar annotate -noId -name ESP6500_ ../references/ESP/ESP6500.merged.vcf.gz ../result/0000.ann.dbsnp.1kg.clinvar.vcf > ../result/0000.ann.dbsnp.1kg.clinvar.ESP.vcf
#dbNSFP

# 반복해서 ESP6500, 1000G 등 주석도 붙여주기
##DB 여기서 우선 찾아보기 웬만한건 있어보임 
https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&invt=AbuMiw
bgzip result/0000.ann.vcf
tabix -p vcf result/0000.ann.vcf.gz
