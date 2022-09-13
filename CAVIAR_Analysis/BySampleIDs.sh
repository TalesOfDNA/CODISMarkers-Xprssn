#!/bin/bash
#SBATCH -n 10
#SBATCH --mem=32G
#SBATCH -t 5:00:00 

module load tabix/0.2.6 
module load bcftools/1.9

input_dir="/users/mbanuel1/data/mbanuel1/GyrmekSNPSTRFiles/Gyrmek_SNP_STR/"
samples="/users/mbanuel1/data/mbanuel1/GyrmekSNPSTRFiles/Sample_IDs/"
temp_dir="/users/mbanuel1/data/mbanuel1/GyrmekSNPSTRFiles/Temp/"
out_dir="/users/mbanuel1/data/mbanuel1/GyrmekSNPSTRFiles/ByPopulation/"

bcftools view -Oz -S $samples/D3S1358_FIN_SampleIDs.txt $input_dir/1kg.snp.str.chr3.vcf.gz > $temp_dir/D3S1358_FIN.vcf.gz
tabix -p vcf $temp_dir/D3S1358_FIN.vcf.gz
tabix $temp_dir/D3S1358_FIN.vcf.gz 3:45482231-45682294 > $temp_dir/D3S1358_FIN_100kb.vcf
awk '!match($3, "STR")' $temp_dir/D3S1358_FIN_100kb.vcf > $out_dir/D3S1358_FIN_100kbVariants.vcf

bcftools view -Oz -S $samples/D3S1358_GBR_SampleIDs.txt $input_dir/1kg.snp.str.chr3.vcf.gz > $temp_dir/D3S1358_GBR.vcf.gz
tabix -p vcf $temp_dir/D3S1358_GBR.vcf.gz
tabix $temp_dir/D3S1358_GBR.vcf.gz 3:45482231-45682294 > $temp_dir/D3S1358_GBR_100kb.vcf
awk '!match($3, "STR")' $temp_dir/D3S1358_GBR_100kb.vcf > $out_dir/D3S1358_GBR_100kbVariants.vcf

bcftools view -Oz -S $samples/D3S1358_TSI_SampleIDs.txt $input_dir/1kg.snp.str.chr3.vcf.gz > $temp_dir/D3S1358_TSI.vcf.gz
tabix -p vcf $temp_dir/D3S1358_TSI.vcf.gz
tabix $temp_dir/D3S1358_TSI.vcf.gz 3:45482231-45682294 > $temp_dir/D3S1358_TSI_100kb.vcf
awk '!match($3, "STR")' $temp_dir/D3S1358_TSI_100kb.vcf > $out_dir/D3S1358_TSI_100kbVariants.vcf

bcftools view -Oz -S $samples/D12S391_CEU_SampleIDs.txt $input_dir/1kg.snp.str.chr12.vcf.gz > $temp_dir/D12S391_CEU.vcf.gz
tabix -p vcf $temp_dir/D12S391_CEU.vcf.gz
tabix $temp_dir/D12S391_CEU.vcf.gz 12:12349954-12550029 > $temp_dir/D12S391_CEU_100kb.vcf
awk '!match($3, "STR")' $temp_dir/D12S391_CEU_100kb.vcf > $out_dir/D12S391_CEU_100kbVariants.vcf

bcftools view -Oz -S $samples/D12S391_GBR_SampleIDs.txt $input_dir/1kg.snp.str.chr12.vcf.gz > $temp_dir/D12S391_GBR.vcf.gz
tabix -p vcf $temp_dir/D12S391_GBR.vcf.gz
tabix $temp_dir/D12S391_GBR.vcf.gz 12:12349954-12550029 > $temp_dir/D12S391_GBR_100kb.vcf
awk '!match($3, "STR")' $temp_dir/D12S391_GBR_100kb.vcf > $out_dir/D12S391_GBR_100kbVariants.vcf

bcftools view -Oz -S $samples/D22S1045_TSI_SampleIDs.txt $input_dir/1kg.snp.str.chr22.vcf.gz > $temp_dir/D22S1045_TSI.vcf.gz
tabix -p vcf $temp_dir/D22S1045_TSI.vcf.gz
tabix $temp_dir/D22S1045_TSI.vcf.gz 22:37436327-37636377 > $temp_dir/D22S1045_TSI_100kb.vcf
awk '!match($3, "STR")' $temp_dir/D22S1045_TSI_100kb.vcf > $out_dir/D22S1045_TSI_100kbVariants.vcf
