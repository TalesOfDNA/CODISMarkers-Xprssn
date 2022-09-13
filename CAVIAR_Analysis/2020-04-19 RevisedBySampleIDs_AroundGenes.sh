#!/bin/bash
#SBATCH -n 10
#SBATCH --mem=32G
#SBATCH -t 5:00:00 

module load tabix/0.2.6 
module load bcftools/1.9

input_dir="/users/mbanuel1/data/mbanuel1/GyrmekSNPSTRFiles/Gyrmek_SNP_STR/"
samples="/users/mbanuel1/data/mbanuel1/GyrmekSNPSTRFiles/Sample_IDs/"
temp_dir="/users/mbanuel1/data/mbanuel1/GyrmekSNPSTRFiles/Temp/"
out_dir="/users/mbanuel1/data/mbanuel1/GyrmekSNPSTRFiles/100kbFromGenes/"

bcftools view -Oz -S $samples/D3S1358_FIN_SampleIDs.txt $input_dir/1kg.snp.str.chr3.vcf.gz > $temp_dir/D3S1358_FIN.vcf.gz
tabix -p vcf $temp_dir/D3S1358_FIN.vcf.gz
tabix $temp_dir/D3S1358_FIN.vcf.gz 3:45329998-45690913  > $temp_dir/LARS2_FIN_100kb.vcf
#awk '!match($3, "STR")' $temp_dir/LARS2_FIN_100kb.vcf > $out_dir/LARS2_FIN_100kbVariants.vcf

bcftools view -Oz -S $samples/D3S1358_GBR_SampleIDs.txt $input_dir/1kg.snp.str.chr3.vcf.gz > $temp_dir/D3S1358_GBR.vcf.gz
tabix -p vcf $temp_dir/D3S1358_GBR.vcf.gz
tabix $temp_dir/D3S1358_GBR.vcf.gz 3:45329998-45690913  > $temp_dir/LARS2_GBR_100kb.vcf
#awk '!match($3, "STR")' $temp_dir/LARS2_GBR_100kb.vcf > $out_dir/LARS2_GBR_100kbVariants.vcf

bcftools view -Oz -S $samples/D3S1358_TSI_SampleIDs.txt $input_dir/1kg.snp.str.chr3.vcf.gz > $temp_dir/D3S1358_TSI.vcf.gz
tabix -p vcf $temp_dir/D3S1358_TSI.vcf.gz
tabix $temp_dir/D3S1358_TSI.vcf.gz 3:45329998-45690913  > $temp_dir/LARS2_TSI_100kb.vcf
#awk '!match($3, "STR")' $temp_dir/LARS2_TSI_100kb.vcf > $out_dir/LARS2_TSI_100kbVariants.vcf

bcftools view -Oz -S $samples/D12S391_CEU_SampleIDs.txt $input_dir/1kg.snp.str.chr12.vcf.gz > $temp_dir/D12S391_CEU.vcf.gz
tabix -p vcf $temp_dir/D12S391_CEU.vcf.gz
tabix $temp_dir/D12S391_CEU.vcf.gz 12:12408342-12610001 > $temp_dir/LOH2CR2_CEU_100kb.vcf
#awk '!match($3, "STR")' $temp_dir/LOH2CR2_CEU_100kb.vcf > $out_dir/LOH2CR2_CEU_100kbVariants.vcf

bcftools view -Oz -S $samples/D12S391_GBR_SampleIDs.txt $input_dir/1kg.snp.str.chr12.vcf.gz > $temp_dir/D12S391_GBR.vcf.gz
tabix -p vcf $temp_dir/D12S391_GBR.vcf.gz
tabix $temp_dir/D12S391_GBR.vcf.gz 12:12408342-12610001 > $temp_dir/LOH2CR2_GBR_100kb.vcf
#awk '!match($3, "STR")' $temp_dir/LOH2CR2_GBR_100kb.vcf > $out_dir/LOH2CR2_GBR_100kbVariants.vcf

bcftools view -Oz -S $samples/D12S391_CEU_SampleIDs.txt $input_dir/1kg.snp.str.chr12.vcf.gz > $temp_dir/D12S391_CEU.vcf.gz
tabix -p vcf $temp_dir/D12S391_CEU.vcf.gz
tabix $temp_dir/D12S391_CEU.vcf.gz 12:12410013-12719840  > $temp_dir/LOH2CR1_CEU_100kb.vcf
#awk '!match($3, "STR")' $temp_dir/LOH2CR1_CEU_100kb.vcf > $out_dir/LOH2CR1_CEU_100kbVariants.vcf

bcftools view -Oz -S $samples/D22S1045_TSI_SampleIDs.txt $input_dir/1kg.snp.str.chr22.copy.vcf.gz > $temp_dir/D22S1045_TSI.vcf.gz
tabix -p vcf $temp_dir/D22S1045_TSI.vcf.gz
tabix $temp_dir/D22S1045_TSI.vcf.gz 22:37421878-37671094 > $temp_dir/IL2RB_TSI_100kb.vcf
#awk '!match($3, "STR")' $temp_dir/IL2RB_TSI_100kb.vcf > $out_dir/IL2RB_TSI_100kbVariants.vcf
