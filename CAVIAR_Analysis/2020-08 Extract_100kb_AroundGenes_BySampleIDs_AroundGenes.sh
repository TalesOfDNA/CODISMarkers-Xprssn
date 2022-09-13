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
tabix $temp_dir/D3S1358_FIN.vcf.gz 3:45329998-45690913  > $out_dir/LARS2_FIN_100kb.vcf

bcftools view -Oz -S $samples/D3S1358_GBR_SampleIDs.txt $input_dir/1kg.snp.str.chr3.vcf.gz > $temp_dir/D3S1358_GBR.vcf.gz
tabix -p vcf $temp_dir/D3S1358_GBR.vcf.gz
tabix $temp_dir/D3S1358_GBR.vcf.gz 3:45329998-45690913  > $out_dir/LARS2_GBR_100kb.vcf

bcftools view -Oz -S $samples/D3S1358_TSI_SampleIDs.txt $input_dir/1kg.snp.str.chr3.vcf.gz > $temp_dir/D3S1358_TSI.vcf.gz
tabix -p vcf $temp_dir/D3S1358_TSI.vcf.gz
tabix $temp_dir/D3S1358_TSI.vcf.gz 3:45329998-45690913  > $out_dir/LARS2_TSI_100kb.vcf

bcftools view -Oz -S $samples/CSF1PO_FIN_SampleIDs.txt $input_dir/1kg.snp.str.chr5.vcf.gz > $temp_dir/CSF1PO_FIN.vcf.gz
tabix -p vcf $temp_dir/CSF1PO_FIN.vcf.gz
tabix $temp_dir/CSF1PO_FIN.vcf.gz  5:149332854-149592935 > $out_dir/CSF1R_FIN_100kb.vcf

bcftools view -Oz -S $samples/CSF1PO_GBR_SampleIDs.txt $input_dir/1kg.snp.str.chr5.vcf.gz > $temp_dir/CSF1PO_GBR.vcf.gz
tabix -p vcf $temp_dir/CSF1PO_GBR.vcf.gz
tabix $temp_dir/CSF1PO_GBR.vcf.gz  5:149332854-149592935 > $out_dir/CSF1R_GBR_100kb.vcf

bcftools view -Oz -S $samples/D18S51_YRI_SampleIDs.txt $input_dir/1kg.snp.str.chr18.vcf.gz > $temp_dir/D18S51_YRI.vcf.gz
tabix -p vcf $temp_dir/D18S51_YRI.vcf.gz
tabix $temp_dir/D18S51_YRI.vcf.gz  18:60894959-61134743  > $out_dir/KDSR_YRI_100kb.vcf