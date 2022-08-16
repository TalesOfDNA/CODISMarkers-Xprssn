import matplotlib.pyplot as plt
import numpy as np
from math import exp, sqrt
import random
import sys
import statsmodels.api as sm
np.set_printoptions(threshold=sys.maxsize)

def regression_analysis(str_genos, exp):
	""" Performs Regression_Analysis 
		for STR genotypes (X) and simulated gene expression (Y)
		and returns p-value and correlation coefficient
	"""
	mod = sm.OLS(exp, str_genos)
	fii = mod.fit()
	p_value = fii.summary2().tables[1]['P>|t|'] 
	r = fii.rsquared

	return p_value, r

def simulate_phenotype(adx_str_fname, adx_snps_fname, eur_str_fname, eur_snps_fname, n_causal):
	"""
	Takes STRs and SNPs genotypes from simulation
	select 4 causal SNPs in LD of at least 30% with STR
	simulates gene expression under additive model
	fits regression model

	Input:
		For EUR and ADX populations:
			str 1D array 
			snps 2D array (SAMPLES as Rows, SNPs as columns)
			pve (SNP effect size)

	Output:
		boolean (1 significant, 0 not significant)
	"""

	# Load files
	# SNPs files have SAMPLES as rows, SNPs as columns
	adx_str = np.loadtxt(adx_str_fname)
	adx_snps = np.loadtxt(adx_snps_fname)
	eur_str = np.loadtxt(eur_str_fname)
	eur_snps = np.loadtxt(eur_snps_fname)

	# Remove SNPs where there is no variation within populations (cols with same value)
	monomorphic_sites = np.all(adx_snps == adx_snps[0], axis=0)
	adx_snps = adx_snps[:, ~monomorphic_sites]

	monomorphic_sites = np.all(eur_snps == eur_snps[0], axis=0)
	eur_snps = eur_snps[:, ~monomorphic_sites]

	out_directory = "../Reps_Output/CSF1PO_PowerAnalysis/MAF_0.05_Histograms/"

	# Calculate Minor allele frequency (MAF) ADX
	alt_allele_freq = np.sum(adx_snps, axis = 0)/(2*adx_snps.shape[0]) 
	MAF_adx = []
	for freq in alt_allele_freq:
		if freq < 0.5:
			MAF_adx.append(freq)
		else:
			MAF_adx.append(1-freq)
	
	MAF_adx = np.array(MAF_adx)
	# Checking if any SNPs satisfy MAF threshold
	MAF_within_adx = np.any(MAF_adx > 0.05)
	if MAF_within_adx:
		maf_exists = 1
	else:
		maf_exists = 0

	# Calculate Minor allele frequency (MAF) EUR
	alt_allele_freq = np.sum(eur_snps, axis = 0)/(2*eur_snps.shape[0])
	MAF_eur = []
	for freq in alt_allele_freq:
		if freq < 0.5:
			MAF_eur.append(freq)
		else:
			MAF_eur.append(1-freq)
	
	MAF_eur = np.array(MAF_eur)
	# Checking if any SNPs satisfy MAF threshold
	MAF_within_eur = np.any(MAF_eur > 0.05)
	if MAF_within_eur:
		maf_exists_eur = 1
	else:
		maf_exists_eur = 0

	# print(np.sum(MAF > 0.05))
	# # print("MAF")
	# # print(MAF)

	# # Subset by SNPs that have at least 30% correlation with STR 
	# cor = []
	# for col in range(0, adx_snps.shape[1]):
	# 	cor.append(abs(np.corrcoef(adx_snps[:, col], adx_str))[0,1])

	# cor_eur = []
	# for col in range(0, eur_snps.shape[1]):
	# 	cor_eur.append(abs(np.corrcoef(eur_snps[:, col], eur_str))[0,1])

	# # Make numpy arrays 
	# cor = np.array(cor)
	# cor_eur = np.array(cor_eur)

	# # Subset SNPs with cor > 0.30 with STR
	# mask = np.where(cor > 0.30)
	# mask_eur = np.where(cor_eur > 0.30)

	# #! In case we don't have any SNPs with LD > 0.30
	# #! Choose Causal Randomly from full set
	# if len(mask[0]) == 0:
	# 	indices = range(0, adx_snps.shape[1])
	# 	causal_index_adx = random.choices(indices, k=4)
	# else:
	# 	causal_index_adx = random.choices(mask[0], k=4)
	
	# # EUR
	# if len(mask_eur[0]) == 0:
	# 	indices = range(0, eur_snps.shape[1])
	# 	causal_index_eur = random.choices(indices, k=4)
	# else:
	# 	causal_index_eur = random.choices(mask_eur[0], k=4)

	# # Choose 4 causal variants
	# #causal_index = random.choices(mask[0], k=4)
	# causal_snps = adx_snps[:, causal_index_adx]

	# #causal_index = random.choices(mask_eur[0], k=4)
	# causal_snps_eur = eur_snps[:, causal_index_eur]

	# # Get dimensions of matrix to set up for simulating an additive phenotype model
	# # Standardize the values so that the variance of genetic and environmental components sum to 1
	# ind, nsnps = causal_snps.shape[0], causal_snps.shape[1]
	# Xmean = np.mean(causal_snps, axis=0)
	# Xsd = np.std(causal_snps, axis=0, ddof=1)
	# Xmarginal = (causal_snps - Xmean)/Xsd

	# ind_eur, nsnps_eur = causal_snps_eur.shape[0], causal_snps_eur.shape[1]
	# Xmean_eur = np.mean(causal_snps_eur, axis=0)
	# Xsd_eur = np.std(causal_snps_eur, axis=0, ddof=1)
	# Xmarginal_eur = (causal_snps_eur - Xmean_eur)/Xsd_eur

	# # Sample betas from normal distribution mean=0, std=1
	# betas = np.random.normal(0, 1, size=nsnps)
	# betas_eur = np.random.normal(0, 1, size=nsnps_eur)

	# ### Simulate marginal (additive) effects ###
	# y_marginal = np.matmul(Xmarginal, betas)
	# betas = betas * sqrt(pve/np.var(y_marginal, ddof=1))
	# y_marginal= np.matmul(Xmarginal, betas)

	# y_marginal_eur = np.matmul(Xmarginal_eur, betas_eur)
	# betas_eur = betas_eur * sqrt(pve/np.var(y_marginal_eur, ddof=1))
	# y_marginal_eur = np.matmul(Xmarginal_eur, betas_eur)

	# """ NOTE: Each effect size for the marginal, epistatic, and random error effects are drawn from a standard normal distribution. Meaning beta ~ MVN(0,I) and epsilon ~ MVN(0,I). We then scale both the additive genetic effects so that collectively they explain a fixed proportion of genetic variance. Namely, the additive effects make up the narrow-sense heritability (or h^2). Once we obtain the final effect sizes for all causal SNPs, we draw errors to achieve the target 1-h^2.
	# """

	# ### Simulate residual error drawing noise from normal distribution ###
	# y_err = np.random.normal(0, 1, size=ind)
	# y_err = y_err * sqrt((1-pve)/np.var(y_err, ddof=1))

	# y_err_eur = np.random.normal(0, 1, size=ind_eur)
	# y_err_eur = y_err_eur * sqrt((1-pve)/np.var(y_err_eur, ddof=1))

	# #! Simulate gene expression 
	# # Simulate the continuous phenotypes
	# y_adx = y_marginal + y_err
	# y_eur = y_marginal_eur + y_err_eur

	# #! Regression standardizing STR and without standardizing
	# # Regression of EUR samples
	# Xmean = np.mean(eur_str)
	# Xsd = np.std(eur_str, ddof=1)
	# eur_norm = (eur_str - Xmean)/Xsd
	# pval_eur_norm, rsqrd = regression_analysis(eur_norm, y_eur)
	# pval_eur, rsqrd = regression_analysis(eur_str, y_eur)

	# # Regression of ADX samples
	# Xmean = np.mean(adx_str)
	# Xsd = np.std(adx_str, ddof=1)
	# adx_norm = (adx_str - Xmean)/Xsd
	# pval_adx_norm, rsqrd = regression_analysis(adx_norm, y_adx)
	# pval_adx, rsqrd = regression_analysis(adx_str, y_adx)
	
	# #! Check for significance in association
	# if pval_eur_norm.to_list()[0] < 0.05:
	# 	eur_norm = 1
	# else:
	# 	eur_norm = 0

	# if pval_eur.to_list()[0] < 0.05:
	# 	eur = 1
	# else:
	# 	eur = 0
	
	# if pval_adx_norm.to_list()[0] < 0.05:
	# 	adx_norm = 1
	# else:
	# 	adx_norm = 0
	
	# if pval_adx.to_list()[0] < 0.05:
	# 	adx = 1
	# else:
	# 	adx = 0
	
	return MAF_adx, maf_exists, MAF_eur, maf_exists_eur

#! Main 
def main():
	""" Assign Arguments """
	sample_size = int(sys.argv[1])
	STR = str(sys.argv[2])

	# Reps and snp effect size
	reps = list(range(1, 1001))
	pve_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
	n_causal = 4

	directory = "../Reps_Output/" + STR + "_Corrected/"

	out_directory = "../Reps_Output/" + STR + "_PowerAnalysis/MAF_0.05_Histograms/"

	file_name = out_directory + str(sample_size) + "_PVE_PowerAnalysis.txt"

	#f = open(file_name, 'w')
	#f.write("PVE\tNormEUR\tEUR\tNormADX\tADX\n")

	# # power analysis with diff pve values
	# for pve in pve_values:
	# 	eurNorm_allReps = []
	# 	adxNorm_allReps = []
	# 	eur_allReps = []
	# 	adx_allReps = []
	MAF_adx_allReps = []
	MAF_eur_allReps = []
	maf_exists_allReps = []
	maf_eur_exists_allReps = []
	for rep in reps:
		# Set up file names
		adx_str_fname = directory + str(rep) + "_" + str(sample_size) + "_" + "adx_str.txt"
		adx_snps_fname = directory + str(rep) + "_" + str(sample_size) + "_" + "adx_snps.txt"
		eur_str_fname = directory + str(rep) + "_" + str(sample_size) + "_" + "eur_str.txt"
		eur_snps_fname = directory + str(rep) + "_" + str(sample_size) + "_" + "eur_snps.txt"

		# Simulate phenotype
		#eur_norm, eur, adx_norm, adx, MAF = simulate_phenotype(adx_str_fname, adx_snps_fname, eur_str_fname, eur_snps_fname, pve, n_causal)
		MAF_adx, maf_exists, MAF_eur, maf_exists_eur = simulate_phenotype(adx_str_fname, adx_snps_fname, eur_str_fname, eur_snps_fname, n_causal)
		# eurNorm_allReps.append(eur_norm)
		# adxNorm_allReps.append(adx_norm)
		# eur_allReps.append(eur)
		# adx_allReps.append(adx)
		MAF_adx_allReps.append(MAF_adx)
		maf_exists_allReps.append(maf_exists)
		MAF_eur_allReps.append(MAF_eur)
		maf_eur_exists_allReps.append(maf_exists_eur)

	# power_eur = sum(eur_allReps)/1000
	# power_adx = sum(adx_allReps)/1000
	# power_eur_norm = sum(eurNorm_allReps)/1000
	# power_adx_norm = sum(adxNorm_allReps)/1000

	# Make single list from array of arrays
	adx_allMAF = np.concatenate(MAF_adx_allReps).ravel()
	eur_allMAF = np.concatenate(MAF_eur_allReps).ravel()
	print("Num Reps where No SNPs meet MAF threshold in ADX")
	print(1000 - sum(maf_exists_allReps))

	print("Num Reps where No SNPs meet MAF threshold in EUR")
	print(1000 - sum(maf_eur_exists_allReps))

	plt.hist(adx_allMAF)
	title = "Histogram of Minor Allele Frequencies (MAF) \nAdmixed pop | 200kb window\nSample size " + str(sample_size) + " diploid individuals"
	plt.title(title)
	plt.xlabel("MAF")
	plt.xlim([0,1])
	plt.tight_layout()  # This line avoids title being cut off
	figure_name = "Histogram_MAF_" + str(sample_size) + "_ADX.png"
	plt.savefig(out_directory + figure_name)
	plt.clf()

	plt.hist(eur_allMAF)
	title = "Histogram of Minor Allele Frequencies (MAF) \nEuropean pop | 200kb window \nSample size " + str(sample_size) + " diploid individuals"
	plt.title(title)
	plt.xlabel("MAF")
	plt.xlim([0,1])
	plt.tight_layout() # This line avoids title being cut off
	figure_name = "Histogram_MAF_" + str(sample_size) + "_EUR.png"
	plt.savefig(out_directory + figure_name)
	plt.clf()

	# f.write(str(pve) + "\t" + str(power_eur_norm) + "\t" + 
	# 				str(power_eur) + "\t" +  str(power_adx_norm) + 
	# 				"\t" + str(power_adx) + "\n")

	#f.close()

if __name__ == '__main__':
    main()


