import msprime
import numpy as np
from math import exp, sqrt

def copy_number_matrix(strs):
    """
	From Tskit documentation
    Returns the copy number matrix from the specified tree sequence
    simulated under a MicrosatMutationModel.
    """
    C = np.zeros((strs.num_sites, strs.num_samples), dtype=int)
    for var in strs.variants():
        alleles = np.array([int(allele) for allele in var.alleles])
        C[var.site.id] = alleles[var.genotypes]
    return C


def simulate_phenotype(snps, strs, s_len):
	"""
	Takes STRs and SNPs mutations from ts output in simulation;
	selects STR and 4 causal SNPs in LD of at least 30% to STR
	simulates gene expression under additive model
	"""
	str_list = []
	snp_list = []

	# Get STR positions
	for str_ in strs.variants():
		str_list.append(str_.site.position)

	for snp_ in snps.variants():
		snp_list.append(snp_.genotypes)
	
	# Select STR index by min distance from middle of seq
	diff = abs(np.asarray(str_list) - int(s_len/2))
	STR_pos = np.argmin(diff)

	# Function for getting actual genotypes
	C = copy_number_matrix(strs)

	CSF1PO1 = C[STR_pos]

	# Subset by SNPs that have at least 30% correlation with STR 
	cor = []
	for row in range(0, len(snp_list)):
		cor.append(abs(np.corrcoef(snp_list[row], CSF1PO1)[0,1]))

	pool_causal = [i for i, v in enumerate(cor) if(v > .30)]
	causal_index = random.choices(pool_causal, k=2)
	causal_snps = [snp_list[i] for i in causal_index]

	# Transposing to get data in ROWS = INDVIDUALS and COLS = SNPs
	# No Header or row names in matrix for phenotype simulation
	causal_snps1 = np.transpose(causal_snps)

	""" 
	Convert STR and SNPs into diploid by adding two adjacent haplotypes
	"""
	# We add adjacent rows for diploid genotypes
	causal_snps = (causal_snps1[::2] + causal_snps1[1::2])
	# print("CAusal after mean")
	# print(causal_snps)
	
	CSF1PO = []
	for ind in range(0, len(CSF1PO1), 2):
		CSF1PO.append((CSF1PO1[ind] + CSF1PO1[ind+1])/2)
	# print("CSF1PO after mean")
	# print(CSF1PO)
	
	# Get dimensions of matrix to set up for simulating an additive phenotype model
	# Standardize the values so that the variance of genetic and environmental components sum to 1
	ind, nsnps = causal_snps.shape[0], causal_snps.shape[1]
	Xmean = np.mean(causal_snps, axis=0)
	Xsd = np.std(causal_snps, axis=0, ddof=1)
	Xmarginal = (causal_snps - Xmean)/Xsd

	# Sample betas from normal distribution mean=0, std=1
	betas = np.random.normal(0, 1, size=nsnps)
	#betas = np.asarray([1.3839820, -0.8688081, 0.2746070, -0.3381371])
	pve=0.6
	### Simulate marginal (additive) effects ###
	y_marginal = np.matmul(Xmarginal, betas)
	betas = betas * sqrt(pve/np.var(y_marginal, ddof=1))
	y_marginal= np.matmul(Xmarginal, betas)

	""" NOTE: Each effect size for the marginal, epistatic, and random error effects are drawn from a standard normal distribution. Meaning beta ~ MVN(0,I) and epsilon ~ MVN(0,I). We then scale both the additive genetic effects so that collectively they explain a fixed proportion of genetic variance. Namely, the additive effects make up the narrow-sense heritability (or h^2). Once we obtain the final effect sizes for all causal SNPs, we draw errors to achieve the target 1-h^2.
	"""

	### Simulate residual error drawing noise from normal distribution ###
	y_err = np.random.normal(0, 1, size=ind)
	y_err = y_err * sqrt((1-pve)/np.var(y_err, ddof=1))

	#! Simulate gene expression 
	# Simulate the continuous phenotypes
	y = y_marginal + y_err
	
	return [np.asarray(CSF1PO), y]
