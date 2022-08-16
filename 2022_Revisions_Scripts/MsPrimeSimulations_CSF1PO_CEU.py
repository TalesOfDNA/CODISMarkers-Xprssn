import msprime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#from IPython.display import SVG
import random
import statsmodels.api as sm
from math import exp, sqrt
import sys
# import demes
# import demesdraw
from ModelDemography import demographic_model_browning
import time
np.set_printoptions(threshold=sys.maxsize)
#from StrSnpHandling import save_genotypes_byPop
# start_time = time.time()
# end_time = time.time() - start_time

#! Helper Functions
""" Parameters for call_simis """
model = msprime.SMM(lo=22, hi=28)
STR_mu = str(0.001569517) #CSF1PO
SNP_mu = str(2e-8)


def run_simis(num_samples, ancestry_seed, mutation_seed, s_len, rep):
	"""
	Runs simulations with given number of samples and ancestry seeds 

	RETURNS:
		Generators for simulated SNPs and STR genotypes
	"""
	
	demographic_model = demographic_model_browning()

	# Call tree (ancestry) simulation using the demographic model
	ts = msprime.sim_ancestry(samples={"ADMIX": num_samples, "EUR": num_samples}, demography=demographic_model, sequence_length=s_len, recombination_rate=2e-8)

	#! For sampling YRI instead later
	# ts = msprime.sim_ancestry(samples={"ADMIX": num_samples, "AFR": num_samples}, demography=demographic_model, sequence_length=s_len, recombination_rate=2e-8)

	# Add mutations to the tree 
	snps = msprime.sim_mutations(ts, rate=SNP_mu, random_seed=ancestry_seed)
	strs = msprime.sim_mutations(ts, rate=STR_mu, model=model, random_seed=mutation_seed)

	return strs, snps

def save_genotypes_byPop(strs, snps, num_samples, s_len, rep):
	"""
	Function takes in STRs and SNPs, combines adjacent haplotypes into diploid genotypes,
	saves output in txt files by population

	Input:
		  treeseq variants for:
		  	 STRs
			 SNPs
	Output: 
		  txt files
		  in rows are samples, cols SNPs
	"""	

	STR_name = "CSF1PO"
	output_dir = "../Reps_Output/CSF1PO_Corrected/"

	#! Fetch STR at middle position and its genotypes, combine adjacent haplotypes into single genotype per diploid sample
	from itertools import islice
	str_pos = int(s_len/2)
	str_info = next(islice(strs.variants(), str_pos, str_pos+1))
	alleles = np.array([int(allele) for allele in str_info.alleles])
	CSF1PO_1 = alleles[str_info.genotypes]

	# We add adjacent rows for diploid genotypes (taking the average of the STR)
	CSF1PO = (CSF1PO_1[::2] + CSF1PO_1[1::2])/2

	# Subset Admixed (first samples) and European Samples 
	CSF1PO_Admx = CSF1PO[:num_samples]
	CSF1PO_EUR = CSF1PO[num_samples:]

	#! Fetch SNPs, combine adjacent haplotypes into diploid genotypes
	# Here SNPs are initially as rows and samples as columns so we transpose the matrix
	# to combine adjacent haploid samples into diploid samples
	snp_list = []
	for snp_ in snps.variants():
		snp_list.append(snp_.genotypes)

	# Transpose to get samples are rows, snps as columns
	snp_list = np.array(snp_list).T

	# We add adjacent rows for diploid genotypes
	snps_diploid = (snp_list[::2] + snp_list[1::2])

	# Subset Admixed and EUR samples
	snps_Admx = snps_diploid[:num_samples]
	snps_EUR = snps_diploid[num_samples:]

	# Prepare file names
	adx_str_fname = output_dir + str(rep) + "_" + str(num_samples) + "_adx_str.txt"
	adx_snp_fname = output_dir + str(rep) + "_" + str(num_samples) + "_adx_snps.txt"
	eur_str_fname = output_dir + str(rep) + "_" + str(num_samples) + "_eur_str.txt"
	eur_snp_fname = output_dir + str(rep) + "_" + str(num_samples) + "_eur_snps.txt"

	# SAMPLES as Rows, SNPs as columns
	np.savetxt(adx_str_fname, CSF1PO_Admx, fmt='%1.2f')
	np.savetxt(adx_snp_fname, snps_Admx, fmt='%1d')
	np.savetxt(eur_str_fname, CSF1PO_EUR, fmt='%1.2f')
	np.savetxt(eur_snp_fname, snps_EUR, fmt='%1d')

def main():

	""" Call helper function to run simulation reps and save outputs """

	# Parameters
	s_len = 200001

	""" Assign Arguments """
	sample_size = int(sys.argv[1])
	rep = int(sys.argv[2])

	num_replicates = 1

	# Set random seeds
	rng = np.random.RandomState(42)
	seeds = rng.randint(1, 2**31, size=(num_replicates, 2))

	#! Run Simis

	strs, snps = run_simis(sample_size, *seeds[0], s_len, rep)
	save_genotypes_byPop(strs, snps, sample_size, s_len, rep)

if __name__ == '__main__':
    main()