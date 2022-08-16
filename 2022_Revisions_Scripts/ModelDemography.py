#! Model of American Admixture from Browning et al 2011
import math
import msprime
from math import exp

def demographic_model_browning():
	
	T_OOA = 920
	demography = msprime.Demography()
	demography.add_population(name="AFR", description="African", initial_size=14474)
	demography.add_population(
		name="EUR",
		description="European",
		initial_size=34039,
		growth_rate=0.0038,
	)
	demography.add_population(
		name="EAS",
		description="East Asian",
		initial_size=45852,
		growth_rate=0.0048,
	)
	demography.add_population(
		name="ADMIX",
		description="Admixed America",
		initial_size=54664,
		growth_rate=0.05,
	)
	demography.add_population(
		name="OOA",
		description="Bottleneck out-of-Africa",
		initial_size=1861,
	)
	demography.add_population(
		name="AMH", description="Anatomically modern humans", initial_size=14474
	)
	demography.add_population(
		name="ANC",
		description="Ancestral equilibrium",
		initial_size=7310,
	)
	demography.set_symmetric_migration_rate(["AFR", "EUR"], 2.5e-5)
	demography.set_symmetric_migration_rate(["AFR", "EAS"], 0.78e-5)
	demography.set_symmetric_migration_rate(["EUR", "EAS"], 3.11e-5)

	demography.add_admixture(
		12,
		derived="ADMIX",
		ancestral=["AFR", "EUR", "EAS"],
		proportions=[1 / 6, 2 / 6, 3 / 6],
	);

	demography.add_population_split(T_OOA, derived=["EUR", "EAS"], ancestral="OOA")
	demography.add_symmetric_migration_rate_change(
		time=T_OOA, populations=["AFR", "OOA"], rate=15e-5
	)
	demography.add_population_split(2040, derived=["OOA", "AFR"], ancestral="AMH")
	demography.add_population_split(5920, derived=["AMH"], ancestral="ANC")

	return demography

def demographic_model1(sample_size):

	mu=1.25e-8 # mutation rate per bp
	rho=1e-8 # recombination rate per bp
	nbp = 1e8 # generate 100 Mb
	N0=7310 # initial population size
	Thum=5920 # time (gens) of advent of modern humans 
	Naf=14474 # size of african population
	Tooa=2040 # number of generations back to Out of Africa 
	Nb=1861 # size of out of Africa population
	mafb=1.5e-4 # migration rate Africa and Out-of-Africa 
	Teu=920 # number generations back to Asia-Europe split 
	Neu=1032; Nas=554 # bottleneck population sizes 
	mafeu=2.5e-5; mafas=7.8e-6; meuas=3.11e-5 # mig. rates 
	reu=0.0038 # growth rate per generation in Europe 
	ras=0.0048 # growth rate per generation in Asia 
	Tadmix=12 # time of admixture
	Nadmix=30000 # initial size of admixed population 
	radmix=.05 # growth rate of admixed population

	# pop0 is Africa, pop1 is Europe, pop2 is Asia, pop3 is admixed
	refsamplesize = 0
	admsamplesize = sample_size

	demography = msprime.Demography()
	pop_config = [
	msprime.PopulationConfiguration(sample_size=refsamplesize,initial_size=Naf,growth_rate=0.0), msprime.PopulationConfiguration(sample_size=refsamplesize,initial_size=Neu*exp(reu*Teu),growth_rate=reu), msprime.PopulationConfiguration(sample_size=refsamplesize,initial_size=Nas*exp(ras*Teu),growth_rate=ras), msprime.PopulationConfiguration(sample_size=admsamplesize,initial_size=Nadmix*exp(radmix*Tadmix),growth_rate=radmix)] 
	mig_mat = [[0,mafeu,mafas,0],[mafeu,0,meuas,0], [mafas,meuas,0,0],[0,0,0,0]]
	# Admixture event, 1/6 Africa, 2/6 Europe, 3/6 Asia
	admixture_event = [ msprime.MassMigration(time=Tadmix,source=3,destination=0,proportion=1.0/6.0), msprime.MassMigration(time=Tadmix+0.0001,source=3,destination=1,proportion=2.0/5.0), msprime.MassMigration(time=Tadmix+0.0002,source=3,destination=2,proportion=1.0)]

	# Asia and Europe split
	eu_event = [
	msprime.MigrationRateChange(time=Teu,rate=0.0), msprime.MassMigration(time=Teu+0.0001,source=2,destination=1,proportion=1.0), msprime.PopulationParametersChange(time=Teu+0.0002,initial_size=Nb,growth_rate=0.0,population_id=1), msprime.MigrationRateChange(time=Teu+0.0003,rate=mafb,matrix_index=(0,1)), msprime.MigrationRateChange(time=Teu+0.0003,rate=mafb,matrix_index=(1,0))]
	# Out of Africa event
	ooa_event = [
	msprime.MigrationRateChange(time=Tooa,rate=0.0), msprime.MassMigration(time=Tooa+0.0001,source=1,destination=0,proportion=1.0)]

	# initial population size
	init_event = [ msprime.PopulationParametersChange(time=Thum,initial_size=N0,population_id=0)]

	events = admixture_event + eu_event + ooa_event + init_event

	return demography