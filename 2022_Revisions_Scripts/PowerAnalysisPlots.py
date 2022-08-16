import matplotlib.pyplot as plt
import numpy as np
from math import exp, sqrt
import random
import sys
import statsmodels.api as sm
np.set_printoptions(threshold=sys.maxsize)

#! Main 
def main():
	""" Assign Arguments - ONLY TAKES AS ARGUMENT STR NAME"""

	STR = str(sys.argv[1])

	# Variables and lists
	sample_sizes = list(range(20, 220, 20))
	pve_values = list(np.around(np.arange(0.05, 0.55, 0.05), 3))
	#[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

	eur_power = []
	adx_power = []

	directory = "../Reps_Output/" + STR + "_PowerAnalysis/LD_MAF/"
	out_directory = "../Reps_Output/" + STR + "_PowerAnalysis/PowerPlots/"

	for sample in sample_sizes:
		# Set up file name
		fname = directory + str(sample) + "_PVE_PowerAnalysis.txt"

		eur_power_byPVE = []
		adx_power_byPVE = []
		# Open file and read lines. Save info.
		with open(fname, 'r') as f:
			for line in f:
				if "PVE" in line:
					continue
				else:
					split_line = line.split('\t')
					eur, adx = float(split_line[1]), float(split_line[3])
					eur_power_byPVE.append(eur)
					adx_power_byPVE.append(adx)

		eur_power.append(eur_power_byPVE)
		adx_power.append(adx_power_byPVE)


	# Numpy arrays
	eur_power = np.array(eur_power)
	adx_power = np.array(adx_power)

	plt.clf()

	colors = ['black', 'royalblue', 'navy', 'c', 'grey', 'lightcoral', 'orchid', 'gold', 'indianred', 'violet']

	fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(12,6)) #sharey='row',
	plt.xticks(sample_sizes)
	# plt.yticks(fontsize=13)
	

	# Plot EUR
	for index, PVE in enumerate(pve_values):
		plt.xticks(sample_sizes)
		ax1.set_ylim([0,1])
		plot_eur = eur_power[:,index]
		ax1.plot(sample_sizes, plot_eur, '-o', markersize=4, linewidth=1, label=PVE, color=colors[index])
		ax1.set(title =  STR + " YRI population")
		ax1.set(ylabel="Power")
		ax1.set(xlabel="Individuals sampled")

	# Plot ADMIXED
	for index, PVE in enumerate(pve_values):
		plot_adx = adx_power[:,index]
		ax2.set_ylim([0,1])
		#plt.scatter(sample_sizes, plot_eur, color=colors[index])
		ax2.plot(sample_sizes, plot_adx, '-o', markersize=4, linewidth=1, label=PVE, color=colors[index])
		# Make room on below
		ax2.set(title = STR + " AMR population")
		ax2.set(ylabel="Power")
		ax2.set(xlabel="Individuals sampled")

	fig.subplots_adjust(bottom=0.2)
	fig.legend(pve_values, loc='lower center', bbox_transform=fig.transFigure, title="PVE", ncol=len(pve_values))	
	#fig.suptitle(f'Power Analysis for ' + STR, fontsize=18)
	#fig.supxlabel('Individuals Sampled')
	#fig.supylabel('Power')

	plt.savefig(out_directory + STR + "_Power_TESTER.png", dpi=500)


if __name__ == '__main__':
    main()


