#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Saturday August 26, 2020

Run CAVIAR on Significant SNPs and STRs

@author: mayrabanuelos
"""

import sys, re, os
from subprocess import Popen, PIPE, DEVNULL

extension = '.txt'
zfiles = list(f for f in os.listdir("ZScores/") if f.endswith('Zscores' + extension))
z_dir = "ZScores/"
ld_dir = "LD_Matrices/"
out_dir = "CAVIAR_output/"

for f in zfiles:
# Set up file names
	f_name = f.split('_')[1]
	population = f.split('_')[2]
	#population = temp.split('.')[0]

	## Pull LD file
	LD_file = ld_dir + f_name + "_" + population + "_LD.txt"
	out_file = out_dir + f_name + "_" + population + "_4"
	z_file = z_dir + f
	# Run CAVIAR
	cmd = "CAVIAR -o %s -l %s -z %s -c 4 -f 1"%(out_file, LD_file, z_file)
	p = Popen(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
	output = p.communicate()[0]
#posteriors = list(f for f in os.listdir("../Output/CAVIAR_output/") if f.endswith('_post' + extension))
#posteriors = list(f for f in os.listdir("../Output/CAVIAR_output/") if f.endswith('6_post'))

# for f in posteriors:
#	# Sort posterior probabilities, Caviar Score and Retain top 25
#	post_file = "CAVIAR_output/" + f
#	top = out_dir + f + "_sorted"
#	cmd = "sort -k3 -n -r %s | awk 'NR>=0&&NR<=25' > %s"%(post_file, top)
#	p = Popen(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)