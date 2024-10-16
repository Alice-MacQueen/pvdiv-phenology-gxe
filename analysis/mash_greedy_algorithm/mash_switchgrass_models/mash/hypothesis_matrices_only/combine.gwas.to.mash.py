import pandas as pd 
import numpy as np 
import sys


####Create the true matrices of effects for all three populations

populations = ['gulf','midwest','all']
sites = ['MI','MO','NE','OK','SD','TX1','TX2','TX3']
traits = ['flowering','greenup']

for trait in traits:
	for pop in populations:
		for effects in ['strong','random']:
			print(trait,pop,effects)
			first = pd.read_csv('/stor/work/Kirkpatrick/sps2972/switchgrass/gwas/gwasoutputs/' + sites[0] + '/' + trait + '.wc.' + pop + '.' + effects + '.results.csv')
			betas = first[['estim']].rename(columns={'estim':sites[0]})
			stderr = first[['std.err']].rename(columns={'std.err':sites[0]})
			for site in sites[1:]:
				temp = pd.read_csv('/stor/work/Kirkpatrick/sps2972/switchgrass/gwas/gwasoutputs/' + site + '/' + trait + '.wc.' + pop + '.' + effects + '.results.csv')
				betas[site] = temp[['estim']]
				stderr[site] = temp[['std.err']]


			betas = betas.dropna()
			stderr = stderr.dropna()
			if betas.shape[0] != stderr.shape[0]:

				idx = betas.index.intersection(stderr.index)
				betas = betas.loc[idx]
				stderr = stderr.loc[idx]
			print(betas.shape)
			print(stderr.shape)
			betas.to_csv('combined_mash_effects/' + trait + '.' + pop + '.betas.' + effects + '.txt',sep = '\t',header = False, index = False)
			stderr.to_csv('combined_mash_stderrs/' + trait + '.' + pop + '.stderrs.' + effects + '.txt',sep = '\t',header = False, index = False)
			
