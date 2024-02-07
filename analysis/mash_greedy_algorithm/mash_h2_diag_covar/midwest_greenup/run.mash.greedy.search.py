import pandas as pd 
import numpy as np 
import glob
import sys
import os
import argparse
from scipy import stats

############
#need to add arguments for:
#pheno
#strong/random
#cohort from which effects were derived


parser = argparse.ArgumentParser()
parser.add_argument('-p','--phenotype',type=str,help='phenotype to analyze: greenup or flowering')
parser.add_argument('-e','--effects',type=str,help='type of effects to analyze: random or strong')
parser.add_argument('-c','--cohort',type=str,help='cohort to analyze: gulf, midwest, or all')

args = parser.parse_args()
pheno = args.phenotype
effects = args.effects
cohort = args.cohort


#test with greenup, all, random
canonical_singleton = glob.glob('matrices/' + pheno + '.singleton.*')
canonical_identity = glob.glob('matrices/' + pheno + '.identity.*')
canonical_simple_het = glob.glob('matrices/' + pheno + '.simple_het_*.*')
canonical_equal = glob.glob('matrices/' + pheno + '.equal_effects.*')

datadriven = glob.glob('matrices/' + pheno + '.' + cohort + '.' + effects + '.*PCA*')

hypothesis = glob.glob('matrices/' + pheno + '.' + cohort + '*.U.mat.csv')

# matrices/data_driven.flowering.midwest.strong.PCA_2.csv
all_hypothesis_matrices = hypothesis + datadriven + canonical_singleton + canonical_equal + canonical_identity + canonical_simple_het
steps = len(all_hypothesis_matrices)
mostlikely=[]
testers = all_hypothesis_matrices

final_likelihoods = pd.DataFrame(np.zeros((len(all_hypothesis_matrices),2)),columns = ['matrix','likelihood'])
counter = 1

while len(testers) > 0:

	for matrix in testers:

		outputfile = matrix.split('/')[1].split('.')
		outputfile = '.'.join(outputfile[:-1])

		if os.path.isfile('likelihoods/step' + str(counter) + '/' +  outputfile + '.' + cohort + '.' + effects + '.txt'):
			print('skipped')
			pass
		else:
			print('Rscript mash_greedy.R -m ' + matrix + ',' + ','.join(mostlikely) + ' -p ' + pheno + ' -c ' + cohort + ' -e ' + effects + ' -s ' + str(counter))
			os.system('Rscript mash_greedy.R -m ' + matrix + ',' + ','.join(mostlikely) + ' -p ' + pheno + ' -c ' + cohort + ' -e ' + effects + ' -s ' + str(counter) + '\n')

	likelihood_files = glob.glob('likelihoods/step' + str(counter) + '/' + pheno + '.*.' + cohort + '.' + effects + '.txt')
	step_likelihoods = pd.DataFrame(np.zeros((len(testers),2)),columns = ['matrix','likelihood'])

	for e,f in enumerate(likelihood_files):
		temploglik = [i.strip() for i in open(f)]
		step_likelihoods.loc[e] = [f.split('/')[-1],float(temploglik[0])]


	print(step_likelihoods)
	
	step_likelihoods['matrix'] = step_likelihoods['matrix'].str.replace('.' + cohort + '.' + effects + '.txt','.csv')
	print(step_likelihoods)
	matrixlabel = step_likelihoods.iloc[step_likelihoods['likelihood'].idxmax()]['matrix']
	# matrixlabel = '.'.join(matrixlabel[:3])
	mostlikely.append('matrices/' + matrixlabel)

	final_likelihoods.iloc[counter - 1] = [matrixlabel,step_likelihoods.iloc[step_likelihoods['likelihood'].idxmax()]['likelihood']]
	testers.remove('matrices/' + matrixlabel)


	if counter >= 2:
		reduced_ll = final_likelihoods.iloc[counter - 2]['likelihood']
		full_ll = final_likelihoods.iloc[counter - 1]['likelihood']
		LR_statistic = -2*(reduced_ll-full_ll)

		#calculate p-value of test statistic using 1 degrees of freedom
		p_val = stats.chi2.sf(LR_statistic, 1)
		print(counter,p_val)
		if p_val > 0.05:
			break
	
	
	print(counter)

	counter += 1
        
	print(final_likelihoods)
final_likelihoods.to_csv('likelihood_paths/' + pheno + '.' + cohort + '.' + effects + '.stop.condition.csv',index = False)

# likelihoods = [i.strip() for i in likelihood_files]
# print(likelihoods)
