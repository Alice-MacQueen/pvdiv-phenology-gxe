import pandas as pd 
import numpy as np 
import glob
import sys
import os
import argparse
############
#need to add arguments for:
#pheno
#strong/random
#cohort from which effects were derived


parser = argparse.ArgumentParser()
parser.add_argument('-p','--phenotype',type=str,help='phenotype to analyze: greenup or flowering')
parser.add_argument('-e','--effects',type=str,help='type of effects to analyze: random or strong')
parser.add_argument('-c','--cohort',type=str,help='cohort to analyze: gulf, midwest, or all')
parser.add_argument('-b','--bootstrap',type=str,help = 'bootstrap to analyze in this run')

args = parser.parse_args()
pheno = args.phenotype
effects = args.effects
cohort = args.cohort
bootstrap = args.bootstrap

likelihood_path = pd.read_csv('likelihood_paths/' + pheno + '.' + cohort + '.' + effects + '.csv')['matrix'].tolist()
# #test with greenup, all, random
all_hypothesis_matrices = glob.glob('hypothesis_matrices/' + pheno + '.*')
# steps = len(all_hypothesis_matrices)
# mostlikely=[]
# testers = all_hypothesis_matrices

# final_likelihoods = pd.DataFrame(np.zeros((len(all_hypothesis_matrices),2)),columns = ['matrix','likelihood'])
counter = 1

# while len(testers) > 0:
for i,j in enumerate(likelihood_path):
	thispath = likelihood_path[:i+1]
	final = thispath[-1]
	print(final)
	matrices = ['hypothesis_matrices/' + i + '.U.mat.csv' for i in thispath]
	
	outputfile = final
	
	if os.path.isfile('bootstrap_likelihoods/step' + str(counter) + '/' +  outputfile + '.' + cohort + '.' + effects + '.' + str(bootstrap) + '.txt'):
		print('skipped',i)
		pass
	else:
		print('Rscript mash_bootstrap.R -m ' + ','.join(matrices) + ' -p ' + pheno + ' -c ' + cohort + ' -e ' + effects + ' -s ' + str(counter) + ' -b ' + str(bootstrap))
		os.system('Rscript mash_bootstrap.R -m ' + ','.join(matrices) + ' -p ' + pheno + ' -c ' + cohort + ' -e ' + effects + ' -s ' + str(counter) + ' -b ' + str(bootstrap) + '\n')
	counter +=1
# 	likelihood_files = glob.glob('likelihoods/step' + str(counter) + '/' + pheno + '.*.' + cohort + '.' + effects + '.txt')
# 	step_likelihoods = pd.DataFrame(np.zeros((len(testers),2)),columns = ['matrix','likelihood'])

# 	for e,f in enumerate(likelihood_files):
# 		temploglik = [i.strip() for i in open(f)]
# 		step_likelihoods.loc[e] = [f.split('/')[-1],float(temploglik[0])]


	
# 	matrixlabel = step_likelihoods.iloc[step_likelihoods['likelihood'].idxmax()]['matrix'].split('.')
# 	matrixlabel = '.'.join(matrixlabel[:3])
# 	mostlikely.append('hypothesis_matrices/' + matrixlabel + '.U.mat.csv')

# 	final_likelihoods.iloc[counter - 1] = [matrixlabel,step_likelihoods.iloc[step_likelihoods['likelihood'].idxmax()]['likelihood']]
# 	testers.remove('hypothesis_matrices/' + matrixlabel + '.U.mat.csv')
	
# 	print(counter)

# 	counter += 1
        

# final_likelihoods.to_csv('likelihood_paths/' + pheno + '.' + cohort + '.' + effects + '.csv',index = False)

# # likelihoods = [i.strip() for i in likelihood_files]
# # print(likelihoods)
