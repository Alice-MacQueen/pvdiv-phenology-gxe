import pandas as pd 
import numpy as np 
import sys

newfile = open('pod.setup.canonical.data.driven.sh','w')
for trait in ['greenup','flowering']:
	for effect in ['random','strong']:
		for cohort in ['gulf','midwest','all']:
				newfile.write('Rscript mash_make_data_driven_mats.R -p ' + trait + ' -c ' + cohort + ' -e ' + effect + '\n')

for trait in ['greenup','flowering']:
        newfile.write('Rscript mash_make_canonical_mats.R -p ' + trait + '\n')

newfile = open('pod.greedy.sh','w')
for trait in ['greenup','flowering']:
        for effect in ['random','strong']:
                for cohort in ['gulf','midwest','all']:
                                newfile.write('python run.mash.greedy.search.py -p ' + trait + ' -c ' + cohort + ' -e ' + effect + '\n')
