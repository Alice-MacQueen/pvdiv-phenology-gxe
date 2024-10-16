import pandas as pd 
import numpy as np 
import sys

newfile = open('pod.greedy.sh','w')
for trait in ['greenup','flowering']:
	for effect in ['random','strong']:
		for cohort in ['gulf','midwest','all']:
				newfile.write('python run.mash.greedy.search.py -p ' + trait + ' -c ' + cohort + ' -e ' + effect + '\n')

newfile = open('pod.bootstrap.sh','w')
for trait in ['greenup','flowering']:
        for effect in ['random','strong']:
                for cohort in ['gulf','midwest','all']:
                        for bootstrap in range(1,101):
                                newfile.write('python run.mash.greedy.bootstrap.py -p ' + trait + ' -c ' + cohort + ' -e ' + effect + ' -b ' + str(bootstrap) + '\n')
