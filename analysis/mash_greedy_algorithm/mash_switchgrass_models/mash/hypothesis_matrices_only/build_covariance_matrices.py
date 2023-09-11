import pandas as pd 
import numpy as np 
import sys
###################
# The following class object is written to generate the covariance matrices hypothesized in the
# mash framework of MacQueen et al. There are three types of matrices included in the analysis
# the first set are the hypothesis based matrices calculated from data collected for individuals
# at each of the eight sites that were analyzed. The second set are those that are data driven 
# and use the PCS of the XXX matrix to generate hypothesized covariance relationships between the 
# different sites. Finally, the canonical matrices are generated using the mash software and represent
# site specific effects, no effects, equal site effects, and a heterogenous effect case.
###################

###################

class switchgrass_covariance(object):

	def __init__(self,pheno):

		#read in the relevant metadata from all the individuals and sites
		self.metadata = pd.read_csv('metadata.txt',sep = "\t")

		#read in sites file that maps four letter site codes to labels used in manuscript
		self.sites = pd.read_csv('sites.txt',sep = '\t')
		self.site_codes = pd.unique(self.sites['manu_site'].sort_values())
		self.SITES = pd.unique(self.sites['SITE'].sort_values())
		self.sitemapper = {i:j for i,j in zip(self.sites['SITE'],self.sites['manu_site'])}
		#Number of folds that were used in the analysis, here 5
		self.folds = [i for i in range(1,6)]

		#Read in all possible ids to be used for generating covariance matrices
		self.allids = pd.read_csv('../gwas/plink_files/Pvirgatum_V5_GWAS_363g_Gulf_and_Midwest_subset_maf0.05.fam', delim_whitespace = True, header = None)
		self.allids = self.allids[1].to_list()
		
		#read in the phenotypes
		self.phenotypes = pd.read_csv(pheno + '_phenos.txt',delim_whitespace = True)
		self.phenotypes = self.phenotypes[(self.phenotypes['SITE'].isin(self.SITES)) & (self.phenotypes['SUBPOP'].isin(['Gulf','Midwest']))].drop_duplicates(keep='first')		
		self.phenotypes['manu_site'] = self.phenotypes['SITE'].map(self.sitemapper)
		
		#read in the list of cohort ids
		self.cohorts = ['all','midwest','gulf']


		###Need to read in the Weather related data (stored in class object as self.weather_data)
		self.weather_data = pd.read_csv('Weather_related_phe_for_asreml_h2_' + pheno + '.txt', sep = '\t')
		self.weather_data['manu_site'] = self.weather_data['SITE'].map(self.sitemapper)

		self.greenupdata = pd.read_csv('Weather_related_phe_for_asreml_h2_greenup.txt', sep = '\t')
		self.greenupdata['manu_site'] = self.weather_data['SITE'].map(self.sitemapper)
		self.greenupphenos = pd.read_csv('greenup_phenos.txt',delim_whitespace = True)
		self.greenupphenos = self.greenupphenos[(self.greenupphenos['SITE'].isin(self.SITES)) & (self.greenupphenos['SUBPOP'].isin(['Gulf','Midwest']))].drop_duplicates(keep='first')		
		self.greenupphenos['manu_site'] = self.greenupphenos['SITE'].map(self.sitemapper)

		#hypothesis driven matrices generating functions
		self.generate_hypothesis_based_cgdd_tave(pheno)
		self.generate_hypothesis_based_GR50(pheno)
		self.generate_hypothesis_based_dyln_change_sec(pheno)
		self.generate_hypothesis_based_gr2fl(pheno)
		self.generate_hypothesis_based_rainfall(pheno)
		self.generate_hypothesis_based_dyln_fl50(pheno)

	def generate_hypothesis_based_cgdd_tave(self,pheno):

		# Enumeration of different covariance matrices, which can be calculated once for 
		# the midwest, gulf, and all populations

		#Cumulative growing days for the n days prior to greenup - variables named cgdd_12c_Xd
			#X=5,10,18

		#Average temperature for the n days prior to greenup - variables named tave_Xd
			#X = 5,10,18
		if pheno == 'greenup':
			for cohort in self.cohorts:
				for matrixtype in ['cgdd_12c_5d','cgdd_12c_10d','cgdd_12c_18d','tave_5d','tave_10d','tave_18d',]:
					covar_inputdf = pd.DataFrame(np.zeros((len(self.allids),len(self.site_codes))), index = self.allids, columns = self.site_codes)
					covar_inputdf = covar_inputdf.reset_index()
					for site in self.site_codes:
						#read in the ids that were used in the training set for this site-fold combination
						site_ids = pd.read_csv('../gwas/wholecohortfiles/' + site + '/' + cohort + '.ids.txt', header = None, sep = '\t')[1]	
						sitephenos = self.phenotypes[(self.phenotypes['manu_site'] == site)]
						
						#take the values for the particular type of matrix that is being generated and set them in the
						#covar_inputdf matrix in the corresponding columns
						sitepheno_dict = {i:j for i,j in zip(sitephenos['PLANT_ID'],sitephenos[matrixtype])}
						covar_inputdf[site] = covar_inputdf['index'].map(sitepheno_dict)
					
					covar_inputdf = covar_inputdf.set_index('index')

					#get the correlation between all observations that are observed in each pair
					corrdf = covar_inputdf.cov()
					# variances = covar_inputdf.var(axis = 0)

					# variances = variances/np.max(variances)
					# np.fill_diagonal(corrdf.values,[i for i in variances])

					# corrdf = np.abs(corrdf)
					corrdf.to_csv('hypothesis_matrices/greenup.' + cohort + '.' + matrixtype + '.U.mat.csv')

			### When the SVD is fixed will updated with the same loop as above for both the Gulf and Midwest populations
		
					
	def generate_hypothesis_based_GR50(self,pheno):
		if pheno == 'flowering':
			for cohort in self.cohorts:
				#Correlations in greenup named 50% greenup
				daylength = self.weather_data[['PLANT_ID','SUBPOP','GR50','manu_site']]

				covar_inputdf = pd.DataFrame(np.zeros((len(self.allids),len(self.site_codes))), index = self.allids, columns = self.site_codes)
				covar_inputdf = covar_inputdf.reset_index()
				for site in self.site_codes:
					#read in the ids that were used in the training set for this site-fold combination
					site_ids = pd.read_csv('../gwas/wholecohortfiles/' + site + '/' + cohort + '.ids.txt', header = None, sep = '\t')[1]	
					sitephenos = self.phenotypes[(self.phenotypes['manu_site'] == site)]
					
					#take the values for the particular type of matrix that is being generated and set them in the
					#covar_inputdf matrix in the corresponding columns
					sitepheno_dict = {i:j for i,j in zip(sitephenos['PLANT_ID'],sitephenos['GR50'])}
					covar_inputdf[site] = covar_inputdf['index'].map(sitepheno_dict)

				covar_inputdf =covar_inputdf.set_index('index')
				#get the correlation between all observations that are observed in each pair
				# corrdf = covar_inputdf.cov()
				corrdf = covar_inputdf.cov()
				# variances = covar_inputdf.var(axis = 0)
				# variances = variances/np.max(variances)
				# np.fill_diagonal(corrdf.values,[i for i in variances])

				corrdf.to_csv('hypothesis_matrices/flowering.' + cohort + '.GR50.U.mat.csv')

	def generate_hypothesis_based_dyln_change_sec(self,pheno):
		#Correlations in daylength change on the day of flowering (daylength_change), but only for flowering
		if pheno == 'flowering':
			for cohort in self.cohorts:
				daylength = self.weather_data[['PLANT_ID','SUBPOP','dyln_change_sec','manu_site']]

				covar_inputdf = pd.DataFrame(np.zeros((len(self.allids),len(self.site_codes))), index = self.allids, columns = self.site_codes)
				covar_inputdf = covar_inputdf.reset_index()
				for site in self.site_codes:
					#read in the ids that were used in the training set for this site-fold combination
					site_ids = pd.read_csv('../gwas/wholecohortfiles/' + site + '/' + cohort + '.ids.txt', header = None, sep = '\t')[1]	
					sitephenos = self.phenotypes[(self.phenotypes['manu_site'] == site) & (self.phenotypes['PLANT_ID'].isin(site_ids))]


					#take the values for the particular type of matrix that is being generated and set them in the
					#covar_inputdf matrix in the corresponding columns
					sitepheno_dict = {i:j for i,j in zip(sitephenos['PLANT_ID'],sitephenos['dyln_change_sec'])}
					covar_inputdf[site] = covar_inputdf['index'].map(sitepheno_dict)

				covar_inputdf =covar_inputdf.set_index('index')
				#get the correlation between all observations that are observed in each pair
				corrdf = covar_inputdf.cov()
				# variances = covar_inputdf.var(axis = 0)
				# variances = variances/np.max(variances)
				# np.fill_diagonal(corrdf.values,[i for i in variances])

				corrdf.to_csv('hypothesis_matrices/flowering.' + cohort + '.dyln_change_sec.U.mat.csv')

	def generate_hypothesis_based_dyln_fl50(self,pheno):
		#Correlations in daylength change on the day of flowering (daylength_change), but only for flowering
		if pheno == 'flowering':
			for cohort in self.cohorts:
				daylength = self.weather_data[['PLANT_ID','SUBPOP','dyln_fl50','manu_site']]

				covar_inputdf = pd.DataFrame(np.zeros((len(self.allids),len(self.site_codes))), index = self.allids, columns = self.site_codes)
				covar_inputdf = covar_inputdf.reset_index()
				for site in self.site_codes:
					#read in the ids that were used in the training set for this site-fold combination
					site_ids = pd.read_csv('../gwas/wholecohortfiles/' + site + '/' + cohort + '.ids.txt', header = None, sep = '\t')[1]	
					sitephenos = self.phenotypes[(self.phenotypes['manu_site'] == site) & (self.phenotypes['PLANT_ID'].isin(site_ids))]
					## & (self.phenotypes['PLANT_ID'].isin(site_ids)


					#take the values for the particular type of matrix that is being generated and set them in the
					#covar_inputdf matrix in the corresponding columns
					sitepheno_dict = {i:j for i,j in zip(sitephenos['PLANT_ID'],sitephenos['dyln_fl50'])}
					covar_inputdf[site] = covar_inputdf['index'].map(sitepheno_dict)

				covar_inputdf =covar_inputdf.set_index('index')
				#get the correlation between all observations that are observed in each pair
				corrdf = covar_inputdf.cov()
				# variances = covar_inputdf.var(axis = 0)
				# variances = variances/np.max(variances)
				# np.fill_diagonal(corrdf.values,[i for i in variances])
				corrdf.to_csv('hypothesis_matrices/flowering.' + cohort + '.dyln_fl50.U.mat.csv')


	def generate_hypothesis_based_gr2fl(self,pheno):
		#Correlations in GDD between greenup and flowering
		if pheno == 'flowering':
			for cohort in self.cohorts:
				daylength = self.weather_data[['PLANT_ID','SUBPOP','cgdd_12c_gr2fl','manu_site']]

				covar_inputdf = pd.DataFrame(np.zeros((len(self.allids),len(self.site_codes))), index = self.allids, columns = self.site_codes)
				covar_inputdf = covar_inputdf.reset_index()
				for site in self.site_codes:
					#read in the ids that were used in the training set for this site-fold combination
					site_ids = pd.read_csv('../gwas/wholecohortfiles/' + site + '/' + cohort + '.ids.txt', header = None, sep = '\t')[1]	
					sitephenos = self.phenotypes[(self.phenotypes['manu_site'] == site) & (self.phenotypes['PLANT_ID'].isin(site_ids))]

					#take the values for the particular type of matrix that is being generated and set them in the
					#covar_inputdf matrix in the corresponding columns
					sitepheno_dict = {i:j for i,j in zip(sitephenos['PLANT_ID'],sitephenos['cgdd_12c_gr2fl'])}
					covar_inputdf[site] = covar_inputdf['index'].map(sitepheno_dict)

				covar_inputdf = covar_inputdf.set_index('index')
				#get the correlation between all observations that are observed in each pair
				corrdf = covar_inputdf.cov()
				# variances = covar_inputdf.var(axis = 0)
				# variances = variances/np.max(variances)
				# np.fill_diagonal(corrdf.values,[i for i in variances])
				
				corrdf.to_csv('hypothesis_matrices/flowering.' + cohort + '.cgdd_12c_gr2fl.U.mat.csv')


	def generate_hypothesis_based_rainfall(self,pheno):

		#Correlations in rainfall in the n days prior to flowering (Rainfalln)
			#1,3,5

		if pheno == 'flowering':
			for cohort in self.cohorts:
				for matrixtype in ['crain_1d','crain_3d','crain_5d']:
					covar_inputdf = pd.DataFrame(np.zeros((len(self.allids),len(self.site_codes))), index = self.allids, columns = self.site_codes)
					covar_inputdf = covar_inputdf.reset_index()
					for site in self.site_codes:
						#read in the ids that were used in the training set for this site-fold combination
						site_ids = pd.read_csv('../gwas/wholecohortfiles/' + site + '/' + cohort + '.ids.txt', header = None, sep = '\t')[1]	
						sitephenos = self.phenotypes[(self.phenotypes['manu_site'] == site)]

						#take the values for the particular type of matrix that is being generated and set them in the
						#covar_inputdf matrix in the corresponding columns
						sitepheno_dict = {i:j for i,j in zip(sitephenos['PLANT_ID'],sitephenos[matrixtype])}
						covar_inputdf[site] = covar_inputdf['index'].map(sitepheno_dict)
					
					covar_inputdf =covar_inputdf.set_index('index')
					#get the correlation between all observations that are observed in each pair
					corrdf = covar_inputdf.cov()
					# variances = covar_inputdf.var(axis = 0)
					# variances = variances/np.max(variances)
					# np.fill_diagonal(corrdf.values,[i for i in variances])

					corrdf.to_csv('hypothesis_matrices/flowering.' + cohort + '.' + matrixtype + '.U.mat.csv')

x = switchgrass_covariance('greenup')
x

y = switchgrass_covariance('flowering')
y