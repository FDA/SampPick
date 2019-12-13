'''
** SampPick Sample Picking Algorithm ***

This code is available at https://www.github.com/fda/SampPick

This software runs a simulated annealing algorithm on a user supplied collection of samples. 
The software will iteratively subsitute sample individuals to create an sample cohort of size N (useer specified). 
The cohorts will be analyzed by comparing the frequency distribution of the alleles as compared to a background distribution of interest. 
We have supplied a background frequency distribution for alleles found in the North American population. 
In addition, a reference sample could be added with the command -r so that you can add a premade cohort to compare
 your optimized cohort against .


The sofware will ouput three files. The first is a bar chart of alleles and their frequencies in the background distribution,
 the optimized sample( and optionally the reference) The second is a csv file of the two or three frequency distributions 
 described above. The third is the list of samples chosen for the optimized samples. This list will present the sample and
  their respective DRB alleles

Additionally, there is a second mode accessible through the flag --basic. It will allow for the simple output of a bar chart
 of a target frequency and a reference population as well as the associated frequencies.

Please cite : 

McGill JR, Yogurtcu ON, Verthelyi D, Yang H, Sauna ZE (2019) SampPick: Selection of a Cohort of Subjects Matching a Population HLA Distribution. Front. Immunol. 10:2894. doi: 10.3389/fimmu.2019.02894


'''



from random import sample
from random import random
from collections import Counter
import math
import matplotlib.pyplot as plt
import pandas as pd
import argparse

#These arguements are passed from the command line running the script with -h willl bring up a help screen
parser =argparse.ArgumentParser( formatter_class=argparse.ArgumentDefaultsHelpFormatter,usage='\n%(prog)s [options] available',description='A tool for the analysis and optimization of Sample Distributions',epilog='Please see the README file for information about parameter values')
parser.add_argument('-s','--size',metavar='',default=50,type=int,help='Sample size desired (Integer)')
parser.add_argument('-r','--reference',metavar='',default='',type=str, help='CSV file to compare with (see README for help)')
parser.add_argument('-i','--iterations',metavar='',type=int,default=10000,help='Number of iterations for Simulated Annealing (Integer)')
parser.add_argument('-a','--alpha',metavar='',type=float,default=0.0007,help='Temperature decrease Alpha value (Float)')
parser.add_argument('-c','-changes',metavar='',default=1,type=int,help='Number of changes made on each iteration of the algorithm (Integer)')
parser.add_argument('-t','--target',metavar='',type=str,default='NA_DRB1_frequencies.csv',help='Filename of csv file (see README for help)')
parser.add_argument('-o','--outfile', metavar='',type=str,default='Results',help='Name for the csv and Png files created')
parser.add_argument('-n','--name',metavar='',type=str,default='Optimized Sample',help='Name for the Reference Sample')
parser.add_argument('available',type=str,help='Filename of csv file see Readme for help')
parser.add_argument('--basic',action='store_true', help='just analyze the available sample\nin its entirety')
parser.add_argument('--verbose',action='store_true',default=False, help='print output every 1000 iterations')

pas=parser.parse_args()



'''
This function calculates the Kullback-Leibler Distance between two probability distributions
The inputs are intended to be dictionaries with the keys being the name of an allele (can be a different feature)
and the value being a frequency.

It is important to avoid items in either dictionary with a value of 0 to avoid errors stemming from either division by 0 or 
taking the log of 0
'''
def rel_entropy(sample_guess2,background_distribution2):
		re=[]
		for a in background_distribution2:
			if sample_guess2[a]!=0:
				re.append(sample_guess2[a]*math.log(sample_guess2[a]/background_distribution2[a]))
		return(sum(re))

'''
This function calculates an average frequency between the background and sample.
It then takes the average of the distances between the background and average and sample and average.
The lambda shown is changeable but is generally set to 0.5
'''

def get_jensen_score_sample(sample_guess,background_distribution,lam=0.5):
	
	average_dist={key:(sample_guess[key]+background_distribution[key])/2 for key in background_distribution}
	
	js_score=lam*rel_entropy(sample_guess,average_dist)+(1-lam)*rel_entropy(background_distribution,average_dist)
	jsd=math.sqrt(js_score)
	return(jsd)

'''
Basic mode just takes the available sample and graphs the distribution of alleles
'''
if pas.basic:
	'''
	inputs:
	target_distribution - a csv file with alleles in column one and frequencies in column 2. It is expected to have a header
	samples_available - a csv file with a header and three columns. A unique sample identifier, allele 1, allele 2
	outfilename - a name for the output file. 'Results' will be used by default
	'''

	def simulated_annealing(target_distribution=pas.target,samples_available=pas.available,outfilename=pas.outfile):
		#Create a dictionary for the target population allele frequencies
		target_distribution={l.split(',')[0]:float(l.split(',')[1].strip()) for l in open(target_distribution).readlines()[1:] if float(l.split(',')[1].strip())!=0}
		
		#These are mainly used to keep the keys in order since the dictionary is unordered 
		dfkeys=target_distribution.keys()
		df_target_values=[target_distribution[k] for k in dfkeys]
		
		#This pandas data frame will be used for storing all frequencies
		#If a reference Sample is used, a three column data frame will be created 
		if len(pas.reference):
			drb_df=pd.DataFrame(index=dfkeys,columns=['Background Distribution',pas.name,'Reference'])
			#ref lines  flattens a list of lists per subject to a simple list of alleles in the population available.
			ref_lines=[i.strip().split(',')[1:3] for i in open(pas.reference).readlines()[1:]]
			ref_lines=[i for k in ref_lines for i in k]
			
			#A Counter is used to gather frequencies of alleles
			
			ref_distribution=Counter(ref_lines)
			#It is scaled so that the counts are actually frequencies
			for allele in ref_distribution:
				ref_distribution[allele]/=float(len(ref_lines))
			#Then it is turned into a dictionary
			ref_distribution={allele:ref_distribution[allele] for allele in ref_distribution}
			#Add zeros for alleles not found in the subjects but found in the target Population
			#These are needed for getting average values later
			for allele in target_distribution:
				if allele not in ref_distribution:
					ref_distribution.update({allele:0})
			#Add reference distribution values to data frame
			drb_df['Reference']=[ref_distribution[k] for k in dfkeys]


		# If no reference population is used, the data frame only needs two columns
		else:
			drb_df=pd.DataFrame(index=dfkeys,columns=['Background Distribution',pas.name])
		
		#ref lines  flattens a list of lists per subject to a simple list of alleles in the population available.

		ref_lines=[i.strip().split(',')[1:3] for i in open(samples_available).readlines()[1:]]
		ref_lines=[i for k in ref_lines for i in k]
		
		#A Counter is used to gather frequencies of alleles
		ref_distribution=Counter(ref_lines)
		

		#Add zeros for alleles not found in the subjects but found in the target population
		#These are needed for getting average values later
		for allele in target_distribution:
			if allele not in ref_distribution:
				ref_distribution.update({allele:0})
		
		#Then it is turned into a dictionary
		ref_distribution={allele:ref_distribution[allele] for allele in ref_distribution}
		
		
		#It is scaled so that the counts are actually frequencies
		for allele in ref_distribution:
			ref_distribution[allele]/=float(len(ref_lines))

		#Add available samples and background values dictionary to data frame		
		drb_df[pas.name]=[ref_distribution[k] for k in dfkeys]
		
		drb_df['Background Distribution']=df_target_values
		drb_df.to_csv(outfilename+'_frequencies.csv')
		# Bar chart of frequencies. It takes out alleles which have a frequency of less than 0.01
		# This will eliminate hundreds of alleles that are not relevant to see. It should get down to 
		# 30 or so alleles for DRB1 alleles

		ax=drb_df[drb_df.sum(axis=1)>.01].sort_index(axis=0).plot(kind='bar',cmap='Accent',width=.85)
		ax.set_ylabel('Frequency in Population')
		plt.tight_layout()
		fig=ax.get_figure()
		
		# Saved using the name from the inputs
		fig.savefig(outfilename+'.png',height=1200,width=800,dpi=600)
		print("Finished")
	
		
else:

	'''
	This is the full version of the algorithm. It will create a N person sample of the available samples that is "closest"
	to a reference population with respect to allele frequencies

	inputs:
	N - the size of the sample you would like to draw from the available sample
	target_distribution - a csv file with alleles in column one and frequencies in column 2. It is expected to have a header
	samples_available - a csv file with a header and three columns. A unique sample identifier, allele 1, allele 2
	changes - the number of subjects which will be changed at  each iteration of the algorithm. 
			- one should be good in most cases but any number up to N should work
	number_of iteration - see README
	temperature_decrease_alpha - see README
	outfilename - a name for the output file. 'Results' will be used by default
	'''	
	def simulated_annealing(N=pas.size,target_distribution=pas.target,samples_available=pas.available,changes=pas.c,number_of_iterations=pas.iterations,temperature_decrease_alpha=pas.alpha,outfilename=pas.outfile):
		
		starting_temperature=1.0
		scores=[]

		#Create a dictionary for the target population allele frequencies
		
		target_distribution={l.split(',')[0]:float(l.split(',')[1].strip()) for l in open(target_distribution).readlines()[1:] if float(l.split(',')[1].strip())!=0}
		
		#These are mainly used to keep the keys in order since the dictionary is unordered 
		dfkeys=target_distribution.keys()
		df_target_values=[target_distribution[k] for k in dfkeys]
		
		#This pandas data frame will be used for storing all frequencies
		#If a reference Sample is used, a three column data frame will be created 
		if len(pas.reference):
			drb_df=pd.DataFrame(index=dfkeys,columns=['Background Distribution',pas.name,'Reference'])
			#ref lines  flattens a list of lists per subject to a simple list of alleles in the population available.
			ref_lines=[i.strip().split(',')[1:3] for i in open(pas.reference).readlines()[1:]]
			ref_lines=[i for k in ref_lines for i in k]
			
			#A Counter is used to gather frequencies of alleles
			
			ref_distribution=Counter(ref_lines)
			#It is scaled so that the counts are actually frequencies
			for allele in ref_distribution:
				ref_distribution[allele]/=float(len(ref_lines))
			#Then it is turned into a dictionary
			ref_distribution={allele:ref_distribution[allele] for allele in ref_distribution}
			#Add zeros for alleles not found in the subjects but found in the target Population
			#These are needed for getting average values later
			for allele in target_distribution:
				if allele not in ref_distribution:
					ref_distribution.update({allele:0})
			#Add reference distribution values to data frame
			drb_df['Reference']=[ref_distribution[k] for k in dfkeys]


		# If no reference population is used, the data frame only needs two columns
		else:
			drb_df=pd.DataFrame(index=dfkeys,columns=['Background Distribution',pas.name])
		
		# Add target distribution values to the dataframe	
		drb_df['Background Distribution']=df_target_values
		
		#create a dictionary of subjects which uses their unique identifier as a key and a tuple for their two alleles
		subjects={line.split(',')[0]:(line.split(',')[1],line.split(',')[2].strip()) for line in open(samples_available).readlines()[1:]}
		
		#Just error checking the sample size being less than the samples available
		if N>=len(subjects.keys()):
			print("\n\nPlease choose a sample size less than the total number available.\n\n")
			quit()
		
		#Choose N subjects from the available sample
		current_sample=sample(subjects.keys(),N)
		
		#Create a dictionary of alleles and their respective frequencies in the sample chosen
		alleles=[]
		for subj in current_sample:
			alleles+=subjects[subj]
		current_distribution=Counter(alleles)
		current_distribution={allele:current_distribution[allele] for allele in current_distribution}
		
		#This is scaled to frequencies
		for allele in current_distribution:
			current_distribution[allele]/=float(N*2)
		
		#Add zeros for alleles not found in the subjects but found in the target Population
		#These are needed for getting average values later
		for allele in target_distribution:
			if allele not in current_distribution:
				current_distribution.update({allele:0})
		
		#Calculate the JSD score for the original sample
		current_score=get_jensen_score_sample(current_distribution,target_distribution)
		
		#scores keeps track of the JSD score throughout the algorithm
		scores.append(current_score)
		
		#Iterations is set to 0. It will increment by 1 till it reaches number_of_iterations
		i=0
		
		# current temp will start at 1 and exponentiall decay by a factor of temperature_decrease_alpha
		current_temp=starting_temperature
		
		
		while(i<number_of_iterations):
			i+=1
			
			# Choose a number/numbers from 1 to N. 
			new_subject_position=sample(range(N),changes)
			
			# Create a pool of samples which are not in the current sample
			new_subject_pool=[subject for subject in subjects.keys() if subject not in current_sample]
			
			# Create a copy of the current sample. Since you may reject the change
			# it is important to retain the sample you had at the start of this 
			# iteration.
			new_sample=current_sample[:]
			
			# Replace the subjects at the positions above with subjects from the pool of unused subjects
			for new_subject_number in new_subject_position:
				new_sample[new_subject_number]=sample(new_subject_pool,1)[0]
			
			# This is the same as line 196 - 213 but for the new sample
			alleles=[]
			for subj in new_sample:
				alleles+=subjects[subj]
			sample_distribution=Counter(alleles)
			for allele in sample_distribution:
				sample_distribution[allele]/=float(N*2)
			sample_distribution={allele:sample_distribution[allele] for allele in sample_distribution}
			for allele in target_distribution:
				if allele not in sample_distribution:
					sample_distribution.update({allele:0})
			new_score=get_jensen_score_sample(sample_distribution,target_distribution)
			
			# Generate a random number from unif(0,1) distribution
			r_num=random()
			
			# This if statement will replace the current sample with the new sample if
			# either of two conditions are met. If the new score is an improvement over
			# the current sample from the start of the iteration or if the random number 
			# generated is less than the current temperature which is decreasing.
			if new_score<current_score or r_num<current_temp:
				current_sample=new_sample[:]
				current_score=new_score
		
			#scores keeps track of the JSD score throughout the algorithm
			scores.append(current_score)
			
			# THe temperature is decreased		
			current_temp*=(1-temperature_decrease_alpha)
			
			# Just keeping track of progress 
			if pas.verbose:
				if i%1000==0:
					print('new:\t%f\tcurrent:\t%f\trandom:\t%f\ttemp:\t%f\tleft:\t%d\r'%(new_score,current_score,r_num,current_temp,number_of_iterations-i))
		'''
		Write a csv file with the donors chosen in the final iteration
		with their id and their alleles		
		'''
		outfile=open(outfilename+'_donors_selected.csv','w')
		outfile.write('Donor,Allele_1,Allele_2\n')
		for current_donor in current_sample:
			outfile.write('%s,%s,%s\n'%(current_donor,subjects[current_donor][0],subjects[current_donor][1]))
		outfile.close()	
		
		# Also saves the frequencies of alleles in both the background and chosen sample ( and reference sample
		# if provided)
		drb_df[pas.name]=[sample_distribution[k] for k in dfkeys]
		drb_df.to_csv(outfilename+'_frequencies.csv')
		
		# Bar chart of frequencies. It takes out alleles which have a frequency of less than 0.01
		# This will eliminate hundreds of alleles that are not relevant to see. It should get down to 
		# 30 or so alleles for DRB1 alleles
		ax=drb_df[drb_df.sum(axis=1)>.01].sort_index(axis=0).plot(kind='bar',cmap='Accent',width=.85)
		ax.set_ylabel('Frequency in Population')
		plt.tight_layout()
		fig=ax.get_figure()
		fig.savefig(outfilename+'.png',height=1200,width=800,dpi=600)
		print("Finished")

if __name__=='__main__':
	simulated_annealing()
	
