# SampPick
** SampPick Sample Picking Algorithm ***

This software runs a simulated annealing algorithm on a user supplied collection of samples. The software will iteratively subsitute sample individuals to create an sample cohort of size N (useer specified). The cohorts will be analyzed by comparing the frequency distribution of the alleles as compared to a background distribution of interest. We have supplied a background frequency distribution for alleles found in the North American population. In addition, a reference sample could be added with the command -r so that you can add a premade cohort to compare your optimized cohort against. This reference sample should be in the same format as the available sample. (See section 2)


The sofware will ouput three files. The first is a bar chart of alleles and their frequencies in the background distribution, the optimized sample( and optionally the reference) The second is a csv file of the two or three frequency distributions described above. The third is the list of samples chosen for the optimized samples. This list will present the sample and their respective DRB alleles

Additionally, there is a second mode accessible through the flag --basic. It will allow for the simple output of a bar chart of a target frequency and a reference population as well as the associated frequencies.

Please cite : 

McGill JR, Yogurtcu ON, Verthelyi D, Yang H, Sauna ZE (2019) SampPick: Selection of a Cohort of Subjects Matching a Population HLA Distribution. Front. Immunol. 10:2894. doi: 10.3389/fimmu.2019.02894

Instructions for Use.

1) Requirements
	
	a) Python (tested on 2.7, and 3.7) 
	 
	b) numPy - http://www.numpy.org (tested on v1.16.4)
	
	c) pandas - http://pandas.pydata.org (tested on v0.25.1)
	
	d) matplotlib - http://www.matplotlib.org (tested on v3.1.1)
	
	e) This program has only been tested on a UNIX or UNIX-like machine

2) Getting ready to Run/Files needed:

	a) SampPick.py - This is the entire program

	b) a background distribution - this is selected with the -t option. if it is blank, the software will use the supplied background distribution for North America (NA_DRB1_frequencies.csv). Additionally this file can be used as a reference for the format of your own file.

	c) available samples and reference sample. This file is required. It is expected to be a csv file with the first line being a header. Each subsequent line should be 

	Unique Subject Identifier , DRB1 allele 1 , DRB1 allele 2

	Included is the list of available samples used in the paper above as a reference.


3) Running SampPick
	
	The program is run with the command :

	python SampPick.py [options] available_sample_file

	Required:
	available           Filename of csv  see above for help

	Options:
	  
	  -h, --help          show help message and exit
	  
	  -s , --size         Sample size desired (Integer) (default: 50)
	  
	  -r , --reference    CSV file to compare optmized sample to 
	                      
	  -i , --iterations   Number of iterations for Simulated Annealing (Integer) (default: 10000)

	  -a , --alpha        Temperature decrease Alpha value (Float) (default: 0.0007)

	  -c , -changes       Number of changes made on each iteration of the algorithm (Integer) (default: 1)
	  
	  -t , --target       Filename of csv file (see README for help) (default: NA_DRB1_frequencies.csv)
	 
	  -o , --outfile      Name for the csv and Png files created (default: Results)

	  -n , --name         Name for the Reference Sample (default: Optimized Sample)

	  --basic             just analyze the available sample in its entirety (default: False)
	  
	  --verbose           print output every 1000 iterations (default: False)

	  The user options are entirely optional. It should be noted that the number of iterations and temperature decrease alpha values are linked. It is important that (1-alpha)^x approaches 0 at the correct speed for the algorithm to run correctly. A quick way to inspect your values is to plot them in R # plot(1:(iteration),(1-alpha)**(1:iterations)) and approximate the shape of the curve produced with the default values.
	
4) Modifications for non standard use.
   
	A) Using a biallelic trait other than HLA-DRB1 Make sure that the first line of both the background frequencies and the available donors have header lines. The alleles can be changed to anything the user wants as long as the alleles for each user have entries in the background frequencies.

	B) For a mono-allelic trait, simply repeating the first allele column in the available patients. 

	C) If SampPick is used to optimize a population for multiple alleles at once, the user has the option of calculating the joint frequency of each combination of alleles. It is up to the user to decide whether the frequency of a combination of alleles is simply the product of the two frequencies (Assuming independence of the alleles) or if there is some linkage disequilibrium that needs to be accounted for. After the background distribution is calculated, the user can then proceed as in B above, using each combination of alleles as a mono-allelic trait. 

4) Questions regarding this algorithm to select a cohort matching a target distribution should be sent to joseph.mcgill@fda.hhs.gov or zuben.sauna@fda.hhs.gov

