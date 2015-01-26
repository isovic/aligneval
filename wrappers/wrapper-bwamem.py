#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
import subprocess;
import multiprocessing;

import fastqparser;
import measurement;

ALIGNER_PATH = SCRIPT_PATH + '/../aligners/bwa';



def Run(reads_file, reference_file, machine_name, machine_suffix, output_path):
	parameters = '';
	# suffix = 'params_' + machine_suffix;
	suffix = machine_suffix;
	num_threads = multiprocessing.cpu_count();
	
	if (('illumina' in machine_name.lower()) or ('roche' in machine_name.lower())):
		parameters = '-t %s' % str(num_threads);
	elif (('pacbio' in machine_name.lower())):
		parameters = '-t %s -x pacbio' % str(num_threads);
	elif (('real_nanopore' in machine_name.lower())):
		parameters = '-t %s -x ont2d' % str(num_threads);
		
	elif (('10_percent' in machine_name.lower())):
		parameters = '-t %s -x pacbio' % str(num_threads);
	elif (('20_percent' in machine_name.lower())):
		parameters = '-t %s -x pacbio' % str(num_threads);
	elif (('30_percent' in machine_name.lower())):
		parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Mapped: 47.30%
	elif (('40_percent' in machine_name.lower())):
		parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Mapped: 47.30%
	elif (('50_percent' in machine_name.lower())):
		parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Mapped: 47.30%
		
	elif (('real_nanopore' in machine_name.lower())):
		parameters = '-t %s -x ont2d -v 3' % str(num_threads); # Mapped: 47.30%

	# 	parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Mapped: 47.30%
	# elif (('debug' in machine_name.lower())):
	# 	parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Mapped: 47.30%
	# elif (('experimental' in machine_name.lower())):
	# 	#parameters = '-t %s -c1000 -k11 -W40 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # oko 14.5% mapirano
	# 	#parameters = '-t %s -c1000 -k10 -W40 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # oko 17.5% mapirano
	# 	#parameters = '-t %s -c1000 -k9 -W40 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # oko 15.5% mapirano
	# 	#parameters = '-t %s -c1000 -k10 -W0 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Skrsi se
	# 	#parameters = '-t %s -c1000 -k10 -W20 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Skrsi se
	# 	#parameters = '-t %s -c1000 -k10 -W30 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Mapped: 47.30%
	# 	parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Mapped: 47.30%
	# 	#parameters = '-t %s -c1000 -k10 -W30 -r 5 -A1 -B1 -O1 -E1 -L0 -v 3' % str(num_threads); # Mapped: 43.66%
	# 	#parameters = '-t %s -c1000 -k10 -W30 -r 5 -A1 -B2 -O1 -E1 -L0 -v 3' % str(num_threads); # Mapped: 29.51%
	# 	#parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B2 -O3,5 -E2,1 -L0 -v 3' % str(num_threads); # Mapped: 33.56%
	
	else:			# default
		parameters = '-t %s' % str(num_threads);	# The default.			-t %s -h 1

#def Run(reads_file, reference_file, machine_name, output_path):
	#num_threads = multiprocessing.cpu_count();
	
	#parameters = '';
	#if (machine_name == 'Illumina-single_end' or machine_name == 'Roche454-single_end'):
		#parameters = '-t %s' % str(num_threads);
		#suffix = 'params_ngs';
	#elif (machine_name == 'PacBio-cov20'):
		#parameters = '-t %s -x pacbio' % str(num_threads);
		#suffix = 'params_pacbio';
	#elif (machine_name == 'OxfordNanopore-pbsim-10_percent'):
		#parameters = '-t %s -x pacbio' % str(num_threads);
		#suffix = 'params_10perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-20_percent'):
		#parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#suffix = 'params_20perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent'):
		#parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#suffix = 'params_30perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-40_percent'):
		##parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads);
		#suffix = 'params_40perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-50_percent'):
		#parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#suffix = 'params_50perc';

	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent_observed_ratio'):
		#parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#suffix = 'params_30perc_obs';
	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent_observed_ratio_lastal'):
		#parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#suffix = 'params_30perc_obsl';
	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent_observed_ratio_graphmap'):
		#parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#suffix = 'params_30perc_obsg';
	#elif (machine_name == 'OxfordNanopore-pbsim-40_percent_observed_ratio_lastal'):
		#parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#suffix = 'params_40perc_obsl';
	#elif (machine_name == 'OxfordNanopore-pbsim-40_percent_observed_ratio_graphmap'):
		#parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#suffix = 'params_40perc_obsg';

	#elif (machine_name == 'real_nanopore'):
		#parameters = '-t %s -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0' % str(num_threads);
		#suffix = 'params_real_nanopore';
	
	##-k INT        minimum seed length [19]
	##-W INT        discard a chain if seeded bases shorter than INT [0]
	##-c INT        skip seeds with more than INT occurrences [500]
	##-r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
	## 
	##-A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
	##-B INT        penalty for a mismatch [4]
	##-O INT[,INT]  gap open penalties for deletions and insertions [6,6]
	##-E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
	##-L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
	## 
	##-x STR        read type. Setting -x changes multiple parameters unless overriden [null]
                     ##pacbio: -k17 -W40 -r10 -A2 -B5 -O2 -E1 -L0
                     ##pbread: -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -N25 -FeaD.001
	## -h INT        if there are <INT hits with score >80% of the max score, output all in XA [5]

	#elif (machine_name == 'real_experimental'):
		##parameters = '-t %s -c1000 -k11 -W40 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # oko 14.5% mapirano
		##parameters = '-t %s -c1000 -k10 -W40 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # oko 17.5% mapirano
		##parameters = '-t %s -c1000 -k9 -W40 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # oko 15.5% mapirano
		##parameters = '-t %s -c1000 -k10 -W0 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Skrsi se
		##parameters = '-t %s -c1000 -k10 -W20 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Skrsi se
		##parameters = '-t %s -c1000 -k10 -W30 -r10 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Mapped: 47.30%
		#parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B5 -O2 -E1 -L0 -v 3' % str(num_threads); # Mapped: 47.30%
		##parameters = '-t %s -c1000 -k10 -W30 -r 5 -A1 -B1 -O1 -E1 -L0 -v 3' % str(num_threads); # Mapped: 43.66%
		##parameters = '-t %s -c1000 -k10 -W30 -r 5 -A1 -B2 -O1 -E1 -L0 -v 3' % str(num_threads); # Mapped: 29.51%
		##parameters = '-t %s -c1000 -k10 -W30 -r 5 -A2 -B2 -O3,5 -E2,1 -L0 -v 3' % str(num_threads); # Mapped: 33.56%
		
		#suffix = 'params_real_experimental';
		
	#else:
		#parameters = '-t %s' % str(num_threads);	# The default.			-t %s -h 1
		#suffix = 'params_default';
	
	mapper_name = 'BWAMEM-%s' % suffix;
	
	reads_basename = os.path.splitext(os.path.basename(reads_file))[0];
	sam_file = '%s/%s.sam' % (output_path, mapper_name);
	memtime_file = '%s/%s.memtime' % (output_path, mapper_name);
	memtime_file_index = '%s/%s-index.memtime' % (output_path, mapper_name);
	
	print '[Testing BWA-MEM] Creating index...';
	command = '%s %s/bwa index %s' % (measurement.MeasureCommand(memtime_file_index), ALIGNER_PATH, reference_file);
	print '[Testing BWA-MEM] ', command;
	subprocess.call(command, shell=True);
	print ' ';
	
	print '[Testing BWA-MEM] Running BWA-MEM...';
	command = '%s %s/bwa mem %s %s %s > %s' % (measurement.MeasureCommand(memtime_file), ALIGNER_PATH, parameters, reference_file, reads_file, sam_file);
	print '[Testing BWA-MEM] ', command;
	subprocess.call(command, shell=True);
	print ' ';
	
	print '[Testing BWA-MEM] BWA-MEM wrapper script finished processing.';
	


if __name__ == "__main__":
	if (len(sys.argv) != 5):
		print 'Usage:';
		print '\t%s <reads_file> <reference_file> <machine_name> <output_path>' % sys.argv[0];
		exit(1);
	
	reads_file = sys.argv[1];
	reference_file = sys.argv[2];
	machine_name = sys.argv[3];
	output_path = sys.argv[4];
	
	Run(reads_file, reference_file, machine_name, machine_name, output_path);
