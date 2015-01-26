#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
import subprocess;
import multiprocessing;

import fastqparser;
import measurement;

ALIGNER_PATH = SCRIPT_PATH + '/../aligners/blasr';



def Run(reads_file, reference_file, machine_name, machine_suffix, output_path):
	parameters = '';
	# suffix = 'params_' + machine_suffix;
	suffix = machine_suffix;
	num_threads = multiprocessing.cpu_count();
	
	if (('illumina' in machine_name.lower()) or ('roche' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1 -minMatch 7' % str(num_threads);
	elif (('pacbio' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1' % str(num_threads);
	elif (('real_nanopore' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1' % str(num_threads);
		
	elif (('10_percent' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1' % str(num_threads);
	elif (('20_percent' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1' % str(num_threads);
	elif (('30_percent' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1 -minMatch 7' % str(num_threads);
	elif (('40_percent' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1 -minMatch 7' % str(num_threads);
	elif (('50_percent' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1 -minMatch 7' % str(num_threads);

	elif (('debug' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1 -minMatch 7' % str(num_threads);
	elif (('experimental' in machine_name.lower())):
		parameters = '-nproc %s -sam -bestn 1 -minMatch 7' % str(num_threads);
		
	else:			# default
		parameters = '-nproc %s -sam -bestn 1' % str(num_threads);	# Maybe use: -bestn 1 ?

#def Run(reads_file, reference_file, machine_name, output_path):
	#num_threads = multiprocessing.cpu_count();
	
	#parameters = '';
	#if (machine_name == 'Illumina-single_end' or machine_name == 'Roche454-single_end'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_ngs';
	#elif (machine_name == 'PacBio-cov20'):
		#parameters = '-nproc %s -sam' % str(num_threads);
		#suffix = 'params_pacbio';
	#elif (machine_name == 'OxfordNanopore-pbsim-10_percent'):
		#parameters = '-nproc %s -sam' % str(num_threads);
		#suffix = 'params_10perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-20_percent'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_20perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_30perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-40_percent'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_40perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-50_percent'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_50perc';

	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent_observed_ratio'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_30perc_obs';
	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent_observed_ratio_lastal'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_30perc_obsl';
	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent_observed_ratio_graphmap'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_30perc_obsg';
	#elif (machine_name == 'OxfordNanopore-pbsim-40_percent_observed_ratio_lastal'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_40perc_obsl';
	#elif (machine_name == 'OxfordNanopore-pbsim-40_percent_observed_ratio_graphmap'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_40perc_obsg';

	#elif (machine_name == 'real_nanopore'):
		#parameters = '-nproc %s -sam -minMatch 7' % str(num_threads);
		#suffix = 'params_real_nanopore';
	#else:
		#parameters = '-nproc %s -sam' % str(num_threads);	# Maybe use: -bestn 1 ?
		#suffix = 'params_default';
	
	mapper_name = 'BLASR-%s' % suffix;
	
  #Options for Refining Hits.
   #-sdpTupleSize K (11)
               #Use matches of length K to speed dynamic programming alignments.  This controls
               #accuracy of assigning gaps in pairwise alignments once a mapping has been found,
               #rather than mapping sensitivity itself.
   #-scoreMatrix "score matrix string" 
               #Specify an alternative score matrix for scoring fasta reads.  The matrix is 
               #in the format 
                  #ACGTN
                #A abcde
                #C fghij
                #G klmno
                #T pqrst
                #N uvwxy . The values a...y should be input as a quoted space separated 
               #string: "a b c ... y". Lower scores are better, so matches should be less 
               #than mismatches e.g. a,g,m,s = -5 (match), mismatch = 6. 
   #-affineOpen value (10) 
               #Set the penalty for opening an affine alignment.
   #-affineExtend a (0)
               #Change affine (extension) gap penalty. Lower value allows more gaps.
 #Alignment Output.
   #-bestn n (10)
               #Report the top 'n' alignments.
	
	reads_basename = os.path.splitext(os.path.basename(reads_file))[0];
	sam_file = '%s/%s.sam' % (output_path, mapper_name);
	memtime_file = '%s/%s.memtime' % (output_path, mapper_name);
	memtime_file_index = '%s/%s-index.memtime' % (output_path, mapper_name);
	
	print '[Testing BLASR] Creating index...';
	command = '%s %s/alignment/bin/sawriter %s.blasrsa %s' % (measurement.MeasureCommand(memtime_file_index), ALIGNER_PATH, reference_file, reference_file);
	print '[Testing BLASR] ', command;
	subprocess.call(command, shell=True);
	print ' ';
	
	print '[Testing BLASR] Running BLASR...';
	command = '%s %s/alignment/bin/blasr %s %s %s -sa %s.blasrsa -out %s' % (measurement.MeasureCommand(memtime_file), ALIGNER_PATH, reads_file, reference_file, parameters, reference_file, sam_file);
	print '[Testing BLASR] ', command;
	subprocess.call(command, shell=True);
	print ' ';
	
	print '[Testing BLASR] BLASR wrapper script finished processing.';
	


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
