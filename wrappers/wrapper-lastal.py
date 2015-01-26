#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
import subprocess;
import multiprocessing;

import fastqparser;
import measurement;

#LAST_PATH = '/home/ivan/work/eclipse-workspace/aligner-comparison/last-475/src';
ALIGNER_PATH = SCRIPT_PATH + '/../aligners/last-475/src';



def GetSAMHeader(reference_file):
	[headers, seqs, quals] = fastqparser.read_fastq(reference_file);
	
	line = '';
	
	i = 0;
	while i < len(headers):
		line += '@SQ\tSN:%s\tLN:%d\n' % (headers[i], len(seqs[i]));
		i += 1;
	
	return line;

def Run(reads_file, reference_file, machine_name, machine_suffix, output_path):
	parameters = '';
	# suffix = 'params_' + machine_suffix;
	suffix = machine_suffix;
	
	if ((machine_name.lower() == 'illumina' or machine_name.lower() == 'roche')):
		parameters = '-v ';
	elif ((machine_name.lower() == 'pacbio')):
		parameters = '-v -q 1 -r 1 -a 1 -b 1';
	elif ((machine_name.lower() == 'real_nanopore')):
		parameters = '-v -q 1 -r 1 -a 1 -b 1';

	elif ((machine_name.lower() == 'q1t1')):
		parameters = '-v -q 1 -r 1 -a 1 -b 1 -T 1';

	elif ((machine_name.lower() == 'simulated')):
		parameters = '-v -q 1 -r 1 -a 1 -b 1 -T 1';

	# elif (('pacbio' in machine_name.lower())):
	# 	parameters = '-v -q 1 -r 1 -a 1 -b 1';
	# elif (('10_percent' in machine_name.lower())):
	# 	parameters = '-v -q 1 -r 1 -a 1 -b 1';
	# elif (('20_percent' in machine_name.lower())):
	# 	parameters = '-v -q 1 -r 1 -a 1 -b 1';
	# elif (('30_percent' in machine_name.lower())):
	# 	parameters = '-v -q 1 -r 1 -a 1 -b 1';
	# elif (('40_percent' in machine_name.lower())):
	# 	parameters = '-v -q 1 -r 1 -a 1 -b 1';
	# elif (('50_percent' in machine_name.lower())):
	# 	parameters = '-v -q 1 -r 1 -a 1 -b 1';
	# elif (('real_nanopore' in machine_name.lower())):
	# 	parameters = '-v -q 1 -r 1 -a 1 -b 1';
	elif ((machine_name.lower() == 'debug')):
		parameters = '-v -q 1 -r 1 -a 1 -b 1';
		
	elif ((machine_name.lower() == 'q1')):
		parameters = '-v -q 1 -r 1 -a 1 -b 1';
	elif ((machine_name.lower() == 'q2')):
		parameters = '-v -q 2 -r 1 -a 1 -b 1';
	
	else:			# default
		parameters = '-v ';	# The default.

#def Run(reads_file, reference_file, machine_name, output_path):
	#parameters = '';
	#if (machine_name == 'Illumina-single_end' or machine_name == 'Roche454-single_end'):
		#parameters = '-v';
		#suffix = 'params_ngs';
	#elif (machine_name == 'PacBio-cov20'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_pacbio';
	#elif (machine_name == 'OxfordNanopore-pbsim-10_percent'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_10perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-20_percent'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_20perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_30perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-40_percent'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_40perc';
	#elif (machine_name == 'OxfordNanopore-pbsim-50_percent'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_50perc';

	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent_observed_ratio'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_30perc_obs';
	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent_observed_ratio_lastal'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_30perc_obsl';
	#elif (machine_name == 'OxfordNanopore-pbsim-30_percent_observed_ratio_graphmap'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_30perc_obsg';
	#elif (machine_name == 'OxfordNanopore-pbsim-40_percent_observed_ratio_lastal'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_40perc_obsl';
	#elif (machine_name == 'OxfordNanopore-pbsim-40_percent_observed_ratio_graphmap'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_40perc_obsg';

	#elif (machine_name == 'real_nanopore'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_real_nanopore';
	#elif (machine_name == 'q1'):
		#parameters = '-v -q 1 -r 1 -a 1 -b 1';
		#suffix = 'params_q1';
	#elif (machine_name == 'q2'):
		#parameters = '-v -q 2 -r 1 -a 1 -b 1';
		#suffix = 'params_q2';
	#else:
		#parameters = '-v';	# The default.
		#suffix = 'params_default';
	
	mapper_name = 'LAST-%s' % suffix;

	reads_fasta = reads_file;
	reads_basename = os.path.splitext(os.path.basename(reads_file))[0];
	maf_file = '%s/%s.maf' % (output_path, mapper_name);
	sam_file = '%s/%s.sam' % (output_path, mapper_name);
	memtime_file = '%s/%s.memtime' % (output_path, mapper_name);
	memtime_file_index = '%s/%s-index.memtime' % (output_path, mapper_name);
	memtime_file_maftosam = '%s/%s-maftosam.memtime' % (output_path, mapper_name);



	if (reads_file[-1] == 'q'):
		print '[Testing LAST-475] Converting FASTQ to FASTA...';
		reads_fasta = reads_file[0:-1] + 'a';
		#command = '%s/fastq2fasta.py %s %s' % (TOOLS_PATH, reads_file, reads_fasta);
		#print command;
		#subprocess.call(command, shell=True);
		fastqparser.convert_to_fasta(reads_file, reads_fasta);
		print ' ';
	
	reference_db_file = reference_file + '.db';
	print '[Testing LAST-475] Running Lastdb...';
	command = '%s %s/lastdb %s %s' % (measurement.MeasureCommand(memtime_file_index), ALIGNER_PATH, reference_db_file, reference_file);
	print '[Testing LAST-475] ', command;
	subprocess.call(command, shell=True);
	print ' ';
	
	print '[Testing LAST-475] Running Lastal...';
	if (machine_name.lower() == 'simulated'):
		# command = '%s %s/lastal %s %s %s | %s/last-split > %s' % (measurement.MeasureCommand(memtime_file), ALIGNER_PATH, parameters, reference_db_file, reads_fasta, ALIGNER_PATH, maf_file);
		# print '[Testing LAST-475] Running with last-split!';
		command = '%s %s/lastal %s %s %s > %s' % (measurement.MeasureCommand(memtime_file), ALIGNER_PATH, parameters, reference_db_file, reads_fasta, maf_file);
	else:
		command = '%s %s/lastal %s %s %s > %s' % (measurement.MeasureCommand(memtime_file), ALIGNER_PATH, parameters, reference_db_file, reads_fasta, maf_file);
		print '[Testing LAST-475] Not using last-split!';
	print '[Testing LAST-475] ', command;
	subprocess.call(command, shell=True);
	print ' ';

	print '[Testing LAST-475] Converting output MAF to SAM file format...';
#	command = 'cat "%s" > %s' % (GetSAMHeader(reference_file), sam_file);
#	print command;
#	subprocess.call(command, shell=True);
	fp = open(sam_file, 'w');
	fp.write(GetSAMHeader(reference_file));
	fp.close();
	command = '%s %s/../scripts/maf-convert.py sam %s >> %s' % (measurement.MeasureCommand(memtime_file_maftosam), ALIGNER_PATH, maf_file, sam_file);
	print '[Testing LAST-475] ', command;
	subprocess.call(command, shell=True);
	print ' ';
	
	print '[Testing LAST-475] Lastal wrapper script finished processing.';


if __name__ == "__main__":
	if (len(sys.argv) != 5):
		print 'Usage:';
		print '\t%s <reads_file> <reference_file> <machine_name> <output_path>' % sys.argv[0];
		print ' ';
		print '\t<machine_name>\tJust a symbolic name for preset parameters (preset within this script.';
		print '\t\t\tUse "q1" or "q2" for the options commonly used with Oxford Nanopore data.';
		exit(1);
	
	reads_file = sys.argv[1];
	reference_file = sys.argv[2];
	machine_name = sys.argv[3];
	output_path = sys.argv[4];
	
	Run(reads_file, reference_file, machine_name, machine_name, output_path);
