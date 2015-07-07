#! /usr/bin/python

import re;

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
sys.path.append(SCRIPT_PATH + '/../src');

import subprocess;
import multiprocessing;

import basicdefines;

import fastqparser;

ALIGNER_URL = 'https://github.com/thegenemyers/DALIGNER.git';
ALIGNER_DB_URL = 'https://github.com/thegenemyers/DAZZ_DB.git';
ALIGNER_PATH = SCRIPT_PATH + '/../aligners/DALIGNER/DALIGNER/';
ALIGNER_DB_PATH = SCRIPT_PATH + '/../aligners/DALIGNER/DAZZ_DB/';
BIN = 'daligner';
MAPPER_NAME = 'DALIGNER';



def convert_reads_to_pacbio_format(reads_file, daligner_reads_file):
	try:
		fp_in = open(reads_file, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % reads_file);
		exit(0);

	try:
		fp_out = open(daligner_reads_file, 'w');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % daligner_reads_file);
		exit(0);

	current_read = 0;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		current_read += 1;

		### Check if the read is already formatted like PacBio.
		if (header.count('/') == 2 and 'RQ) in header'):
			fp_out.write('\n'.join(read) + '\n');
			continue;

		trimmed_header = header.replace('_', ' ').split()[0];
		pacbio_header = '%s/%d/0_%d RQ=0.850' % (trimmed_header, current_read, len(read[1]));
		read[0] = '%s%s' % (read[0][0], pacbio_header); ### Keep the first char of the header line.
		read[1] = re.sub("(.{500})", "\\1\n", read[1], 0, re.DOTALL);	### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
		if (len(read) == 4):
			read[3] = re.sub("(.{500})", "\\1\n", read[3], 0, re.DOTALL);	### Wrap the qual line, because DALIGNER has a 9998bp line len limit.
		fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();

# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#	reads_file			Path to a FASTA/FASTQ file containing reads.
#	reference_file		Path to a reference genome FASTA file.
#	machine_name		A symbolic name to specify a set of parameters for a specific sequencing platform.
#	output_path			Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#	output_suffix		A custom suffix that can be added to the output filename.
def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
	parameters = '';
	num_threads = multiprocessing.cpu_count() / 2;

	if ((machine_name.lower() == 'illumina') or (machine_name.lower() == 'roche')):
		parameters = '-t %s' % str(num_threads);

	elif ((machine_name.lower() == 'pacbio')):
		# parameters = '-t %s -x pacbio' % str(num_threads);
		parameters = '-vd';

	elif ((machine_name.lower() == 'nanopore')):
		parameters = '-t %s -x ont2d' % str(num_threads);

	elif ((machine_name.lower() == 'debug')):
		parameters = '-t %s' % str(num_threads);

	else:			# default
		parameters = '-t %s' % str(num_threads);



	if (output_suffix != ''):
		output_filename = '%s-%s' % (MAPPER_NAME, output_suffix);
	else:
		output_filename = MAPPER_NAME;
	
	reads_basename = os.path.splitext(os.path.basename(reads_file))[0];
	sam_file = '%s/%s.sam' % (output_path, output_filename);
	memtime_file = '%s/%s.memtime' % (output_path, output_filename);
	memtime_file_index = '%s/%s-index.memtime' % (output_path, output_filename);

	print reads_file;

	index_file = reference_file + '.dam';
	
	# Run the indexing process, and measure execution time and memory.
	if (True or (not os.path.exists(index_file))):
		sys.stderr.write('[%s wrapper] Copying reference to satisfy the extension requirements...\n' % (MAPPER_NAME));
		command = 'cp %s %s.fasta' % (reference_file, reference_file);
		subprocess.call(command, shell=True);
		sys.stderr.write('\n');

		sys.stderr.write('[%s wrapper] Generating index...\n' % (MAPPER_NAME));
		command = '%s %s/fasta2DAM %s %s' % (basicdefines.measure_command(memtime_file_index), ALIGNER_DB_PATH, index_file, reference_file);
		sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
		subprocess.call(command, shell=True);
		sys.stderr.write('\n\n');
	else:
		sys.stderr.write('[%s wrapper] Reference index already exists. Continuing.\n' % (MAPPER_NAME));
		sys.stderr.flush();

	# return '';

	daligner_reads_file = '%s-daligner.fasta' % (os.path.splitext(reads_file)[0]);
	if (True or (not os.path.exists(daligner_reads_file))):
		sys.stderr.write('[%s wrapper] Modifying the reads file to have PacBio headers...\n' % (MAPPER_NAME));
		# command = 'cp %s %s.fasta' % (reads_file, reads_file);
		# subprocess.call(command, shell=True);
		convert_reads_to_pacbio_format(reads_file, daligner_reads_file);
		sys.stderr.write('\n');

	sys.stderr.write('[%s wrapper] Converting the reads file into a DB file...\n' % (MAPPER_NAME));
	command = '%s %s/fasta2DB %s.db %s' % (basicdefines.measure_command(memtime_file_index), ALIGNER_DB_PATH, daligner_reads_file, daligner_reads_file);
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n\n');

	# Run the alignment process, and measure execution time and memory.
	sys.stderr.write('[%s wrapper] Running %s...\n' % (MAPPER_NAME, MAPPER_NAME));
	command = '%s %s/HPCmapper %s %s %s > run-daligner.sh' % (basicdefines.measure_command(memtime_file), ALIGNER_PATH, parameters, reference_file, daligner_reads_file);
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n\n');
	
	sys.stderr.write('[%s wrapper] %s wrapper script finished processing.\n' % (MAPPER_NAME, MAPPER_NAME));

	return sam_file


# This is a standard interface for setting up the aligner. It should assume that the aligner
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
	sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (MAPPER_NAME, MAPPER_NAME));
	sys.stderr.write('[%s wrapper] Creating a folder for all %s repos...\n' % (MAPPER_NAME, MAPPER_NAME));
	command = 'mkdir -p %s/%s' % (basicdefines.ALIGNERS_PATH_ROOT_ABS, MAPPER_NAME);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('[%s wrapper] Cloning git repository.\n' % (MAPPER_NAME));
	command = 'cd %s/%s; git clone %s' % (basicdefines.ALIGNERS_PATH_ROOT_ABS, MAPPER_NAME, ALIGNER_URL);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('[%s wrapper] Running make.\n' % (MAPPER_NAME));
	command = 'cd %s; make' % (ALIGNER_PATH);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('[%s wrapper] Cloning git repository.\n' % (MAPPER_NAME));
	command = 'cd %s/%s; git clone %s' % (basicdefines.ALIGNERS_PATH_ROOT_ABS, MAPPER_NAME, ALIGNER_DB_URL);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');	

	sys.stderr.write('[%s wrapper] Running make.\n' % (MAPPER_NAME));
	command = 'cd %s; make' % (ALIGNER_DB_PATH);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	# sys.stderr.write('[%s wrapper] Checking out commit "eb428d7d31ced059ad39af2701a22ebe6d175657" for reproducibility purposes.\n' % (MAPPER_NAME));
	# command = 'cd %s; git checkout eb428d7d31ced059ad39af2701a22ebe6d175657' % (ALIGNER_PATH);
	# subprocess.call(command, shell='True');
	# sys.stderr.write('\n');

	sys.stderr.write('[%s wrapper] All installation steps finished.\n' % (MAPPER_NAME));
	sys.stderr.write('\n');



def verbose_usage_and_exit():
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s mode [<reads_file> <reference_file> <machine_name> <output_path> [<output_suffix>]]\n' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\t- mode          - either "run" or "install". If "install" other parameters can be ommitted.\n');
	sys.stderr.write('\t- machine_name  - "illumina", "roche", "pacbio", "nanopore" or "default".\n');
	sys.stderr.write('\t- output_suffix - suffix for the output filename.\n');

	exit(0);

if __name__ == "__main__":
	if (len(sys.argv) < 2 or len(sys.argv) > 7):
		verbose_usage_and_exit();

	if (sys.argv[1] == 'install'):
		download_and_install();
		exit(0);

	elif (sys.argv[1] == 'run'):
		if (len(sys.argv) < 6):
			verbose_usage_and_exit();

		reads_file = sys.argv[2];
		reference_file = sys.argv[3];
		machine_name = sys.argv[4];
		output_path = sys.argv[5];
		output_suffix = '';

		if (len(sys.argv) == 7):
			output_suffix = sys.argv[6];
		run(reads_file, reference_file, machine_name, output_path, output_suffix);

	else:
		verbose_usage_and_exit();
