#! /usr/bin/python

import matplotlib.pyplot as plt
import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
import sys
sys.path.append(SCRIPT_PATH + '/src');
sys.path.append(SCRIPT_PATH + '/wrappers');
import subprocess;

import evalalignments4;
from basicdefines import *;
from dataset_specification import *

if __name__ == "__main__":
	setup_measure_command();

	aligner_wrappers = find_files(WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py');
	for wrapper in aligner_wrappers:
		wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0];
		command = 'import %s' % (wrapper_basename);
		exec(command);
	
	use_blast = False;
        blast_yes_no = raw_input("BLAST may take a long time to execute on some datasets while consuming a lot of disk space as well. Do you want to run it? [y/n] ");
        if (blast_yes_no == 'y'):
		use_blast = True;

	use_marginalign = False;
        marginalign_yes_no = raw_input("marginAlign may take a long time to execute on some datasets while consuming a lot of disk space as well. Do you want to run it? [y/n] ");
        if (marginalign_yes_no == 'y'):
		use_marginalign = True;

	sys.stderr.write('\n');
	sys.stderr.write('Using BLAST: %s\n' % ('yes' if (use_blast == True) else 'no'));
	sys.stderr.write('Using marginAlign: %s\n' % ('yes' if (use_marginalign == True) else 'no'));
	sys.stderr.write('\n');

	num_processed_datasets = 0;
	num_datasets = len(simulated_datasets) * len(genomes);

	machine_num = 0;
	for simulated_dataset in simulated_datasets:
		for reference_name in genomes:
			sys.stderr.write('[%d/%d] Starting simulated_dataset = "%s", reference_name = "%s"...\n' % ((num_processed_datasets + 1), num_datasets, simulated_dataset, reference_name));
			reference_file = '%s/%s.fa' % (REFERENCE_GENOMES_ROOT_ABS, reference_name);
			reads_fastq = '%s/%s/%s/reads.fq' % (READS_SIMULATED_ROOT_ABS, simulated_dataset, reference_name);
			reads_fasta = '%s/%s/%s/reads.fa' % (READS_SIMULATED_ROOT_ABS, simulated_dataset, reference_name);
			output_path = '%s/reads-simulated/%s/%s' % (EVALUATION_PATH_ROOT_ABS, simulated_dataset, reference_name);
			
			if not os.path.exists(output_path):
				sys.stderr.write('Creating folder "%s".\n' % (output_path));
				os.makedirs(output_path);

			###########################################
			###########################################
			machine_suffix = 'v1';
			# machine_suffix = 'release-v1';
			# for wrapper in aligner_wrappers:
			# 	wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0];
			# 	command = '%s.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);' % (wrapper_basename);
			# 	exec(command);

			wrapper_graphmap.run(reads_fasta, reference_file, 'default', output_path, 'default-' + machine_suffix);
			wrapper_graphmap.run(reads_fasta, reference_file, 'gotoh', output_path, 'gotoh-' + machine_suffix);
			wrapper_graphmap.run(reads_fasta, reference_file, 'anchor', output_path, 'anchor-' + machine_suffix);

			wrapper_lastal.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);
			wrapper_blasr.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);
			wrapper_bwamem.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);
			wrapper_daligner.run('align', reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);

			if (use_blast == True):
				wrapper_blast.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);

			if (use_marginalign == True):
				machine_name = machine_names[machine_num];
				if (machine_names[machine_num] == 'nanopore' and ('-1d' in simulated_dataset)):
					machine_name += '1d';
				elif (machine_names[machine_num] == 'nanopore' and ('-2d' in simulated_dataset)):
					machine_name += '2d';
				wrapper_marginalign.run('run', reads_fastq, reference_file, machine_name, output_path, machine_suffix);
				wrapper_marginaligngraphmap.run('run', reads_fastq, reference_file, machine_name, output_path, machine_suffix);
				wrapper_marginaligngraphmap.run('run', reads_fastq, reference_file, 'anchor', output_path, 'anchor-' + machine_suffix);
			
			###########################################
			###########################################

			sys.stderr.write('[%d/%d] Finished simulated_dataset = "%s", reference_name = "%s"...\n' % ((num_processed_datasets + 1), num_datasets, simulated_dataset, reference_name));
			num_processed_datasets += 1;
			sys.stderr.write('====================================================\n');
				
		machine_num += 1;

	sys.stderr.write('Done!\n');
	sys.stderr.write('\n');
