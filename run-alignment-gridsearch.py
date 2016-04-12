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

	num_processed_datasets = 0;
	num_datasets = len(simulated_datasets) * len(genomes);

	### Grid search to estimate scalability with error rates and read lengths.
	machine_num = 0;
	for simulated_dataset in simulated_datasets_grid:
		for reference_name in genomes_grid:
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
			machine_suffix = 'v2';
			wrapper_graphmap.run(reads_fasta, reference_file, 'default', output_path, 'default-' + machine_suffix);
			wrapper_graphmap.run(reads_fasta, reference_file, 'gotoh', output_path, 'gotoh-' + machine_suffix);
			wrapper_graphmap.run(reads_fasta, reference_file, 'anchor', output_path, 'anchor-' + machine_suffix);
			sys.stderr.write('[%d/%d] (Grid) Finished simulated_dataset = "%s", reference_name = "%s"...\n' % ((num_processed_datasets + 1), num_datasets, simulated_dataset, reference_name));
			num_processed_datasets += 1;
			sys.stderr.write('====================================================\n');
		machine_num += 1;

	sys.stderr.write('Done!\n');
	sys.stderr.write('\n');
