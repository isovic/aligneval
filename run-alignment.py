#! /usr/bin/python

import matplotlib.pyplot as plt

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
# GOLDEN_PATH = SCRIPT_PATH + '/../golden-bundle';
import sys
sys.path.append(SCRIPT_PATH + '/src');
sys.path.append(SCRIPT_PATH + '/wrappers');

import subprocess;

import evalalignments4;
from basicdefines import *;



if __name__ == "__main__":
	aligner_wrappers = find_files(WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py');
	for wrapper in aligner_wrappers:
		wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0];
		command = 'import %s' % (wrapper_basename);
		exec(command);

	machine_names = [];
	simulated_datasets = [];

	# simulated_datasets.append('Illumina-1k-single_end');								machine_names.append('illumina');
	# simulated_datasets.append('PacBio-1k');												machine_names.append('pacbio');
	# simulated_datasets.append('OxfordNanopore-pbsim-observed_last-2d-1k');				machine_names.append('nanopore');
	# simulated_datasets.append('OxfordNanopore-pbsim-observed_graphmap-2d-1k');			machine_names.append('nanopore');
	simulated_datasets.append('OxfordNanopore-pbsim-observed_last-1d-1k');				machine_names.append('nanopore');
	# simulated_datasets.append('OxfordNanopore-pbsim-observed_graphmap-1d-1k');			machine_names.append('nanopore');
	# simulated_datasets.append('OxfordNanopore-pbsim-observed_marginalign-2d-1k');		machine_names.append('nanopore');

	# simulated_datasets.append('Illumina-0k-single_end');								machine_names.append('illumina');
	# simulated_datasets.append('OxfordNanopore-pbsim-observed_marginalign-2d-0k');		machine_names.append('nanopore');

	genomes = [];
	# genomes.append('neisseria_meningitidis');
	# genomes.append('escherichia_coli');
	genomes.append('saccharomyces_cerevisiae');
	# genomes.append('caenorhabditis_elegans');
	# genomes.append('hg19_v38-chr3');

	num_processed_datasets = 0;
	num_datasets = len(simulated_datasets) * len(genomes);

	machine_num = 0;
	for simulated_dataset in simulated_datasets:
		for reference_name in genomes:
	# for dataset in datasets:
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
			# for wrapper in aligner_wrappers:
			# 	wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0];
			# 	command = '%s.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);' % (wrapper_basename);
			# 	exec(command);

			# if (machine_names[machine_num] == 'nanopore'):
			# 	wrapper_graphmap.run(reads_fasta, reference_file, 'nanopore', output_path, machine_suffix);
			# else:
			# wrapper_graphmap.run(reads_fasta, reference_file, 'nanopore', output_path, machine_suffix);
			# wrapper_graphmap.run(reads_fasta, reference_file, 'nanopore', output_path, 'test1-' + machine_suffix);

			wrapper_graphmap.run(reads_fasta, reference_file, 'test', output_path, 'test-' + machine_suffix);
			# 	wrapper_graphmap.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_names[machine_num] + '-' + machine_suffix);
			# wrapper_lastal.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);
			# wrapper_blasr.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);
			# wrapper_bwamem.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);
			# wrapper_blast.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);
			###########################################
			###########################################

			sys.stderr.write('[%d/%d] Finished simulated_dataset = "%s", reference_name = "%s"...\n' % ((num_processed_datasets + 1), num_datasets, simulated_dataset, reference_name));
			num_processed_datasets += 1;
			sys.stderr.write('====================================================\n');
				
		machine_num += 1;

	sys.stderr.write('Done!\n');
	sys.stderr.write('\n');
