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

import wrapper_graphmap;
import wrapper_lastal;
import wrapper_bwamem;
import wrapper_blasr;

if __name__ == "__main__":
	sam_suffixes = [];
	machines = [];

	# machines.append('OxfordNanopore-pbsim-observed_last-2d-10k');	sam_suffixes.append('real_nanopore');
	# machines.append('OxfordNanopore-pbsim-observed_graphmap-2d-10k');	sam_suffixes.append('real_nanopore');
	# machines.append('OxfordNanopore-pbsim-observed_last-1d-10k');	sam_suffixes.append('real_nanopore');
	# machines.append('OxfordNanopore-pbsim-observed_graphmap-1d-10k');	sam_suffixes.append('real_nanopore');

	machines.append('Illumina-10k-single_end');	sam_suffixes.append('illumina');
	# machines.append('PacBio-10k');	sam_suffixes.append('pacbio');
	
	genomes = [];
	# genomes.append('escherichia_coli');
	genomes.append('neisseria_meningitidis');
	# genomes.append('saccharomyces_cerevisiae');
	# genomes.append('caenorhabditis_elegans');
	# genomes.append('hg19_v38-chr3');

	num_processed_datasets = 0;
	num_datasets = len(machines) * len(genomes);
	
	machine_num = 0;
	for machine_name in machines:
		for reference_name in genomes:
	# for dataset in datasets:
			sys.stderr.write('[%d/%d] Starting machine_name = "%s", reference_name = "%s"...\n' % ((num_processed_datasets + 1), num_datasets, machine_name, reference_name));
			reference_file = '%s/%s.fa' % (REFERENCE_GENOMES_ROOT_ABS, reference_name);
			reads_fastq = '%s/%s/%s/reads.fq' % (READS_SIMULATED_ROOT_ABS, machine_name, reference_name);
			reads_fasta = '%s/%s/%s/reads.fa' % (READS_SIMULATED_ROOT_ABS, machine_name, reference_name);
			output_path = '%s/reads-simulated/%s/%s' % (EVALUATION_PATH_ROOT_ABS, machine_name, reference_name);
			
			if not os.path.exists(output_path):
				sys.stderr.write('Creating folder "%s".\n' % (output_path));
				os.makedirs(output_path);

			###########################################
			###########################################
			machine_suffix = 'v1';
			wrapper_graphmap.run(reads_fastq, reference_file, sam_suffixes[machine_num], output_path, machine_suffix);
			# wrapper_lastal.run(reads_fastq, reference_file, sam_suffixes[machine_num], output_path, machine_suffix);
			# wrapper_blasr.run(reads_fastq, reference_file, sam_suffixes[machine_num], output_path, machine_suffix);
			# wrapper_bwamem.run(reads_fastq, reference_file, sam_suffixes[machine_num], output_path, machine_suffix);
			###########################################
			###########################################

			evalalignments4.EvaluateAlignmentsFromPath(output_path, machine_suffix);

			# sys.stderr.write('[%d/%d] Finished machine_name = "%s", reference_name = "%s"...\n' % ((num_processed_datasets + 1), num_datasets, machine_name, reference_name));
			num_processed_datasets += 1;
			sys.stderr.write('====================================================\n');
				
			# machine_num += 1;

	sys.stderr.write('Done!\n');
	sys.stderr.write('\n');
