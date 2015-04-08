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

def register_scores(simulated_dataset, reference_name, eval_scores, precision_index, recall_index, ret_results_dataset_header, ret_results_genome_header, ret_results_table):
	ret_results_dataset_header.append(simulated_dataset);
	ret_results_genome_header.append(reference_name);
	if (eval_scores != None and len(eval_scores.keys()) > 0):
		# Add the new scores to the existing table.
		for evaluated_mapper in sorted(eval_scores.keys()):
			# If a mapper has not been evaluated before, fill the table with '-' characters.
			if ((evaluated_mapper in ret_results_table) == False):
				ret_results_table[evaluated_mapper] = ['-'] * (len(ret_results_dataset_header) - 2);
			eval_score = eval_scores[evaluated_mapper];
			ret_results_table[evaluated_mapper].append('%.2f / %.2f' % (eval_score[precision_index], eval_score[recall_index]));
		# Check if there were mappers in other datasets which were not evaluated on this dataset, and fill the current row with '-' characters.
		for previously_evaluated_mapper in sorted(ret_results_table.keys()):
			if ((previously_evaluated_mapper in eval_scores.keys()) == False):
				ret_results_table[previously_evaluated_mapper].append('-');
	else:
		# In this case, none of the mappers was evalated on this dataset. Fill the entire row with '-' characters.
		for previously_evaluated_mapper in sorted(ret_results_table.keys()):
			ret_results_table[previously_evaluated_mapper].append('-');

def convert_results_table(results_dataset_header, results_genome_header, results_table):
	table = [];
	row = [results_dataset_header[0], results_genome_header[0]] + sorted(results_table.keys());
	table.append(row);
	i = 1;
	while (i < len(results_dataset_header)):
		dataset_header = results_dataset_header[i];
		genome_header = results_genome_header[i];
		row = [dataset_header, genome_header];
		for evaluated_mapper in sorted(results_table.keys()):
			row.append(results_table[evaluated_mapper][i - 1]);			# -1 because we are using the results_dataset_header index which has an additional entry at the beginning (namely 'Mapper').
		table.append(row);
		i += 1;
	return table;

def filter_only_select_mappers(filter_mappers, eval_scores):
	filtered_eval_scores = {};
	for mapper in filter_mappers:
		match_keys = [key for key in eval_scores.keys() if (mapper.lower() in key.lower())];
		for match_key in match_keys:
			filtered_eval_scores[match_key] = eval_scores[match_key];
	return filtered_eval_scores;



def write_table(fp, table):
	for row in table:
		fp.write('\t'.join(row) + '\n');



if __name__ == "__main__":
	aligner_wrappers = find_files(WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py');
	for wrapper in aligner_wrappers:
		wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0];
		command = 'import %s' % (wrapper_basename);
		exec(command);

	machine_names = [];
	simulated_datasets = [];

	simulated_datasets.append('Illumina-1k-single_end');								machine_names.append('illumina');
	simulated_datasets.append('PacBio-1k');												machine_names.append('pacbio');
	simulated_datasets.append('OxfordNanopore-pbsim-observed_last-2d-1k');				machine_names.append('nanopore');
	simulated_datasets.append('OxfordNanopore-pbsim-observed_graphmap-2d-1k');			machine_names.append('nanopore');
	simulated_datasets.append('OxfordNanopore-pbsim-observed_last-1d-1k');				machine_names.append('nanopore');
	simulated_datasets.append('OxfordNanopore-pbsim-observed_graphmap-1d-1k');			machine_names.append('nanopore');
	simulated_datasets.append('OxfordNanopore-pbsim-observed_marginalign-2d-1k');		machine_names.append('nanopore');

	genomes = [];
	genomes.append('neisseria_meningitidis');
	genomes.append('escherichia_coli');
	genomes.append('saccharomyces_cerevisiae');
	genomes.append('caenorhabditis_elegans');
	genomes.append('hg19_v38-chr3');

	num_processed_datasets = 0;
	num_datasets = len(simulated_datasets) * len(genomes);

	bp_dist = 50;

	results_bp_dataset_header = ['Dataset'];	results_bp_genome_header = ['Genome'];	results_bp_table = {};
	results_cb_dataset_header = ['Dataset'];	results_cb_genome_header = ['Genome'];	results_cb_table = {};

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

			# eval_scores = evalalignments4.EvaluateAlignmentsFromPath(output_path, machine_suffix, bp_dist=bp_dist);
			eval_scores = evalalignments4.EvaluateAlignmentsFromPath(output_path, 'GraphMap-v1', bp_dist=bp_dist);
			print eval_scores;
			print '';

			### Use this to select only requeired mappers for reporting the results:
			# eval_scores = filter_only_select_mappers(['GraphMap'], eval_scores);

			# Precision and recall for mapping position, allowed within +-bp distance from expected location.
			register_scores(simulated_dataset, reference_name, eval_scores, 0, 1, results_bp_dataset_header, results_bp_genome_header, results_bp_table);

			# Precision and recall for correctly called bases.
			register_scores(simulated_dataset, reference_name, eval_scores, 2, 3, results_cb_dataset_header, results_cb_genome_header, results_cb_table);



			sys.stderr.write('[%d/%d] Finished simulated_dataset = "%s", reference_name = "%s"...\n' % ((num_processed_datasets + 1), num_datasets, simulated_dataset, reference_name));
			num_processed_datasets += 1;
			sys.stderr.write('====================================================\n');
			
		machine_num += 1;

	table_bp = convert_results_table(results_bp_dataset_header, results_bp_genome_header, results_bp_table);
	table_cb = convert_results_table(results_cb_dataset_header, results_cb_genome_header, results_cb_table);
	sys.stdout.write('Precision / Recall for mapping position to within +-bp distance:\n');

	write_table(sys.stdout, table_bp);
	sys.stdout.write('\n');
	sys.stdout.write('Precision / Recall for number of correctly mapped bases:\n');
	write_table(sys.stdout, table_cb);

	### Writing the results to files.
	if not os.path.exists(RESULTS_PATH_ROOT_ABS):
		sys.stderr.write('Creating folder "%s".\n' % (RESULTS_PATH_ROOT_ABS));
		os.makedirs(RESULTS_PATH_ROOT_ABS);

	out_path_bp = '%s/precision_recall-%dbp_distance.csv' % (RESULTS_PATH_ROOT_ABS, bp_dist);
	fp_out_bp = open(out_path_bp, 'w');
	write_table(fp_out_bp, table_bp);
	fp_out_bp.close();

	out_path_cb = '%s/precision_recall-correct_bases.csv' % (RESULTS_PATH_ROOT_ABS);
	fp_out_cb = open(out_path_cb, 'w');
	write_table(fp_out_cb, table_cb);
	fp_out_cb.close();



	sys.stdout.write('\n');
	sys.stdout.write('\n');
	sys.stdout.write('\n');
	sys.stderr.write('Done!\n');
	sys.stderr.write('\n');
