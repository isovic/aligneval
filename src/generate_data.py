#! /usr/bin/python

import os;
import glob;
import subprocess;

from basicdefines import *;
import fastqparser;



def peek(fp, num_chars):
	data = fp.read(num_chars);
	if len(data) == 0:
		return '';
	fp.seek(num_chars * -1, 1);
	return data;

def get_single_read(fp):
	lines = '';
	
	line = fp.readline();
	header = line;
	lines += line;
	next_char = peek(fp, 1);
	
	num_lines = 1;
	while len(next_char) > 0 and next_char != lines[0] or (next_char == '@' and num_lines < 4):
		line = fp.readline();
		lines += line;
		next_char = peek(fp, 1);
		num_lines += 1;
		
	return [header.rstrip(), lines.rstrip()];

def interleave(reads1_path, reads2_path, out_path):
	fp1 = open(reads1_path, 'r');
	fp2 = open(reads2_path, 'r');
	fp_out = open(out_path, 'w');

	num_read_pairs = 0;
	while True:
		[header1, read1] = get_single_read(fp1);
		[header2, read2] = get_single_read(fp2);
		
		if (len(read1) == 0 and len(read2) > 0) or (len(read1) > 0 and len(read2) == 0) or (header1[0:-3] != header2[0:-3]):
			print 'ERROR: Reads mismatch! Wrong input files, or reads not in correct order! Read #%d!' % num_read_pairs;
			print 'Header 1: "%s"' % header1;
			print 'Header 2: "%s"' % header2;
			break;
		if (len(read1) == 0 and len(read2) == 0):
			break;
		
		fp_out.write(read1 + '\n');
		fp_out.write(read2 + '\n');
		
		num_read_pairs += 1;
	
	fp1.close();
	fp2.close();
	fp_out.close();

def CreateFolders(machine_name, genome_filename):
	if not os.path.exists(machine_name + '/' + genome_filename):
		print 'Creating %s output folders.' % machine_name;
		os.makedirs(machine_name + '/' + genome_filename);

def EstimateCoverageForNumReads(genome_path, genome_filename, mean_read_length, num_reads):
	fp_in = None;
	complete_genome_path = genome_path + '/' + genome_filename + '.fa'

	try:
		fp_in = open(complete_genome_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % complete_genome_path;
		return;
	
	total_genome_length = 0;
	
	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(header) == 0):
			break;
		
		seq = read[1];
		total_genome_length += len(seq);
		
	fp_in.close();
	
	coverage = int(float(num_reads) / (float(total_genome_length) / float(mean_read_length)));
	
	#print num_reads;
	#print total_genome_length;
	#print mean_read_length;
	#print (float(total_genome_length) / float(mean_read_length))
	#print float(num_reads)
	#print coverage
	
	return coverage;

def ExtractNReadsFromFile(reads_path_prefix, num_reads_to_extract):
	reads_path = reads_path_prefix + '.fq';
	if (os.path.exists(reads_path) == True):
		complete_reads_path = reads_path_prefix + '-complete_dataset.fq';
		#complete_reads_path = os.path.dirname(reads_path) + '/' + os.path.splitext(os.path.basename(reads_path))[0] + '-complete_dataset' + os.path.splitext(reads_path)[1];
		
		print 'Renaming file "%s" to "%s"...' % (reads_path, complete_reads_path);
		
		os.rename(reads_path, complete_reads_path);
		
		fp_in = open(complete_reads_path, 'r');
		fp_out = open(reads_path, 'w');

		num_read_pairs = 0;
		i = 0;
		while i < num_reads_to_extract:
			[header, read] = get_single_read(fp_in);
			
			if (len(read) == 0):
				break;

			fp_out.write(read + '\n');
			
			i += 1;
		
		fp_in.close();
		fp_out.close();
	else:
		print 'ERROR: Reads file "%s" does not exist!' % (reads_path);

def ExtractNAlignmentsFromSAM(reads_path_prefix, num_alignments_to_extract):
	sam_path = reads_path_prefix + '.sam';
	if (os.path.exists(sam_path) == True):
		complete_sam_path = reads_path_prefix + '-complete_dataset.sam';
		#complete_reads_path = os.path.dirname(reads_path) + '/' + os.path.splitext(os.path.basename(reads_path))[0] + '-complete_dataset' + os.path.splitext(reads_path)[1];
		
		print 'Renaming file "%s" to "%s"...' % (sam_path, complete_sam_path);
		
		os.rename(sam_path, complete_sam_path);
		
		fp_in = open(complete_sam_path, 'r');
		fp_out = open(sam_path, 'w');
		
		current_num_alignments = 0;
		
		print 'num_alignments_to_extract = %d' % num_alignments_to_extract;

		for line in fp_in:
			fp_out.write(line);
			
			if (line.startswith('@') == False and len(line.strip()) > 0):
				current_num_alignments += 1;
			
			if (current_num_alignments >= num_alignments_to_extract):
				break;
		
		print 'current_num_alignments = %d' % current_num_alignments;

		fp_in.close();
		fp_out.close();
	else:
		print 'ERROR: Reads file "%s" does not exist!' % (reads_path);

def Generate454(genome_path, genome_filename, fold_coverage=10, mean_fragsize=1500, std_fragsize=10, machine_name='Roche454', num_reads_to_generate=-1):
	complete_genome_path = genome_path + '/' + genome_filename + '.fa'
	simulator_bin = TOOLS_ROOT_ABS + '/art_bin_VanillaIceCream/art_454';
#	out_file_prefix = machine_name + '/' + genome_filename + '/' + genome_filename;

	is_paired_end = False;
	
	if mean_fragsize==0 or std_fragsize==0:
		machine_name += '-single_end';
		CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
		CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
		
		out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
		#shell_command = simulator_bin + r' -s -r 1403002416 ' + complete_genome_path + r' ' + out_file_prefix + r' ' + str(fold_coverage);
		shell_command = simulator_bin + r' -s ' + complete_genome_path + r' ' + out_file_prefix + r' ' + str(fold_coverage);
		is_paired_end = False;
	else:
		machine_name += '-paired_end';
		CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
		CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
		
		out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
		#shell_command = simulator_bin + r' -s -r 1403002416 ' + complete_genome_path + r' ' + out_file_prefix + r' ' + str(fold_coverage) + r' ' + str(mean_fragsize) + r' ' + str(std_fragsize);
		shell_command = simulator_bin + r' -s ' + complete_genome_path + r' ' + out_file_prefix + r' ' + str(fold_coverage) + r' ' + str(mean_fragsize) + r' ' + str(std_fragsize);
		is_paired_end = True;
	
	print 'Executing command: "%s"' % shell_command;	
	subprocess.call(shell_command, shell=True);

	# if (GENERATE_FIXED_AMOUNT_OF_READS == True):
	# 	ExtractNReadsFromFile(out_file_prefix, NUM_READS_TO_GENERATE);
	# 	ExtractNAlignmentsFromSAM(out_file_prefix, NUM_READS_TO_GENERATE);
	if (num_reads_to_generate > 0):
		ExtractNReadsFromFile(out_file_prefix, num_reads_to_generate);
		ExtractNAlignmentsFromSAM(out_file_prefix, num_reads_to_generate);
	
	if is_paired_end == True:
		print 'Interleaving paired end reads to file "%s"' % (out_file_prefix + '.fq');
		interleave(out_file_prefix + '1.fq', out_file_prefix + '2.fq', out_file_prefix + '.fq');

def GenerateIllumina(genome_path, genome_filename, read_length=100, fold_coverage=10, mean_fragsize=500, std_fragsize=10, machine_name='Illumina', num_reads_to_generate=-1):
	complete_genome_path = genome_path + '/' + genome_filename + '.fa'
	simulator_bin = TOOLS_ROOT_ABS + '/art_bin_VanillaIceCream/art_illumina';
	
	is_paired_end = False;
	
	if mean_fragsize==0 or std_fragsize==0:			# Single end reads
		machine_name += '-single_end';
		CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
		CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
		
		out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
		#shell_command = simulator_bin + r' --rndSeed 1403000281 -i ' + complete_genome_path + r' -o ' + out_file_prefix + r' -l ' + str(read_length) + r' -f ' + str(fold_coverage) + r' -sam';	# Single-end
		shell_command = simulator_bin + r' -i ' + complete_genome_path + r' -o ' + out_file_prefix + r' -l ' + str(read_length) + r' -f ' + str(fold_coverage) + r' -sam';	# Single-end
		is_paired_end = False;
	else:			# Paired end reads
		machine_name += '-paired_end';
		CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
		CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
		
		out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
		#shell_command = simulator_bin + r' --rndSeed 1403000281 -i ' + complete_genome_path + r' -o ' + out_file_prefix + r' -l ' + str(read_length) + r' -f ' + str(fold_coverage) + r' -m ' + str(mean_fragsize) + r' -s ' + str(std_fragsize) + r' -sam';	# Paired-end
		shell_command = simulator_bin + r' -i ' + complete_genome_path + r' -o ' + out_file_prefix + r' -l ' + str(read_length) + r' -f ' + str(fold_coverage) + r' -m ' + str(mean_fragsize) + r' -s ' + str(std_fragsize) + r' -sam';	# Paired-end
		is_paired_end = True;

	print 'Executing command: "%s"' % shell_command;
	subprocess.call(shell_command, shell=True);
	
	# if (GENERATE_FIXED_AMOUNT_OF_READS == True):
	# 	ExtractNReadsFromFile(out_file_prefix, NUM_READS_TO_GENERATE);
	# 	ExtractNAlignmentsFromSAM(out_file_prefix, NUM_READS_TO_GENERATE);
	if (num_reads_to_generate > 0):
		ExtractNReadsFromFile(out_file_prefix, num_reads_to_generate);
		ExtractNAlignmentsFromSAM(out_file_prefix, num_reads_to_generate);
	
	if is_paired_end == True:
		print 'Interleaving paired end reads to file "%s"' % (out_file_prefix + '.fq');
		interleave(out_file_prefix + '1.fq', out_file_prefix + '2.fq', out_file_prefix + '.fq');
	


def GetPBSimRefName(ref_file):
	try:
		fp = open(ref_file, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % ref_file;
	
	header = fp.readline();
	header = header.rstrip();
	fp.close();
	
	if not header.startswith('>'):
		print "ERROR: PBsim's ref file does not start with a FASTA header!";
	
	ref_name = header[1:];
	trimmed_ref_name = ref_name.split()[0];
	
	return [ref_name, trimmed_ref_name];

# --difference-ratio   ratio of differences. substitution:insertion:deletion.
def GeneratePacBio(genome_path, genome_filename, fold_coverage=20, length_mean=3000, length_sd=2300.0, length_min=100, length_max=25000, accuracy_mean=0.78, accuracy_sd=0.02, accuracy_min=0.75, difference_ratio='10:60:30', machine_name='PacBio', num_reads_to_generate=-1):
#	machine_name = 'PacBio';
	CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
	CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
	
	complete_genome_path = genome_path + '/' + genome_filename + '.fa'
	out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
	
	simulator_path = TOOLS_ROOT_ABS + '/pbsim-1.0.3-Linux-amd64';
	simulator_bin = TOOLS_ROOT_ABS + '/pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim';
	
	final_sam_file = out_file_prefix + '.sam';
	fp = open(final_sam_file, 'w');
	fp.close();
	
	final_fastq_file = out_file_prefix + '.fq';
	fp = open(final_fastq_file, 'w');
	fp.close();
# SIMULATOR_PATH=/home/ivan/work/eclipse-workspace/golden-bundle/tools/pbsim-1.0.3-Linux-amd64
# REFERENCE_PATH=/home/ivan/work/eclipse-workspace/golden-bundle/reference-genomes/escherichia_coli.fa
# $SIMULATOR_PATH/Linux-amd64/bin/pbsim --data-type CLR --depth 20 --model_qc $SIMULATOR_PATH/data/model_qc_clr --seed 1234567890 --prefix $OUT_FILE_PREFIX $REFERENCE_PATH
# /last-460/scripts/maf-convert.py sam $MAF_FILE $SAM_FILE

	#random_seed = '1234567890';
	random_seed = '32874638';

	# Data type:
	#	Continuous Long Read (CLR) : long and high error rate.
	#	Circular consensus Read (CCS) : short and low error rate.
	# --difference-ratio   ratio of differences. substitution:insertion:deletion.
	shell_command = simulator_bin + r' --data-type CLR --depth ' + str(fold_coverage) + \
					' --model_qc ' + simulator_path + \
					'/data/model_qc_clr ' + \
					' --length-mean ' + str(length_mean) + \
					' --length-sd ' + str(length_sd) + \
					' --length-min ' + str(length_min) + \
					' --length-max ' + str(length_max) + \
					' --accuracy-mean ' + str(accuracy_mean) + \
					' --accuracy-sd ' + str(accuracy_sd) + \
					' --accuracy-min ' + str(accuracy_min) + \
					' --difference-ratio ' + difference_ratio + \
					' --prefix ' + out_file_prefix + \
					r' ' + complete_genome_path;
				
					#' --seed ' + random_seed + \
	
	print 'Simulating PacBio reads using PBsim';
	print 'Executing command: "%s"' % shell_command;
	
	#exit(1);
	subprocess.call(shell_command, shell=True);

	print ' ';
	
	print 'Converting generated *.maf files to corresponding SAM files.';

	maf_files = glob.glob(out_file_prefix + '*.maf');
	maf_files = sorted(maf_files);
	print maf_files;
	
	for maf_file in maf_files:		
		# Convert the maf file to SAM format.
		sam_file = maf_file[0:-3] + 'sam';
		shell_command = TOOLS_ROOT_ABS + '/last-460/scripts/maf-convert.py sam ' + maf_file + ' > ' + sam_file;
		print 'Converting MAF to SAM ("%s" -> "%s")' % (maf_file, sam_file);
		print 'Executing command: "%s"' % shell_command;
		subprocess.call(shell_command, shell=True);
		
		# Use Bash's sed to replace the 'ref' keyword in the generated maf files with the actual name of the reference sequence, and concatenate all the separate SAM files to one file.
		[reference_name, trimmed_reference_name] = GetPBSimRefName(maf_file[0:-3] + 'ref');
		# Here we escape the special characters so that SED command runs properly if any of these characters should to appear in a FASTA header.
		escape_chars = r'\/()[].*^$';
		reference_name = ''.join([('\\' + char) if char in escape_chars else char for char in reference_name]);
		shell_command = r'cat ' + sam_file + r" | sed 's/^\(.*\)ref/\1" + reference_name + r"/' >> " + final_sam_file
		print 'Replacing PBsim\'s "ref" keyword with actual FASTA header';
		print 'Executing command: "%s"' % shell_command;
		subprocess.call(shell_command, shell=True);
		
		fastq_file = maf_file[0:-3] + 'fastq';
		shell_command = r'cat ' + fastq_file + ' >> ' + final_fastq_file;
		print 'Concatenating FASTQ file to the total reads file ("%s" -> "%s")' % (fastq_file, final_fastq_file);
		print 'Executing command: "%s"' % shell_command;
		subprocess.call(shell_command, shell=True);
		print ' ';
		
		print ' ';

#	if (GENERATE_FIXED_AMOUNT_OF_READS == True):
	if (num_reads_to_generate > 0):
		ExtractNReadsFromFile(out_file_prefix, num_reads_to_generate);
		ExtractNAlignmentsFromSAM(out_file_prefix, num_reads_to_generate);

	#sam_files = glob.glob(out_file_prefix + '*.sam');
	#sam_files = sorted(sam_files);
	#print sam_files;

def GenerateOxfordNanoporeFromObservedStatistics(genome_filename, num_reads_to_generate=-1):
	if (num_reads_to_generate <= 0):
		coverage = 20;
		machine_suffix = '-cov20';
	else:
		num_reads_to_generate = 10000;
		mean_read_length = 5400;
		coverage = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, genome_filename, mean_read_length, num_reads_to_generate) + 1;
		machine_suffix = '-%dk' % (num_reads_to_generate / 1000);

	## The maximum value for length_max parameter (limited by PBsim) is 100000.

	print 'num_reads_to_generate = %d' % num_reads_to_generate;
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4328.42, length_sd=3888.69, length_min=10, length_max=64749.00,
	# 															accuracy_mean=(1.0 - 0.40), accuracy_sd=0.07, accuracy_min=(0.05), difference_ratio='42:18:40',
	# 															machine_name='OxfordNanopore-pbsim-graphmap-observed_mode1' + machine_suffix, num_reads_to_generate=3200);

	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4328.42, length_sd=3888.69, length_min=10, length_max=64749.00,
	# 															accuracy_mean=(1.0 - 0.50), accuracy_sd=0.02, accuracy_min=(0.05), difference_ratio='44:35:21',
	# 															machine_name='OxfordNanopore-pbsim-graphmap-observed_mode2' + machine_suffix, num_reads_to_generate=6800);

	# These are simulations for difference_ratio obtained from LAST's alignments:
	# 1d reads:
	# [CIGAR statistics - individual indels]
	#                       	mean	std	median	min	max
	# Error rate stats:     	0.41	0.05	0.40	0.10	0.60
	# Insertion rate stats: 	0.05	0.02	0.05	0.00	0.23
	# Deletion rate stats:  	0.16	0.05	0.15	0.00	0.49
	# Mismatch rate stats:  	0.20	0.03	0.20	0.03	0.32
	# Match rate stats:     	0.75	0.04	0.75	0.59	0.97
	# Read length stats:    	3629.76	3294.04	2438.00	57.00	31299.00
	# Difference ratio: 51:11:38 (mismatch:insertion:deletion)
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4200, length_sd=3300, length_min=50, length_max=64000,
	# 															accuracy_mean=(1.0 - 0.41), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio='51:11:38',
	# 															machine_name='OxfordNanopore-pbsim-observed_last-1d' + machine_suffix, num_reads_to_generate=10000);
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
																accuracy_mean=(1.0 - 0.40), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio='51:11:38',
																machine_name='OxfordNanopore-pbsim-observed_last-1d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	# 2d reads:
	# [CIGAR statistics - individual indels]
	#                       	mean	std	median	min	max
	# Error rate stats:     	0.31	0.09	0.30	0.02	0.59
	# Insertion rate stats: 	0.05	0.03	0.05	0.00	0.28
	# Deletion rate stats:  	0.09	0.06	0.08	0.00	0.48
	# Mismatch rate stats:  	0.16	0.07	0.16	0.00	0.36
	# Match rate stats:     	0.78	0.07	0.79	0.59	0.99
	# Read length stats:    	2006.14	3015.25	614.00	42.00	28601.00
	# Difference ratio: 55:17:28 (mismatch:insertion:deletion)
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=2000, length_sd=3300, length_min=50, length_max=64000,
	# 															accuracy_mean=(1.0 - 0.31), accuracy_sd=0.09, accuracy_min=(1.0 - 0.59), difference_ratio='55:17:28',
	# 															machine_name='OxfordNanopore-pbsim-observed_last-2d' + machine_suffix, num_reads_to_generate=10000);
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
																accuracy_mean=(1.0 - 0.20), accuracy_sd=0.05, accuracy_min=(1.0 - 0.50), difference_ratio='55:17:28',
																machine_name='OxfordNanopore-pbsim-observed_last-2d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);



	# These are simulations for difference_ratio obtained from GraphMap's alignments:
	# 1d reads:
	# [CIGAR statistics - individual indels]
	#                       	mean	std	median	min	max
	# Error rate stats:     	0.46	0.06	0.49	0.26	0.95
	# Insertion rate stats: 	0.14	0.07	0.12	0.00	0.30
	# Deletion rate stats:  	0.13	0.05	0.13	0.00	0.60
	# Mismatch rate stats:  	0.19	0.03	0.19	0.04	0.70
	# Match rate stats:     	0.67	0.09	0.66	0.23	0.95
	# Read length stats:    	4328.42	3888.69	3361.00	50.00	64749.00
	# Difference ratio: 44:29:27 (mismatch:insertion:deletion)
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
	# 															accuracy_mean=(1.0 - 0.46), accuracy_sd=0.06, accuracy_min=(1.0 - 0.60), difference_ratio='44:29:27',
	# 															machine_name='OxfordNanopore-pbsim-observed_graphmap-1d' + machine_suffix, num_reads_to_generate=10000);
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
																accuracy_mean=(1.0 - 0.40), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio='44:29:27',
																machine_name='OxfordNanopore-pbsim-observed_graphmap-1d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);

	# 2d reads:
	# [CIGAR statistics - individual indels]
	#                       	mean	std	median	min	max
	# Error rate stats:     	0.26	0.07	0.24	0.13	0.45
	# Insertion rate stats: 	0.06	0.02	0.06	0.02	0.25
	# Deletion rate stats:  	0.11	0.05	0.09	0.03	0.30
	# Mismatch rate stats:  	0.09	0.03	0.09	0.03	0.22
	# Match rate stats:     	0.85	0.04	0.85	0.61	0.93
	# Read length stats:    	6167.24	3532.11	6077.50	158.00	30804.00
	# Difference ratio: 37:23:40 (mismatch:insertion:deletion)
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
	# 															accuracy_mean=(1.0 - 0.26), accuracy_sd=0.07, accuracy_min=(1.0 - 0.45), difference_ratio='37:23:40',
	# 															machine_name='OxfordNanopore-pbsim-observed_graphmap-2d' + machine_suffix, num_reads_to_generate=10000);
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
																accuracy_mean=(1.0 - 0.20), accuracy_sd=0.05, accuracy_min=(1.0 - 0.50), difference_ratio='37:23:40',
																machine_name='OxfordNanopore-pbsim-observed_graphmap-2d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);

	

def GenerateNGSData(num_reads_to_generate=-1):
	##### NGS DATA ###
	if (num_reads_to_generate <= 0):
		coverage_nmeni = 20;
		coverage_ecoli = 20;
		coverage_scerevisiae = 20;
		coverage_celegans = 20;
		coverage_hg19v38chr3 = 20;
		machine_suffix = '-cov20';
	else:
		mean_read_length = 100;
		coverage_nmeni = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for neisseria_meningitidis: %d' % coverage_nmeni;
		coverage_ecoli = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for escherichia_coli: %d' % coverage_ecoli;
		coverage_scerevisiae = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for saccharomyces_cerevisiae: %d' % coverage_scerevisiae;
		coverage_celegans = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for caenorhabditis_elegans: %d' % coverage_celegans;
		coverage_hg19v38chr3 = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for hg19_v38-chr3: %d' % coverage_hg19v38chr3;
		machine_suffix = '-%dk' % (num_reads_to_generate / 1000);

	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', read_length=100, fold_coverage=coverage_nmeni, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', read_length=100, fold_coverage=coverage_ecoli, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', read_length=100, fold_coverage=coverage_scerevisiae, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', read_length=100, fold_coverage=coverage_celegans, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', read_length=100, fold_coverage=coverage_hg19v38chr3, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', fold_coverage=coverage_nmeni, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', fold_coverage=coverage_ecoli, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', fold_coverage=coverage_scerevisiae, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', fold_coverage=coverage_celegans, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', fold_coverage=coverage_hg19v38chr3, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	##################

def GeneratePacBioData(num_reads_to_generate=-1):
	##### PACBIO NORMAL DATA #####
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', fold_coverage=80, length_mean=3000, length_sd=2300.0, accuracy_mean=0.78, accuracy_sd=0.02, accuracy_min=0.75, machine_name='PacBio');
	if (num_reads_to_generate <= 0):
		coverage_nmeni = 20;
		coverage_ecoli = 20;
		coverage_scerevisiae = 20;
		coverage_celegans = 20;
		coverage_hg19v38chr3 = 20;
		machine_suffix = '-cov20';
	else:
		mean_read_length = 3000;
		coverage_nmeni = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for neisseria_meningitidis: %d' % coverage_nmeni;
		coverage_ecoli = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for escherichia_coli: %d' % coverage_ecoli;
		coverage_scerevisiae = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for saccharomyces_cerevisiae: %d' % coverage_scerevisiae;
		coverage_celegans = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for caenorhabditis_elegans: %d' % coverage_celegans;
		coverage_hg19v38chr3 = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', mean_read_length, num_reads_to_generate) + 1;
		print 'Coverage for hg19_v38-chr3: %d' % coverage_hg19v38chr3;
		machine_suffix = '-%dk' % (num_reads_to_generate / 1000);
		
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', fold_coverage=coverage_nmeni, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', fold_coverage=coverage_ecoli, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', fold_coverage=coverage_scerevisiae, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', fold_coverage=coverage_celegans, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', fold_coverage=coverage_hg19v38chr3, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	##############################

def GenerateOxfordNanoporeData(num_reads_to_generate=-1):
	##### OXFORD NANOPORE DATA #####
	# --difference-ratio   ratio of differences. substitution:insertion:deletion.
	GenerateOxfordNanoporeFromObservedStatistics('neisseria_meningitidis', num_reads_to_generate=num_reads_to_generate);
	GenerateOxfordNanoporeFromObservedStatistics('escherichia_coli', num_reads_to_generate=num_reads_to_generate);
	GenerateOxfordNanoporeFromObservedStatistics('saccharomyces_cerevisiae', num_reads_to_generate=num_reads_to_generate);
	GenerateOxfordNanoporeFromObservedStatistics('caenorhabditis_elegans', num_reads_to_generate=num_reads_to_generate);
	GenerateOxfordNanoporeFromObservedStatistics('hg19_v38-chr3', num_reads_to_generate=num_reads_to_generate);
	


def GenerateAll():
	num_reads_to_generate = 10000;

	GenerateNGSData(num_reads_to_generate);
	GeneratePacBioData(num_reads_to_generate);
	GenerateOxfordNanoporeData(num_reads_to_generate);



def Main():
	GenerateAll();

if __name__ == '__main__':
	Main();
