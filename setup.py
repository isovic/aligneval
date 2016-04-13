#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys
sys.path.append(SCRIPT_PATH + '/src')
sys.path.append(SCRIPT_PATH + '/wrappers');

import subprocess;

# from evalpaths import *
# import generate_data
# import filesandfolders;
# from basicdefines import *;
import basicdefines;
import generate_data;



def create_folders():
	sys.stderr.write('Generating folders...\n');

	if not os.path.exists(basicdefines.EVALUATION_PATH_ROOT_ABS):
		sys.stderr.write('Creating folder "%s".\n' % basicdefines.EVALUATION_PATH_ROOT_ABS);
		os.makedirs(basicdefines.EVALUATION_PATH_ROOT_ABS);

	if not os.path.exists(basicdefines.TOOLS_ROOT_ABS):
		sys.stderr.write('Creating folder "%s".\n' % basicdefines.TOOLS_ROOT_ABS);
		os.makedirs(basicdefines.TOOLS_ROOT_ABS);

	if not os.path.exists(basicdefines.ALIGNERS_PATH_ROOT_ABS):
		sys.stderr.write('Creating folder "%s".\n' % basicdefines.ALIGNERS_PATH_ROOT_ABS);
		os.makedirs(basicdefines.ALIGNERS_PATH_ROOT_ABS);

	sys.stderr.write('\n');

def unpack_reference_genomes():
	sys.stderr.write('Unpacking reference genomes [~400 MB]\n');
	tar_gz_references = basicdefines.find_files(basicdefines.REFERENCE_GENOMES_ROOT_ABS, '*.tar.gz');
	for file_path in tar_gz_references:
		subprocess.call(('tar -xzvf %s -C %s/' % (file_path, os.path.dirname(file_path))), shell='True');
	sys.stderr.write('\n');

def download_hg19_GRCh38_reference():
	hg19_GRCh38_path = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz';
	sys.stderr.write('Downloading the Human reference genome (GRCh38).');
	command = 'mkdir -p %s/download; cd %s/download; wget %s' % (basicdefines.REFERENCE_GENOMES_ROOT_ABS, basicdefines.REFERENCE_GENOMES_ROOT_ABS, hg19_GRCh38_path);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	hg_archive_filename = os.path.basename(hg19_GRCh38_path);
	sys.stderr.write('Unpacking the Human reference genome (GRCh38).');
	command = 'cd %s/download; gunzip %s' % (basicdefines.REFERENCE_GENOMES_ROOT_ABS, hg_archive_filename);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	hg_archive_filename = os.path.basename(hg19_GRCh38_path);
	sys.stderr.write('Unpacking the Human reference genome (GRCh38).');
	command = 'mv %s/download/hg38.fa %s/' % (basicdefines.REFERENCE_GENOMES_ROOT_ABS, basicdefines.REFERENCE_GENOMES_ROOT_ABS);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

def download_hg19_GRCh37_reference():
	hg19_GRCh37_path = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz';
	sys.stderr.write('Downloading the Human reference genome (GRCh37).\n');
	command = 'mkdir -p %s/download; cd %s/download; wget %s' % (basicdefines.REFERENCE_GENOMES_ROOT_ABS, basicdefines.REFERENCE_GENOMES_ROOT_ABS, hg19_GRCh37_path);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	hg_archive_filename = os.path.basename(hg19_GRCh37_path);
	sys.stderr.write('Unpacking the Human reference genome (GRCh37).\n');
	command = 'cd %s/download; mkdir -p hgchroms; tar -C hgchroms/ -xzvf %s' % (basicdefines.REFERENCE_GENOMES_ROOT_ABS, hg_archive_filename);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');
	
	hg19_with_masking_fa = '%s/download/hgchroms/hg19_with_masking.fa' % (basicdefines.REFERENCE_GENOMES_ROOT_ABS);
	sys.stderr.write('Joining the chromosomes into one multifasta file on path: %s.\n' % (hg19_with_masking_fa));
	command = 'cd %s/download;' % (basicdefines.REFERENCE_GENOMES_ROOT_ABS);
	command += 'cat hgchroms/chr1.fa > %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr2.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr3.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr4.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr5.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr6.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr7.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr8.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr9.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr10.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr11.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr12.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr13.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr14.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr15.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr16.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr17.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr18.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr19.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr20.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr21.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chr22.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chrX.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chrY.fa >> %s;' % (hg19_with_masking_fa);
	command += 'cat hgchroms/chrM.fa >> %s;' % (hg19_with_masking_fa);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	hg19_fa = '%s/hg19.fa' % (basicdefines.REFERENCE_GENOMES_ROOT_ABS);
	sys.stderr.write('Converting the hg19 multifasta (%s) bases to uppercase (%s).\n' % (hg19_with_masking_fa, hg19_fa));
	command = 'cd %s; %s/samscripts/src/fastqfilter.py touppercase %s %s' % (basicdefines.REFERENCE_GENOMES_ROOT_ABS, basicdefines.TOOLS_ROOT_ABS, hg19_with_masking_fa, hg19_fa);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

def unpack_sample_sim_reads():
	sys.stderr.write('Unpacking pre-simulated data [~400 MB]\n');
	tar_gz = basicdefines.find_files(basicdefines.READS_SIMULATED_ROOT_ABS, '*.tar.gz');
	for file_path in tar_gz:
		subprocess.call(('tar -xzvf %s -C %s/' % (file_path, os.path.dirname(file_path))), shell='True');
	sys.stderr.write('\n');

def download_aligners():
	sys.stderr.write('Installing alignment tools.\n');

	aligner_wrappers = basicdefines.find_files(basicdefines.WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py');

	for wrapper in aligner_wrappers:
		wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0];
		command = 'import %s; %s.download_and_install()' % (wrapper_basename, wrapper_basename);
		exec(command);

def setup_tools():
	if (not os.path.exists(basicdefines.TOOLS_ROOT_ABS)):
		os.makedirs(basicdefines.TOOLS_ROOT_ABS);

	sys.stderr.write('Cloning Cgmemtime Git repo. Git needs to be installed.\n');
	command = 'cd %s; git clone https://github.com/isovic/cgmemtime.git' % (basicdefines.TOOLS_ROOT_ABS);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('Cloning Samscripts Git repo. Git needs to be installed.\n');
	command = 'cd %s; git clone https://github.com/isovic/samscripts.git' % (basicdefines.TOOLS_ROOT_ABS);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('Downloading and unpacking the ART next generation sequence simulator.\n');
	command = 'cd %s; wget http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz; tar -xzvf artbinvanillaicecream031114linux64tgz.tgz' % (basicdefines.TOOLS_ROOT_ABS);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('Downloading and unpacking PBsim.\n');
	command = 'cd %s; wget http://pbsim.googlecode.com/files/pbsim-1.0.3-Linux-amd64.tar.gz; tar -xzvf pbsim-1.0.3-Linux-amd64.tar.gz' % (basicdefines.TOOLS_ROOT_ABS);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('Downloading and unpacking LAST aligner. Its scripts are needed to convert from MAF to SAM.\n');
	command = 'cd %s; wget http://last.cbrc.jp/last-534.zip; unzip last-534.zip' % (basicdefines.TOOLS_ROOT_ABS);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');



def verbose_usage_and_exit():
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s [mode]\n' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\tParameter mode specifies which step to execute.\n');
	sys.stderr.write('\t- mode - "all", "folders", "references", "aligners", "tools", "simdata" or "generate-simdata".\n');
	sys.stderr.write('\n');
	exit(0);

if __name__ == '__main__':
	if (len(sys.argv) != 2):
		verbose_usage_and_exit();

	mode = sys.argv[1];
	mode_valid = False;

	if (mode == 'references' or mode == 'simdata'):
		mode_valid = True;
		sys.stderr.write('Running this script will consume large amount of disk space.\n');
		yes_no = raw_input("Do you want to continue? [y/n] ");
		if (yes_no != 'y'):
			sys.stderr.write('Exiting.\n\n');
			exit(0);

	if (mode == 'all' or mode == 'folders'):
		create_folders();
		mode_valid = True;

	if (mode == 'all' or mode == 'references'):
		unpack_reference_genomes();
#		download_hg19_GRCh38_reference();
		sys.stderr.write('Please make sure you ran "./setup tools" prior to running this command.\n\n');
		download_hg19_GRCh37_reference();
		mode_valid = True;

	if (mode == 'all' or mode == 'aligners'):
		download_aligners();
		mode_valid = True;

	if (mode == 'all' or mode == 'tools'):
		setup_tools();
		mode_valid = True;

	if (mode == 'all' or mode == 'simdata'):
		# generate_data.GenerateAll();
		sys.stderr.write('Please make sure you ran "./setup tools" and "./setup references"  prior to running this command.\n\n');
		unpack_sample_sim_reads();
		if (not os.path.exists('%s/saccharomyces_cerevisiae.fa' % basicdefines.REFERENCE_GENOMES_ROOT_ABS)):
			sys.stderr.write('ERROR: Can not continue with setting up the simulated data until reference sequences are unpacked! Run "./setup.py references" first.\n');
			exit(1);
		generate_data.GenerateGridTest(10000);
		mode_valid = True;

	if (mode == 'generate-simdata'):
		sys.stderr.write('Please make sure you ran "./setup tools" and "./setup references"  prior to running this command.\n\n');
		generate_data.GenerateAll();
		mode_valid = True;

	if (mode_valid == False):
		sys.stderr.write('Selected mode not recognized!\n');
		verbose_usage_and_exit();



# sys.stderr.write('Generating simulated data...\n');
# generate_data.GenerateAll();
# sys.stderr.write('Finished generating simulated data!\n');



# sudo apt-get install python-matplotlib
# sudo apt-get install python-tk
# sudo apt-get install python-pip

