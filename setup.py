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

def download_aligners():
	sys.stderr.write('Installing alignment tools.\n');

	aligner_wrappers = basicdefines.find_files(basicdefines.WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py');

	for wrapper in aligner_wrappers:
		wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0];
		command = 'import %s; %s.download_and_install()' % (wrapper_basename, wrapper_basename);
		exec(command);

def setup_tools():
	sys.stderr.write('Cloning Cgmemtime Git repo. Git needs to be installed.\n');
	command = 'cd %s; git clone https://github.com/isovic/cgmemtime.git' % (basicdefines.TOOLS_ROOT_ABS);
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

# tools/last-534/scripts




if __name__ == '__main__':
	sys.stderr.write('Running this script will consume large amount of disk space.\n');
	yes_no = raw_input("Do you want to continue? [y/n] ");

	if (yes_no != 'y'):
		sys.stderr.write('Exiting.\n\n');
		exit(0);

	create_folders();

	# unpack_reference_genomes();

	download_aligners();

	# setup_tools();



# sys.stderr.write('Generating simulated data...\n');
# generate_data.GenerateAll();
# sys.stderr.write('Finished generating simulated data!\n');



# sudo apt-get install python-matplotlib
# sudo apt-get install python-tk
# sudo apt-get install python-pip

