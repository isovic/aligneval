#! /usr/bin/python

import sys;
import utility_sam;

def parse_vcf_positions(vcf_file):
	try:
		fp = open(vcf_file, 'r');
		lines = fp.readlines();
		fp.close();
	except Exception, e:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n');
		sys.stderr.write(str(e) + '\n');
		exit(1);

	positions = [];
	for line in lines:
		if (line[0] == '#'):
			continue;
		split_line = line.strip().split('\t');
		positions.append(int(split_line[1]) - 1);

	return positions;



def verbose_usage_and_exit():
	sys.stderr.write('Takes a VCF file of simulated reference mutations, extracts positions, and checks alignments of simulated reads. For every position, number of correctly placed bases is counted.\n');
	sys.stderr.write('\n');
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s <alignments.sam> <reference.sam> <mutations.vcf>' % sys.argv[0]);
	sys.stderr.write('\n');

	exit(0);

def main():
	if (len(sys.argv) != 4):
		verbose_usage_and_exit();

	query_sam = sys.argv[1];
	reference_sam = sys.argv[2];
	vcf_file = sys.argv[3];

	[hashed_query, num_queries, num_unique_queries] = utility_sam.HashSAMWithFilter(query_sam, {});
	[hashed_reference, num_references, num_unique_references] = utility_sam.HashSAMWithFilter(reference_sam, {});
	positions = parse_vcf_positions(vcf_file);

	out_summary_prefix = os.path.basename(vcf_file);
	accuracy = utility_sam.CountCorrectlyMappedBasesAtPositions(hashed_query, hashed_reference, positions, out_summary_prefix=out_summary_prefix);
	sys.stderr.write('Accuracy: %.2f\n' % accuracy);



if __name__ == "__main__":
	main();