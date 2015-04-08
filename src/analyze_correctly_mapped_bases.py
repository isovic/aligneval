#! /usr/bin/python

import os;
import sys;
import utility_sam;

def CountSameCigarOps(tested_cigar_pos_list, reference_cigar_pos_list):
	num_same_ops = 0;
	num_same_m_ops = 0;
	num_same_i_ops = 0;
	num_same_d_ops = 0;
	num_same_other_ops = 0;
	
	i = 0;
	j = 0;
	
	while (i < len(tested_cigar_pos_list) and j < len(reference_cigar_pos_list)):
		# If the coordinates are not equal, increase one of them and continue.
		if (tested_cigar_pos_list[i][2] < reference_cigar_pos_list[j][2]):
			i += 1;
			continue;
		elif (tested_cigar_pos_list[i][2] > reference_cigar_pos_list[j][2]):
			j += 1;
			continue;
		# If the coordinates are equal, then compare the CIGAR operations.
		else:
			# If the operation and the count are the same, we found a match.
			if (tested_cigar_pos_list[i][0] == reference_cigar_pos_list[j][0] and
			    tested_cigar_pos_list[i][1] == reference_cigar_pos_list[j][1]):
				num_same_ops += 1;
				
				# This counts only match/mismatch operations that are the same.
				if (tested_cigar_pos_list[i][1] in 'M=X'):
					num_same_m_ops += 1;
				# This counts only insertion operations that are the same.
				elif (tested_cigar_pos_list[i][1] == 'I'):
					num_same_i_ops += 1;
				# This counts only deletion operations that are the same.
				elif (tested_cigar_pos_list[i][1] == 'D'):
					num_same_d_ops += 1;
				# This counts all other operations that are the same.
				else:
					num_same_other_ops += 1;
					
				#print tested_cigar_pos_list[i];
			else:
				if (tested_cigar_pos_list[i][1] == 'I' and reference_cigar_pos_list[j][1] != 'I'):
					i += 1;
					continue;
				elif (tested_cigar_pos_list[i][1] != 'I' and reference_cigar_pos_list[j][1] =='I'):
					j += 1;
					continue;
				
				#i = (i + 1) if (tested_cigar_pos_list[i][1] == 'I') else i;
				#j = (j + 1) if (reference_cigar_pos_list[j][1] == 'I') else j;
				#continue;
			
		i += 1;
		j += 1;
		
	return [num_same_ops, num_same_m_ops, num_same_i_ops, num_same_d_ops, num_same_other_ops];

def CountOperations(cigar_list):
	num_m_ops = 0;
	num_i_ops = 0;
	num_d_ops = 0;
	num_other_ops = 0;
	
	for cigar in cigar_list:
		if (cigar[1] in 'M=X'):
			num_m_ops += 1;
		elif (cigar[1] == 'I'):
			num_i_ops += 1;
		elif (cigar[1] == 'D'):
			num_d_ops += 1;
		else:
			num_other_ops += 1;

	num_m_ops = 1 if (num_m_ops == 0) else num_m_ops;
	num_i_ops = 1 if (num_i_ops == 0) else num_i_ops;
	num_d_ops = 1 if (num_d_ops == 0) else num_d_ops;
	num_other_ops = 1 if (num_other_ops == 0) else num_other_ops;
	
	return [num_m_ops, num_i_ops, num_d_ops, num_other_ops];

def CountCorrectlyMappedBasesIndividualSamLine(sam_line, reference_sam_line):
	cigar_pos_list = sam_line.CalcCigarStartingPositions(True);
	counts_sam_line = CountOperations(cigar_pos_list);
	[num_m_ops, num_i_ops, num_d_ops, num_other_ops] = counts_sam_line;

	reference_cigar_pos_list = reference_sam_line.CalcCigarStartingPositions(True);
	counts_reference = CountOperations(reference_cigar_pos_list);
	[ref_num_m_ops, ref_num_i_ops, ref_num_d_ops, ref_num_other_ops] = counts_reference;

	results_same_ops = CountSameCigarOps(cigar_pos_list, reference_cigar_pos_list);
	[num_same_ops, num_same_m_ops, num_same_i_ops, num_same_d_ops, num_same_other_ops] = results_same_ops;

	return [num_same_m_ops, num_m_ops, ref_num_m_ops];

def CountCorrectlyMappedBasesSamLines(sam_line, reference_sam_lines):
	cigar_pos_list = sam_line.CalcCigarStartingPositions(True);
	[num_m_ops, num_i_ops, num_d_ops, num_other_ops] = CountOperations(cigar_pos_list);

	all_counts = [];
	all_counts_dataset = [];

	for sam_reference in reference_sam_lines:
		reference_cigar_pos_list = sam_reference.CalcCigarStartingPositions(True);
		counts_dataset = CountOperations(reference_cigar_pos_list);
		# [ref_num_m_ops, ref_num_i_ops, ref_num_d_ops, ref_num_other_ops] = counts_dataset;

		counts = CountSameCigarOps(cigar_pos_list, reference_cigar_pos_list);
		# [num_same_ops, num_same_m_ops, num_same_i_ops, num_same_d_ops, num_same_other_ops] = counts;

		all_counts.append(counts);
		all_counts_dataset.append(counts_dataset);

	return [all_counts, all_counts_dataset];

def CountMOps(hashed_reference_sam):
	total_ref_num_m_ops = 0;
	i = 0;
	for sam_line_list in hashed_reference_sam.values():
		if ((i % 100) == 0):
			sys.stderr.write('\rLine %d' % (i));
		sam_line = sam_line_list[0];

		num_m_ops = 0;
		num_i_ops = 0;
		num_d_ops = 0;
		num_other_ops = 0;
		cigar_list = sam_line.SplitCigar();
		for cigar in cigar_list:
			if (cigar[1] in 'M=X'):
				num_m_ops += cigar[0];
			elif (cigar[1] == 'I'):
				num_i_ops += cigar[0];
			elif (cigar[1] == 'D'):
				num_d_ops += cigar[0];
			else:
				num_other_ops += cigar[0];

		total_ref_num_m_ops += num_m_ops;
		i += 1;
	sys.stderr.write('\n');
	return total_ref_num_m_ops;

	

def CountCorrectlyMappedBases(sam_file, hashed_reference_sam, out_summary_prefix=''):
	fp_in = None;
	fp_out = None;
	
	out_file = out_summary_prefix + '.csv';
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!\n' % (__name__, sam_file));
		exit(1);
	
	if (out_summary_prefix != ''):
		try:
			fp_out = open(out_file, 'w');
		except IOError:
			sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_file));
			exit(1);

	sys.stderr.write('Starting to count the number of correctly mapped bases in the tested SAM file!\n');

	total_num_same_ops_individual_m = 0;
	total_num_same_m_ops_individual_m = 0;
	total_num_same_i_ops_individual_m = 0;
	total_num_same_d_ops_individual_m = 0;
	total_num_same_other_ops_individual_m = 0;
	dataset_total_ref_num_m_ops = 0;
	dataset_total_ref_num_i_ops = 0;
	dataset_total_ref_num_d_ops = 0;
	dataset_total_ref_num_other_ops = 0;

	i = 0;
	for line in fp_in:
		if ((i % 100) == 0):
			sys.stderr.write('\rLine %d' % (i));

		# if (i > 1000):
		# 	break;
		
		if (len(line.strip()) == 0 or line[0] == '@'):
			i += 1;
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());
		
		if (sam_line.IsMapped() == True):
			sam_reference = None;

			header_key = sam_line.qname;
			# TODO: THIS NEEDS TO BE REMOVED OR IMPLEMENTED SOMEHOW DIFFERENTLY!!
			# The point of this was that, BLASR doesn't conform to the SAM standard, and makes it difficult to
			# uniformly evaluate the results!
			if 'blasr' in sam_file.lower():
				header_key = '/'.join(header_key.split('/')[:-1]);
				if sam_line.clip_count_front != 0 or sam_line.clip_count_back != 0:
					print 'BLASR CIGAR contains clipping! Please revise clipped_pos! Read: "%s".' % sam_line.qname;

			try:
				sam_reference = hashed_reference_sam[header_key][0];
			except:
				sys.stderr.write('\tERROR: Reference SAM does not contain qname "%s"!\n' % (sam_line.qname));
				i += 1;
				continue;
			
			cigar_pos_list_individual_m = sam_line.CalcCigarStartingPositions(True);
			reference_cigar_pos_list_individual_m = sam_reference.CalcCigarStartingPositions(True);
			
			[num_m_ops_individual_m, num_i_ops_individual_m, num_d_ops_individual_m, num_other_ops_individual_m] = CountOperations(cigar_pos_list_individual_m);
			[ref_num_m_ops_individual_m, ref_num_i_ops_individual_m, ref_num_d_ops_individual_m, ref_num_other_ops_individual_m] = CountOperations(reference_cigar_pos_list_individual_m);
			
			dataset_total_ref_num_m_ops += ref_num_m_ops_individual_m;
			dataset_total_ref_num_i_ops += ref_num_i_ops_individual_m;
			dataset_total_ref_num_d_ops += ref_num_d_ops_individual_m;
			dataset_total_ref_num_other_ops += ref_num_other_ops_individual_m;

			[num_same_ops_individual_m, num_same_m_ops_individual_m, num_same_i_ops_individual_m, num_same_d_ops_individual_m, num_same_other_ops_individual_m] = CountSameCigarOps(cigar_pos_list_individual_m, reference_cigar_pos_list_individual_m);
			total_num_same_ops_individual_m += num_same_ops_individual_m;
			total_num_same_m_ops_individual_m += num_same_m_ops_individual_m;
			total_num_same_i_ops_individual_m += num_same_i_ops_individual_m;
			total_num_same_d_ops_individual_m += num_same_d_ops_individual_m;
			total_num_same_other_ops_individual_m += num_same_other_ops_individual_m;

		i += 1;

	perc_num_same_m_ops_individual_m = (float(total_num_same_m_ops_individual_m) / float(dataset_total_ref_num_m_ops)) * 100.0;
	perc_num_same_i_ops_individual_m = (float(total_num_same_i_ops_individual_m) / float(dataset_total_ref_num_i_ops)) * 100.0;
	perc_num_same_d_ops_individual_m = (float(total_num_same_d_ops_individual_m) / float(dataset_total_ref_num_d_ops)) * 100.0;
	
	fp_in.close();

	if (out_summary_prefix != ''):
		fp_out.write('percent_correct_m\tnum_correct_m\tnum_m_ops_in_reference\n');
		fp_out.write('%.2f\t%.2f\t%.2f\n' % (perc_num_same_m_ops_individual_m, total_num_same_m_ops_individual_m, dataset_total_ref_num_m_ops));
		fp_out.close();
	
	sys.stderr.write('\n');

	return [perc_num_same_m_ops_individual_m, total_num_same_m_ops_individual_m, dataset_total_ref_num_m_ops];
	


def CompareCigars(sam_file, sam_reference_file, out_summary_prefix=''):
	fp_in = None;
	fp_out = None;
	
	out_file = out_summary_prefix + '.csv';
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
 		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!\n' % (__name__, sam_file));
		exit(1);
	
	if (out_summary_prefix != ''):
		try:
			fp_out = open(out_file, 'w');
		except IOError:
			sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_file));
			exit(1);
	
	sys.stderr.write('Loading reference SAM file...\n');
	qnames_with_multiple_alignments = {};
	[hashed_reference, num_references, num_unique_references] = utility_sam.HashSAMWithFilter(sam_reference_file, qnames_with_multiple_alignments);

	sys.stderr.write('Starting to process the tested SAM file!\n');

	all_values = [];

	# total_num_same_ops = 0;
	# total_num_same_m_ops = 0;
	# total_num_same_i_ops = 0;
	# total_num_same_d_ops = 0;
	# total_num_same_other_ops = 0;
	total_num_same_ops_individual_m = 0;
	total_num_same_m_ops_individual_m = 0;
	total_num_same_i_ops_individual_m = 0;
	total_num_same_d_ops_individual_m = 0;
	total_num_same_other_ops_individual_m = 0;
	dataset_total_ref_num_m_ops = 0;
	dataset_total_ref_num_i_ops = 0;
	dataset_total_ref_num_d_ops = 0;
	dataset_total_ref_num_other_ops = 0;

	i = 0;
	for line in fp_in:
		sys.stderr.write('\rLine %d' % (i));
		
		if (len(line.strip()) == 0 or line[0] == '@'):
			i += 1;
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());
		
		if (sam_line.IsMapped() == True):
			sam_reference = None;
			try:
				sam_reference = hashed_reference[sam_line.qname][0];
			except:
				sys.stderr.write('\tERROR: Reference SAM does not contain qname "%s"!\n' % (sam_line.qname));
				i += 1;
				continue;
			
			cigar_pos_list = sam_line.CalcCigarStartingPositions();
			reference_cigar_pos_list = sam_reference.CalcCigarStartingPositions();
			cigar_pos_list_individual_m = sam_line.CalcCigarStartingPositions(True);
			reference_cigar_pos_list_individual_m = sam_reference.CalcCigarStartingPositions(True);
			
			[num_m_ops, num_i_ops, num_d_ops, num_other_ops] = CountOperations(cigar_pos_list);
			[ref_num_m_ops, ref_num_i_ops, ref_num_d_ops, ref_num_other_ops] = CountOperations(reference_cigar_pos_list);
			
			[num_m_ops_individual_m, num_i_ops_individual_m, num_d_ops_individual_m, num_other_ops_individual_m] = CountOperations(cigar_pos_list_individual_m);
			[ref_num_m_ops_individual_m, ref_num_i_ops_individual_m, ref_num_d_ops_individual_m, ref_num_other_ops_individual_m] = CountOperations(reference_cigar_pos_list_individual_m);
			
			
			
			dataset_total_ref_num_m_ops += ref_num_m_ops;
			dataset_total_ref_num_i_ops += ref_num_i_ops;
			dataset_total_ref_num_d_ops += ref_num_d_ops;
			dataset_total_ref_num_other_ops += ref_num_other_ops;
			
			# [num_same_ops, num_same_m_ops, num_same_i_ops, num_same_d_ops, num_same_other_ops] = CountSameCigarOps(cigar_pos_list, reference_cigar_pos_list);
			# total_num_same_ops += num_same_ops;
			# total_num_same_m_ops += num_same_m_ops;
			# total_num_same_i_ops += num_same_i_ops;
			# total_num_same_d_ops += num_same_d_ops;
			# total_num_same_other_ops += num_same_other_ops;

			#[num_same_ops_individual_m, num_same_m_ops_individual_m] =
			[num_same_ops_individual_m, num_same_m_ops_individual_m, num_same_i_ops_individual_m, num_same_d_ops_individual_m, num_same_other_ops_individual_m] = CountSameCigarOps(cigar_pos_list_individual_m, reference_cigar_pos_list_individual_m);
			total_num_same_ops_individual_m += num_same_ops_individual_m;
			total_num_same_m_ops_individual_m += num_same_m_ops_individual_m;
			total_num_same_i_ops_individual_m += num_same_i_ops_individual_m;
			total_num_same_d_ops_individual_m += num_same_d_ops_individual_m;
			total_num_same_other_ops_individual_m += num_same_other_ops_individual_m;
			#print 'Number of same CIGAR operations: %d' % (num_same_ops);
			#print 'Number of CIGAR operations in tested SAM: %d' % (len(cigar_pos_list));
			#print 'Number of CIGAR operations in reference SAM: %d' % (len(reference_cigar_pos_list));
			#print ' ';

		i += 1;

	# perc_num_same_m_ops = (float(total_num_same_m_ops) / float(dataset_total_ref_num_m_ops)) * 100.0;
	# perc_num_same_i_ops = (float(total_num_same_i_ops) / float(dataset_total_ref_num_i_ops)) * 100.0;
	# perc_num_same_d_ops = (float(total_num_same_d_ops) / float(dataset_total_ref_num_d_ops)) * 100.0;

	perc_num_same_m_ops_individual_m = (float(total_num_same_m_ops_individual_m) / float(dataset_total_ref_num_m_ops)) * 100.0;
	perc_num_same_i_ops_individual_m = (float(total_num_same_i_ops_individual_m) / float(dataset_total_ref_num_i_ops)) * 100.0;
	perc_num_same_d_ops_individual_m = (float(total_num_same_d_ops_individual_m) / float(dataset_total_ref_num_d_ops)) * 100.0;
	
	fp_in.close();

	if (out_summary_prefix != ''):
		fp_out.write('percent_correct_m\tpercent_correct_i\tpercent_correct_d\n');
		fp_out.write('%.2f\t%.2f\t%.2f\n' % (perc_num_same_m_ops_individual_m, perc_num_same_i_ops_individual_m, perc_num_same_d_ops_individual_m));
		fp_out.close();
	
	sys.stderr.write('\n');

	return perc_num_same_m_ops_individual_m;
	
	

# src/sam_compare_cigars.py /home/ivan/work/eclipse-workspace/golden-bundle/alignments_for_testing/reads-simulated/PacBio-100k/escherichia_coli/graphmap-params_SSW_r-test-drive.sam /home/ivan/work/eclipse-workspace/golden-bundle/reads-simulated/PacBio-100k/escherichia_coli/reads.sam temp/cigcompare.temp
if __name__ == "__main__":
	if (len(sys.argv) < 3 or len(sys.argv) > 4):
		sys.stderr.write('Compares CIGAR string operations from one SAM line to the corresponding line in the reference SAM (the SAM file generated by a read simulator).\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <input_sam_file> <input_sam_reference> [<out_summary_prefix>]\n' % sys.argv[0]);
		exit(1);
	
	sam_file = sys.argv[1];
	sam_reference_file = sys.argv[2];
	out_summary_prefix = '';
	if ((sys.argv) == 4):
		out_summary_prefix = sys.argv[3];
	
	print 'Percent correctly mapped bases: %.2f' % CompareCigars(sam_file, sam_reference_file, out_summary_prefix);
