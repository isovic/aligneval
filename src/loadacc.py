#! /usr/bin/python

# Data structure for holding the info loaded from file. This class is a containter for only one line of
# the input data.
class AccuracyInfo:
	def __init__(self, correct=False, mapq=255, header_def=False, ref_header='*', map_ref_header='*', ref_pos=-1, map_pos=-1, distance=-1, qname='', secondary=False, ref_reverse=False, reverse=False, mapped=False, paired=False, max_allowed_distance=50):
		self.correct = correct;
		self.qname = qname;			# Sequence name, as found in the aligner's output SAM file.
		self.secondary = secondary;	# SAM alignments have a secondary flag, which can be used to note alternative alignments.
		self.ref_reverse = ref_reverse;
		self.reverse = reverse;
		self.mapped = mapped;
		self.paired = paired;
		self.header_def = header_def;
		self.ref_header = ref_header;
		self.map_ref_header = map_ref_header;
		self.ref_pos = ref_pos;
		self.map_pos = map_pos;
		self.distance = distance;
		self.max_allowed_distance = max_allowed_distance;
		self.mapq = mapq;
	
	def FormatLine(self):
		line = 'qname = %s\tmapped = %s\tcorrect = %s\tmapq = %d\theader_def = %s\tref_header = %s\tmap_header = %s\tref_pos = %d\tmap_pos = %d\tdistance = %d\tsecondary = %s\tref_reverse = %s\treverse = %s\tpaired = %s\tallowed_dist = %d' % (self.qname, (self.mapped), str(True if self.correct != 0 else False), self.mapq, str(self.header_def), str(self.ref_header), str(self.map_ref_header), self.ref_pos, self.map_pos, self.distance, str(self.secondary), str(self.ref_reverse), str(self.reverse), str(self.paired), self.max_allowed_distance);
		return line;
	
	def ParseLine(self, line):
		split_line = line.split('\t');
		
		for param in split_line:
			split_param = param.split('=');
			if len(split_param) < 2:
				continue;
			
			param_name = split_param[0].strip();
			param_value = split_param[1].strip();
			
			if param_name != 'qname':
				eval_line = 'self.%s' % param;
			else:
				eval_line = 'self.%s = \'%s\'' % (param_name, param_value);
			
			exec(eval_line);

# Loads and parses all accuracy info from a given path, and returns it as an array of AccuracyInfo objects.
def LoadAccuracies(acc_path):
	try:
		fp = open(acc_path, 'r');
		lines = fp.readlines();
		fp.close();
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % acc_path;
		return [];

	accuracies = [];
	
	for line in lines:
		acc_info = AccuracyInfo();
		acc_info.ParseLine(line);
		accuracies.append(acc_info);
		
	return accuracies;

# Dummy function that needs to do something.
def AnalyzeAccuracies(accuracies):
	print 'Some cleaver analyses are being performed...';
	
	# This line sorts all accuracies by their qname. Thought it might be useful.
	sorted_accuracies = sorted(accuracies, key=lambda x: x.qname, reverse=True)
	num_all_alignments = len(sorted_accuracies);
	
	i = 0;
	for accuracy in sorted_accuracies:
		# Do something cleaver.
		print 'qname[%d] = %s' % (i, accuracy.qname);
		i += 1;
	
	print 'Done already!';
	
	return [];

if __name__ == "__main__":
	path = '../results/reads-simulated/PacBio/escherichia_coli/bbmap.acc'
	accuracies = LoadAccuracies(path);
	print 'In total, %d entries loaded.' % (len(accuracies));
	AnalyzeAccuracies(accuracies);
