#! /usr/bin/python

machine_names = [];
simulated_datasets = [];

simulated_datasets.append('Illumina-1k-single_end');									machine_names.append('illumina');
simulated_datasets.append('PacBio-1k');												machine_names.append('pacbio');
# simulated_datasets.append('OxfordNanopore-pbsim-observed_last-2d-1k');				machine_names.append('nanopore');
# simulated_datasets.append('OxfordNanopore-pbsim-observed_graphmap-2d-1k');			machine_names.append('nanopore');
# simulated_datasets.append('OxfordNanopore-pbsim-observed_last-1d-1k');				machine_names.append('nanopore');
# simulated_datasets.append('OxfordNanopore-pbsim-observed_graphmap-1d-1k');			machine_names.append('nanopore');
# simulated_datasets.append('OxfordNanopore-pbsim-observed_marginalign-2d-1k');			machine_names.append('nanopore');

# simulated_datasets.append('Illumina-0k-single_end');								machine_names.append('illumina');
# simulated_datasets.append('OxfordNanopore-pbsim-observed_marginalign-2d-0k');		machine_names.append('nanopore');

genomes = [];
# genomes.append('neisseria_meningitidis');
# genomes.append('escherichia_coli');
genomes.append('saccharomyces_cerevisiae');
# genomes.append('caenorhabditis_elegans');
# genomes.append('hg19_v38-chr3');
