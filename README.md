## Aligneval - Framework for evaluating third and second generation sequencing alignments.

### Description
Benchmarking sequence alignment tools is not a trivial task, especially taken into account the variety of mappers, sequencing technologies, simulators and the definition of objective measures of quality of alignment.  

This project aims at developing a unified automatic framework for benchmarking genomic sequence mapping/alignment tools.  
Each mapper/aligner is wrapped up with a Python script which has standard interfaces, including functions that automatically install and configure the corresponding tool.  
Datasets are specified in such a way to cover a wide spectrum of genomic sizes, repeats, number of chromosomes, and sequencing platforms. Once the datasets are generated and the tools installed, the entire alignment and benchmarking process can be run automatically in a single batch.

At the moment the benchmark process only works on simulated data, whereas benchmarking on real data will be added in future releases (soon).

Please note: in order to measure the time/memory consumption, Cgroups API is used. To be able to measure these statistics, the user running the benchmark needs to be allowed the privileges. For this reason, sudo is required.

### Usage
Running the ```setup.py``` script lists the possible setup options.  
To install everything (including the generation of simulated datasets), run:  
```./setup.py all```  

If only tools need to be installed, one can run:  
```./setup.py tools```  

This command will look up all wrapper scripts in the ```wrappers``` folder, and will run the ```download_and_install()``` function from each of those.  

Alternatively, one can install each tool independantly:  
```wrappers/wrapper_graphmap.py install```  

Once a tool is installed, its wrapper can be used as a standalone script for running the tool, e.g.:  
```wrappers/wrapper_graphmap.py run reads.fastq reference.fa nanopore out_folder/```  

If the simulated data has been generated, the alignment process of all mappers can be instantiated with:  
```./run-alignment.py```

Once all alignments have finished, they can be evaluated with:  
``` ./run-evaluation.py ```  
  
  
### Results
After ```run-evaluation.py``` script finishes executing, all results will be collected in a CSV format, and stored in the ```results/``` folder, in files with appropriate names.

Individual alignments for each dataset can be found in the ```evaluation/``` folder.
