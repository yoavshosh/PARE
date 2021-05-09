# PARE - Phylogenetic Analysis of RNA Editing

This is a pipeline for the reconstruction of species phylogeny 
and the analysis of the adaptiveness of RNA editing events along the reconstructed phylogeny.

This pipeline contain several modules the should be run sequentialy| with desiered parameters| 
in order to obtain the ancestral sequences and a DB of MSA data (including ancestral sequences) and RNA editing events data.

All scripts were developed in [python3.8](https://www.python.org/downloads/release/python-380/) but should be compatibale also with some python2 enviorments


## Prerequisits

1. Python 3.6 or higher. including [Biopython](https://biopython.org/docs/1.75/api/Bio.html)

2. An MSA software - Currently| [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/) or [Muscle](https://www.ebi.ac.uk/Tools/msa/muscle/) are supported 

3. [pal2nal](http://www.bork.embl.de/pal2nal/) software

4. [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) softwares package


## Data inputs

You should prepare your project directory in a similar manner to the exmp_files directory before executing the pipeline, including the sub-directories and file names
note that some of the files names contain the species name - the substring of the species name is of course changeable and should match your species name

#### 1. Transcriptomes and Editomes:
You can use your own pipeline for the editing detection and transcriptome assambly.
but you will have to adjust the fasta files containing the transcripts and the the tsv files that contain the editing events data for each species
to match the specific headers/fields pattern (see exmp_files).
For the building of this pipeline we have used [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) for the transcriptome assambly and BLASTx(https://blast.ncbi.nlm.nih.gov/Blast.cgi) to infer the ORFS of the assambled Trinity transcripts.
see [This paper](https://www.sciencedirect.com/science/article/pii/S0092867417303446) for more information regarding the differents stages used to obtain the transcriptomes and editomes of the differnet species.

You should also have a version for each of the transcriptome which contains the native (sense-strand) translated region only (see exmp_files).
This is easily generated once you have the transcriptome which also contain data regarding the ORFs and strand direction of each transcript in the file

An example of a transcripts record in the transcriptomes fasta files (Header data format is: "><header_id>	OrfStart	<orf_start>	OrfEnd	<orf_end>	Strand	<stramd_direction>	<protein description>"):

>comp2553489_c0_seq1	OrfStart	3	OrfEnd	287	Strand	-	sp|C9D7C2|CAC1A_APIME Voltage-dependent calcium channel type A subunit alpha-1 OS=Apis mellifera GN=CAC PE=2 SV=1
ATGTCGTTCCTTCCAGACACCGCCACCACGGGGTCCTCGGCCGCCAGAGCCACGGAACTGGCGCAAATCACGATCATGATGAACAGATCGAAGTAGCGGAGGTTGACCACGAAGTGACAGAACTGCCGGATTGGGTTTGTCGGGTAAAAGATGAACATGGAGCTGTAGGGCAACATAGGTTTAGGTCCGGAGAAGCCTCCATCATCGTTCTGGCTGTTGGTGCGTGGTTGTGCGTCCATGGTGCCGCTGCTGGTGGTGCTGGCTCCTAGGACCGCCCCTGGTTTGTG 


An example of the correspondin native coding region:

>apl|comp2553489_c0_seq1
CACAAACCAGGGGCGGTCCTAGGAGCCAGCACCACCAGCAGCGGCACCATGGACGCACAACCACGCACCAACAGCCAGAACGATGATGGAGGCTTCTCCGGACCTAAACCTATGTTGCCCTACAGCTCCATGTTCATCTTTTACCCGACAAACCCAATCCGGCAGTTCTGTCACTTCGTGGTCAACCTCCGCTACTTCGATCTGTTCATCATGATCGTGATTTGCGCCAGTTCCGTGGCTCTGGCGGCCGAGGACCCCGTGGTGGCGGTGTCTGGAAGGAACGAC

(note how here, the name of the species from which the sequence is taken is attached as a pre-fix with "|" to the id)


An example of the editing sites table (headers names are not necessary in the actual files):
|id					|protein|location|mm_type| DNA_A | DNA_T | DNA_G | DNA_C | RNA_A | RNA_T | RNA_G | RNA_C |Trinity| RNA_coverage | DNA_coverage | p_val			   | AA_before | AA_after|type | protein_length| editing_level		  |strand |
|:------------------|:------|-------:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:------------:|:------------:|:-----------------:|:---------:|:---------:|:---:|:-------------:|:--------------------:|:-----:|
|comp1041989_c0_seq1|O13076 |	602  |  AG	 |	7	 |	0	 |	0	 |	0	 |	19	 |	0	 |	3	 |	0	 |	A	 |	22		    |		7	   |7.8125e-06		   |	K		  |	R		|type2|	1200		 |	0.13636363636399998	|	+	|
|comp1045255_c0_seq1|Q80WA4	|	630  |  AG	 |	8	 |	0	 |	0	 |	0	 |	13	 |	0	 |	2	 |	0	 |	A	 |	15		    |		8	   |0.000104094083013  |	K		  |	R		|type2|	1827		 |	0.133333333333		|	+	|
|comp104833_c0_seq1	|Q9EPA7	|	40	 |  AG	 |	0	 |	0	 |	0	 |	0	 |	0	 |	2	 |	0	 |	2	 |	C	 |	4		    |		0	   |5.9920029999999e-06|	T		  |	T		|type2|	213			 |	0.5					|	-	|
|comp113506_c0_seq1	|E9Q3S4	|	57	 |  AG	 |	2	 |	0	 |	0	 |	0	 |	3	 |	0	 |	3	 |	0	 |	A	| 	6	   	    |		2	   |0.00025			   |	E		  |	E		|type2|	204			 |	0.5					|	+	|

Actually, the really mandatory fields here are [id,protein,location,mm_type,editing_level,strand]. If you don't have a table in the exact format above, you should include the fields' headers in the first line of the table file (separated by tabs).

Note that the data in this table (including positions) is relative to the full mRNA sequences (from the transcriptomes), and not relative to the coding sequence which is only a translated region of the transcript, indicated from the strand and the ORFs data in the transcript header.


#### 2. Orthologous gene groups data:
Also you will have to provide data regarding the the orthologous genes in the different transcripts.
You can use [OrthoMCL](https://orthomcl.org/orthomcl/app) to obtain this or use other packages for this task
Anyhow| you should eventually end up with a file containing all the pairs of orthologous genes (transcripts) of each pair of species among the species you have.
see all_best_hits.txt file in the exmp_files folder

An example of the othologous genes table (see also all_best_hits.txt in exmp_files):
|species 1					 |species 2					 |score |
|:--------------------------:|:-------------------------:|:----:|
|apl&#124;comp100003_c0_seq1	|bim&#124;comp128003_c0_seq1	|1.03	|
|apl&#124;comp1000174_c0_seq1	|bim&#124;comp106298_c0_seq1	|0.783	|
|apl&#124;comp10001_c0_seq1		|bim&#124;comp146229_c0_seq1	|0.445	|
|apl&#124;comp100093_c0_seq1		|bim&#124;comp214875_c0_seq1	|0.68	|

The score of each paring here is not needed for the pipline. It is just the output of OrthoMCL used for generating the pairs.
you can fill it with an arbitrary value if you dont have any score for you orthologs pairing 


#### 3. Tree topology
You can have your evolutionary tree topology in a newick format (see tree.txt in exmp_files for an example).
If you need to infer the topology yourself you can use [RAxML-NG](https://github.com/amkozlov/raxml-ng) or other softwares for that.


## The pipeline

On each command, you can use `--help` for more information

#### Multiple MSA procedures of the different orthologous groups. 
You should consider a proteins based MSA as usually, in coding sequences, the protein structure is expected to be conserved along the phylogeny

example command:

`python run_msa_for_super_orthologs.py -parent_dir <project_directory> -o <output_dir_name> -animals <species_for_MSA> -program clu -msa_processes 30 -proteins_msa True`

This command will first create `<output_dir_name>` in which it will stor fasta files contatining the different orthologous groups indicated from all_best_hits.txt.
It will then invoke MSA processes (in this example - 30 in parallel) for each of the orthologous group and will stor the results in `<output_dir_name>/msa_results/`

You can also use `python run_msa_for_super_orthologs.py --help` for more parameters information

You can also use `run_msa_for_super_orthologs_based_proteins_name.py` to group orthologous genes based on the protein name instead of inferring it from the all_best_hits.txt file

#### Converstion of protein MSA to nucleotides MSA
This stage is needed if you invoked the MSAs in protein level in the previous stage

example command:

`python run_pal2nal_for_super_orthologs_to_create_codon_msa.py -native_coding_mrnas <native_coding_mrnas> -o proteins_super_orthologs <output_dir_of_previous_stage>`

This command will create fasta files contatining the mRNA coding sequences of the orthologous gene groups, and will create a subdirectory named `<output_dir_of_previous_stage>\codons_msa_results` in which it will store the converted MSAs in separated folders.


#### Running PAML for Ancestral Sequences Reconstruction (ASR)

Now, when we have all the codon-wise MSA's for all the super orthologs, we shall execute codeML software from the PAML package to perform the ASR

example command:

`python run_paml_codeml_for_super_orthologs.py -msa_dirs_path <path to the codonwise MSAs parent directory> -tree <tree file path> -codeml_processes 50`

This command will invoke 50 codeml processes at once analyzing the MSAs and performing the phylogenetic analysis and ASR

If you are not working on protein coding sequences, you may consider invoking baseml instead using the `run_paml_baseml_for_super_orthologs.py` module


After this process is done, you will have the PAML results in each separate directory of each MSA. you should parse those results to create an MSA with ancestral sequences
as well as some other data files in a more convenient format
use this command:

`python parse_paml_codeml_output.py -msa_dirs_path <path to the codonwise MSAs parent directory> -tree <tree file path>`

#### Building A DB of all MSA columns and editing events data

Now, we should build a DB combining the data from all the MSAs columns + the editing sites data
exmple command:

`python build_full_tree_nucl_matrix.py -parent_dir <your project directory> -super_orthologs_msa_dir <orthologous gene groups sub-dir> -animal <a list of species names> -ancestors <a list of ancestral nodes names> -threads 10`

This command will create one large table combining all the MSAs columns data together with data regarding RNA events in each columns for each species.
It will also create columns containing data regarding the MSA quality (based on agreement with consensus of adjuscent columns) and the ASR quality in each position (based on PAML paresed results)

This is a large table, but the columns names should be self explenatory.


## The Hypotheses Class







