# a-o-tool 
Exploiting orthology and de novo transcriptome assembly to refine target sequence information

The a&o-tool is a Nextflow pipeline which uses a transcriptome assembly of a species of interest and a known orthologous protein sequence from a closely related species to determine the best matching assembled contigs in the species of interest. These contigs are then searched for open reading frames and translated into an amino acid sequence. A multiple sequence alignment of the translated sequence, the orthologous protein and optional sequences from other species is calculated to visually assess the quality of the output. Furthermore, some metrics on the results are provided.

A paper describing the details of the pipeline will be submitted soon.

# Getting started
  - get Nextflow (https://www.nextflow.io/)
  - get docker (https://docs.docker.com/)
  - clone this repository
  - run the a&o-tool with
    ```
    nextflow run seq_refinement_pipeline.nf -with-docker jfsoellner/a-o-tool_dependencies --configFile <your config file>
    ```
  - the Docker image jfsoellner/a-o-tool_dependencies contains all required software and will be pulled automatically
    
# Input
The pipeline requires a config file in json format. To start the pipeline use:

The fields of the config file:
- *assembly*: path to the pre-computed de novo transcriptome assembly in fasta format. If the assembly shall be computed, set to "" and provide paired-end reads via the next two fields.
- *left_fastq*: If no assembly is provided, the R1 reads have to be provided in fastq format. They are used as input for the assembler. If an assembly is provided set to "".
- *right_fastq*: If no assembly is provided, the R2 reads have to be provided in fastq format. They are used as input for the assembler. If an assembly is provided set to "".
- *refined_species*: Name of the species of interest. 
- *ortholog_species*: Name of the closely related species with known protein sequence.
- *human_swissprot*: If you use human as the ortholog_species you can provide a path to a fasta containing all human sequences available in UniprotKB/Swiss-Prot. This file will then be used for a reciprocal BLAST approach. If you do not use human or you do not want to perform this step, set to "".
- *prot_sequences*: Path to the fasta file containing the known orthologous protein sequence as well as optional sequences from related species or, if available, the currently annotated sequence in the species of interest. **The fasta header must be formated as species_xx, e.g. pig_ENSSSCT00000047607.**
- *msa_n_contigs*: Number of matching contigs from the assembly which should be included in the multiple sequence alignment. If left blank, i.e. "", only the first contig will be shown.
- *publishDir*: Name or path to directory where results should be stored.
- *blastDB*: Name for the BLAST database created for the assembled transcriptome.
- *email*: Upon completion the pipeline will report to the e-mail address provided in the config file. If you do not want to receive an email, set it to "".

In the following example config file we provide the assembly as a fasta file because it has been computed prior to the pipeline run.
If the pipeline should trigger the assembly, simply set the assembly field's value to "" and provide the fastq files to be used for assembly. For both cases an example is included in this repository. 

```
{
"assembly": "./demo_data/BinPacker_assembly_pig_kidney.fa",
"left_fastq": "",
"right_fastq": "",
"refined_species": "pig",
"ortholog_species": "human",
"human_swissprot": "./demo_data/human_swissprot.fasta",
"prot_sequences": "./demo_data/Q15436_orthologs.fa",
"msa_n_contigs": "5",
"publishDir": "./Q15436_kidney_pig",
"blastDB": "blastdb_pig_binpacker_kidney",
"email": ""
}
```

# Output
The main outputs of the pipeline are:
- the multiple sequence alignment stored in the publishDir provided in the config file. It is available in CLUSTALW format (msa.aln) and as coloured html (msa.html)
- a textfile reporting the sequence identities and the difference in sequence identity between the original sequence (if available) and the human protein as well as between the refined sequence(s) and the human protein. 

If you performed a reciprocal BLAST search there will be two files each: One for the reciprocal best hit and one for all hits. If no reciprocal hit was found but the human protein was maybe the second best hit, the contig will be shown as hit. In this case, the files for reciprocal best hits will still be created but their name starts with "empty.".

# Example data
We have included example data in `demo_data` which you can use to familiarise yourself with the output of our tool. In this example we know that the human protein transport protein Sec23A is well conserved across several mammals and that the pig sequence available in UniprotKB/Swiss-Prot only covers a short part of the human sequence. Thus, we wonder whether the protein in pig is really shorter or whether we can find evidence for a longer and better matching sequence.    
  
  The data include:   
  - a configuration file for the human protein transport protein Sec23A (121_config_Q15436_pig_kidney.txt)     
  - a fasta file with some orthologs which we want to be displayed in the multiple sequence alignment (Q15436_orthologs.fa)     
  - all human protein sequences from UniprotKB/Swiss-Prot (human_swissprot.fasta)     
  - an assembly from pig (BinPacker_assembly_pig_kidney.fa.gz, needs to be unziped before you can run the pipeline)     

 To run the pipeline use the command  
 ```
 nextflow run seq_refinement_pipeline.nf -with-docker jfsoellner/a-o-tool_dependencies --configFile demo_data/121_config_Q15436_pig_kidney.txt
 ```
    
 Once the pipeline has finished there will be a directory in your current working directory called "Q15436_kidney_pig" which contains the following output: 

![img/example_out.png](https://github.com/Julia-F-S/a-o-tool/blob/master/img/example_out.png)

