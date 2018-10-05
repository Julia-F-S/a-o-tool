# a-o-tool
Exploiting orthology and de novo transcriptome assembly to refine target sequence information

**1) Getting started**
  - get Nextflow (https://www.nextflow.io/)
  - get docker (https://docs.docker.com/)
  - clone this repository
  - run the a&o-tool with
    ```
    nextflow run seq_refinement_pipeline.nf -with-docker jfsoellner/a-o-tool_dependencies --configFile <your config file>
    ```
**2) Example data**  
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
    
 Once the pipeline has finished there will be a directory in your current working directory called "Q15436_kidney_pig" which contains all the output. The 
