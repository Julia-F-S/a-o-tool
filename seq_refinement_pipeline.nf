#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * The aim is to determine the sequence of a poorly annotated query gene via ortholog sequences.
 * 
 * Required input:
 * - --configFile (see git readme for details) 
 * 
 * This pipeline executes the following steps:
 * 1) Assembly based on RNASeq reads using BinPacker (alternatively an existing assembly can be used)
 * 2) Sequence Ortholog Contig Fishing via BLAST (if human is the ortholog, swissprot sequences are used for reciprocal best BLAST)
 * 3) Longest ORFs extraction
 * 4) Protein translation
 * 5) Clustalw\/MUSCLE multiple sequence alignment
 *
 */
 
/*
 * command line input
 */
parameter_json = file(params.configFile)

// validate command line inputs
if( !parameter_json.exists() ) exit 1, "Config file could not be found: params.configFile"

// parse config file for parameters 
new groovy.json.JsonSlurper().parseText(parameter_json.text).each { k, v -> params[k] = v }
refined_species = params.refined_species
ortholog_species = params.ortholog_species
prot_sequences = file(params.prot_sequences)
blastDBname = params.blastDB // name to be used for BLAST database when creating it
pubDir = file("$params.publishDir")
println "$pubDir"
email = params.email

// set default directory to store BLAST databases in
blastDBdir = file("$baseDir/blastdb/assemblies") 
  
//println(msa_input_map)

// check if msa_n_contigs has been provided, if not set default of 1 (= number of contigs to display in the MSA)
if ( params.msa_n_contigs!="" ) {
  msa_n_contigs = params.msa_n_contigs
}
else {
  msa_n_contigs = 1
}

// make publishDir if it doesn't exist
if( !pubDir.exists() ) {
  new File("$pubDir").mkdir()  
}

// use existing assembly if one is provided, otherwise use the provided fastq file for assembly with BinPacker
if ( params.assembly!="" ) {
  assembly = file(params.assembly)
}
else {
  left_fastq = file(params.left_fastq)
  right_fastq = file(params.right_fastq)
}

// check if human Swiss-Prot sequences were provided
if( params.human_swissprot!="" ) {
  human_swissprot = file(params.human_swissprot)
}

/* 
 * Check if an assembly has been provided in the config file. If not: Use input fastq files for assembly with BinPacker (+ check if fastq files exist).
 */ 
if( params.assembly=="" ) {
  println "assembly doesn't exist, starting BinPacker"
  // run BinPacker on the input fastq files
  process binPackerAssembly {
     publishDir "$pubDir", mode: 'copy'

     input:
      file left_fastq
      file right_fastq
 
     output:
      file '**/BinPacker.fa' into assembly

     """ 

     if [ ! -e $left_fastq ] || [ ! -e $right_fastq ]
     then 
       echo "Fastq file(s) containing reads for the assembly could not be found." 
       exit 1 
     else
       BinPacker -s fq -p pair -l $left_fastq -r $right_fastq > tmp
     fi

     """
    }
}
else {
  // use the provided assembly
  println "assembly found, skipping BinPacker run"
}

/*
 * retrieve protein sequence of the ortholog bait
 */ 
Channel
     .fromPath(prot_sequences)
     .splitFasta( record: [id: true, seqString: true, text: true])
     .filter { record ->  record.id =~ ortholog_species  } 
     .map { record -> record.text }
     .set { ortholog_seq_bait }
     
/*
 * 1a) Make BLAST database from assembly
 */
process makeBlastDatabase_assembly {
  publishDir "$blastDBdir", mode: 'copy'


  input:
  file assembly
  val blastDBname

  output:
  file "${blastDBname}.*" into blastDb

  script:
  """  

  makeblastdb -dbtype nucl -in ${assembly} -out ${blastDBname} -parse_seqids -hash_index 

  """
}

/*
 * 1b) Make BLAST database from human Swiss-Prot sequences, if provided
 */
if (params.human_swissprot != "") {
  process makeBlastDatabase_hsSwissprot {
    publishDir "$blastDBdir", mode: 'copy'

    input:
    file human_swissprot

    output:
    file "human_swissprot.*" into blastDb_hs

    script:
    """

    makeblastdb -dbtype prot -in ${human_swissprot} -out human_swissprot -parse_seqids -hash_index 

    """
  }
}

/* 
 * 2) Find best matching contig via ortholog bait: 
 *    tblastn (query: ortholog bait sequence; db: assembly)
 */
process tblastn {
    publishDir "$pubDir", mode: 'copy'
    
    input:
    file blastDb
    file ortholog_seq_bait

    output:
    file "tblastn_output.txt" into blast_result

    """

    tblastn -query ${ortholog_seq_bait} \
    -db $blastDBname \
    -out "tblastn_output.txt" \
    -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' \
    -evalue 0.0001 \
    -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 \
    -max_target_seqs $msa_n_contigs -num_threads 10 

    """
}

blast_result.into{blast_result_getContigs; blast_result_checkRBH}

/* 
 * Get best matching contigs (names)
 */
process getContigNames {
 publishDir "$pubDir", mode: 'copy'
 
 input:
  file blast_result_getContigs
  val msa_n_contigs
 
 output:
  file "contigs.txt" into contigs
   
 """

 head -n $msa_n_contigs ${blast_result_getContigs} | cut  --fields 2 | sort | uniq > contigs.txt

 """
}

/* 
 * Get best matching contigs (sequences)
 */
 process getContigFasta {
  publishDir "$pubDir", mode: 'copy'
 
  input:
   file 'contigs.txt' from contigs
   file blastDb
 
  output:
   file "contigs.fa" into contigFasta
 
  """ 

  blastdbcmd -db $blastDBname -entry_batch contigs.txt > contigs.fa
  
  """
}

contigFasta.into{contigFasta_blastx; contigFasta_filter}

/*
 * If human Swiss-Prot sequences have been provided: align the best hit contigs to Swiss-prot sequences and determine whether they lead to a reciprocal best hit
 * Assumption: the FASTA header of the human sequence contains the UniProt/Swiss-Prot accession number
 */ 

process blastx {
  publishDir "$pubDir", mode: 'copy'

  input:
  file contigFasta_blastx
  file blastDb_hs

  output:
  file "blastx_output.txt" into blastX_result
  when: 
  params.human_swissprot != ""
 
  """

  blastx -query ${contigFasta_blastx} \
  -db human_swissprot \
  -out "blastx_output.txt" \
  -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' \
  -evalue 0.0001 \
  -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 \
  -num_alignments 3 -num_threads 10

  """
}
  
process checkRBH {
  publishDir "$pubDir", mode: 'copy'

  input:
  file blastX_result
  file blast_result_checkRBH

  output:
  file "contigs_reciprocalHit.txt" into contigs_RH
  file "contigs_reciprocalBestHit.txt" into contigs_RBH
  when: 
  params.human_swissprot != ""
 
  """

  #!/opt/conda/envs/ao-tool_v1/bin/Rscript
      
  # get best hit (human Swiss-Prot protein) for the query contig
  blastx <- read.table("${blastX_result}", sep="\t", header=F, stringsAsFactors = F)
  tblastn <- read.table("${blast_result_checkRBH}", sep="\t", header=F, stringsAsFactors = F)
    
  colnames(blastx) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe", "qlen", "slen")
  colnames(tblastn) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe", "qlen", "slen")
  
  currProt <- strsplit(unique(tblastn\$qaccver), "_")[[1]][2] # fasta header is formated as "human_proteinAccno"
  contigs <- unique(tblastn\$saccver)
  
  rbhRes <- data.frame(contig=contigs, hasRBH=rep(NA_integer_, length(contigs)), anyHit=rep(NA_integer_, length(contigs)))
  for (it in seq_along(contigs)) {
    rbhRes[it, "contig"] <- contigs[it]
    
    currX <- blastx[blastx\$qaccver==contigs[it], "saccver"]
    rbhRes[it, "hasRBH"] <- currX[1] == currProt
    rbhRes[it, "anyHit"] <- currProt %in% currX
  }
  
  rbhContigs <- rbhRes[which(rbhRes\$hasRBH==1), "contig"]
  rhContigs <- rbhRes[which(rbhRes\$anyHit==1), "contig"]
  
  write.table(rbhContigs, file="contigs_reciprocalBestHit.txt", col.names=F, row.names=F, quote=F)
  write.table(rhContigs, file="contigs_reciprocalHit.txt", col.names=F, row.names=F, quote=F)

  """
}

/*
 * if human swiss prot sequences were available: filter contigs for those with a reciprocal best BLAST hit
 */ 
if( params.human_swissprot != "" ){

  process filterContigsRBH {
    publishDir "$pubDir", mode: 'copy'

    input:
    file contigs_RBH
    file contigs_RH
    file blastDb

    output:
    file "contigs_reciprocal*.fa" into filteredContigs
    when: 
    params.human_swissprot != ""
 
    """

    blastdbcmd -db $blastDBname -entry_batch "${contigs_RBH}" > contigs_reciprocalBestHit.fa
    blastdbcmd -db $blastDBname -entry_batch "${contigs_RH}" > contigs_reciprocalHit.fa

    """
  }
} else {
  contigFasta_filter.into{filteredContigs}
}

// transform the list of files emitted by the filteredContigs channel into a channel that emits each file object independently
filteredContigs
    .flatMap()
    .set{filteredContigs_1}

/* 
 * find open reading frames in contigs (TransDecoder)
 */
process getORF {
  publishDir "$pubDir", mode: 'copy'
  
  input:
    file x from filteredContigs_1

  output:
    file "*.transdecoder.pep" into longestOrfs

  """

  if [ -s "${x}" ] 
  then
    echo 'passedOuter'
    # find ORFs
    TransDecoder.LongOrfs -t "${x}"  
  
    # make sure ORFs were found before proceeding to prediction step 
    if [ -s "${x}.transdecoder_dir/longest_orfs.pep" ]
    then
      echo 'passedInner'
      # determine coding potential & get best ORF
      TransDecoder.Predict  -t "${x}" --no_refine_starts --single_best_only

      # enhance displayed contig name
      sed -i 's/>/>${refined_species}_refined_/g' *.transdecoder.pep
    else
      echo 'notPassedInner'
      touch "empty.${x}.transdecoder.pep"
    fi
  else
    echo 'notPassedOuter'
    touch "empty.${x}.transdecoder.pep"
  fi

  """ 
}


/* 
 * multiple sequence alignment 
 */
process muscle {
 
  publishDir "$pubDir", mode: 'copy'
  
  input:
    file prot_sequences
    each x from longestOrfs.collect()

  
  output:
    file '*inputMSA_inclORF*.fa' into inputMSAorf
    file '*msa*.aln' into msaOutClustal
    file '*msa*.html' into msaOutHtml
    
  """

  if [ -s "$pubDir/${x.baseName}.pep" ]
  then
    cat "$pubDir/${x.baseName}.pep" > "inputMSA_inclORF_${x.baseName}.fa"
    cat ${prot_sequences} >> "inputMSA_inclORF_${x.baseName}.fa"
  
    muscle -in "inputMSA_inclORF_${x.baseName}.fa" -htmlout "msa_${x.baseName}.html" -clwout "msa_${x.baseName}.aln" 
  else
    touch "empty.inputMSA_inclORF_${x.baseName}.fa"
    touch "empty.msa_${x.baseName}.aln"
    touch "empty.msa_${x.baseName}.html"
  fi

  """ 
}

// transform the list of files emitted by the filteredContigs channel into a channel that emits each file object independently
inputMSAorf
    .flatMap()
    .set{inputMSAorf_1}

/*
 * compute % identity of original and refined sequence compared to bait species and vice versa.
 */
process percIdentity {
  publishDir "$pubDir", mode: 'copy'
  
  input:
    val ortholog_species
    file x from inputMSAorf_1
    
  output:
    file 'percentIdentityStats*.txt' 
    
  """

  #!/opt/conda/envs/ao-tool_v1/bin/python
  
import re
from Bio import SeqIO
from Bio import pairwise2 as pw2
from Bio.SubsMat.MatrixInfo import blosum62
import os

if os.stat("${x}").st_size != 0:
  # read the sequences  
  record_dict = SeqIO.to_dict(SeqIO.parse("${x}", "fasta"))

  # get the species which has been refined
  r = re.compile(".*refined")
  refined = filter(r.match, record_dict.keys())
  species = refined[0].split("_")[0]

  # get original, refined  + human sequences
  seqHs = [value for key, value in record_dict.items() if 'human' in key.lower()]
  seqOrig = [value for key, value in record_dict.items() if (species in key.lower() and 'refined' not in key.lower())]
  seqRef = [value for key, value in record_dict.items() if (species in key.lower() and 'refined' in key.lower())]

  # open file & write header
  f = open("percentIdentityStats_${x}.txt", 'w') 
  f.write("name1 \\t name2 \\t n_matches \\t query_perc_identity \\t target_perc_identity \\t target_minus_query_perc_identity \\t score \\t length_aa_ortholog \\t length_aa_hs \\n")

  # add original values (if original sequence(s) available)
  if len(seqOrig) != 0:
    for oIt in seqOrig:
        for hsIt in seqHs:
            # compute pairwise alignments to calculate % identities + delta
            origGlobal_align = pw2.align.globalds(hsIt, oIt, blosum62, -10, -0.5)
          
            # compute percentage as (matches_in_alignment / len) * 100
            matches = sum(a==b for a, b in zip(origGlobal_align[0][0], origGlobal_align[0][1]))
            pMatchOrig_hs = round(((matches*100) / len(hsIt)), 2) # query %ID
            pMatchOrig_orig = round(((matches*100) / len(oIt)), 2) # target %ID
            pMatchOrig_delta = pMatchOrig_orig - pMatchOrig_hs
            scoreOrig = ((1-(pMatchOrig_hs/100)) * (1-(pMatchOrig_orig/100))) + (abs(pMatchOrig_delta)/100)

            f.write(oIt.name + "\\t" + hsIt.name + "\\t" + str(matches) + "\\t" + str(pMatchOrig_hs) + "\\t" + str(pMatchOrig_orig) + "\\t" + str(pMatchOrig_delta) + "\\t" + str(scoreOrig) + "\\t" + str(len(oIt))  + "\\t" + str(len(hsIt))+ "\\n")

  # iterate over the refined and human sequences to compute pairwise stats for all of them & write them to file
  for refIt in seqRef:
    for hsIt in seqHs:
        seqRefined = str(refIt.seq)
        seqRefined = seqRefined.replace("*", "") # remove termination symbol

        # compute pairwise alignments to calculate % identities + delta
        refinedGlobal_align = pw2.align.globalds(hsIt, seqRefined, blosum62, -10, -0.5)

        # compute percentage as (matches_in_alignment / len) * 100
        matches = sum(a==b for a, b in zip(refinedGlobal_align[0][0], refinedGlobal_align[0][1]))
        pMatchRefined_hs = round(((matches*100) / len(hsIt)), 2) # query %ID
        pMatchRefined_refined = round(((matches*100) / len(seqRefined)), 2) # target %ID
        pMatchRefined_delta = pMatchRefined_refined - pMatchRefined_hs
        scoreRefined = ((1-(pMatchRefined_hs/100)) * (1-(pMatchRefined_refined/100))) + (abs(pMatchRefined_delta)/100)

        f.write(refIt.name + "\\t" + hsIt.name + "\\t" + str(matches) + "\\t" + str(pMatchRefined_hs) + "\\t" + str(pMatchRefined_refined) + "\\t" + str(pMatchRefined_delta) + "\\t" + str(scoreRefined) + "\\t" + str(len(seqRefined))  + "\\t" + str(len(hsIt))+ "\\n")

  f.close()
else:
  f = open("percentIdentityStats_${x}.txt", 'w')
  f.close()
  
  """
}


if ( email!="" ) {
  workflow.onComplete {
      def subject = 'My pipeline execution'
      def recipient = email
  
      ['mail', '-s', subject, recipient].execute() << """

      Pipeline execution summary
      ---------------------------
      Cmd line    : ${workflow.commandLine}
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      publishDir  : $pubDir
      exit status : ${workflow.exitStatus}
      Error report: ${workflow.errorReport ?: '-'}
      """
  }
}


