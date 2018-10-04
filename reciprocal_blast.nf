#!/usr/bin/env nextflow

/*
 * Pipeline for the general validation approach:
 * 1) make BLAST database from assembly
 * 2) tblastn (query: (human Swiss-Prot) protein sequences; db: assembly)
 * 3) get best hits per protein
 * 4) get sequence of best hits
 * 5) blastx (query: best hits; db: human Swiss-Prot protein sequences)
 *
 * Input:
 * - path to the assembly (fasta)
 * - path to the human Swiss-Prot protein sequences (fasta)
 *
 *
 * Based on an example by the nextflow author(s).
 */
 

/*
 * The pipeline inputs parameters which have to be specified as command line options
 */
assembly = file(params.assembly)
params.dbDir = "$baseDir/blastdb"
params.dbName = assembly.baseName
params.swissprot = "$baseDir/human_swissprot.fasta"
params.swissprotDB = "swissprot_human"
params.email = ""

dbName = params.dbName
dbDir = params.dbDir
swissprot = file(params.swissprot)
swissprotDB = params.swissprotDB
tissue = params.tissue
species = params.species
pubDir = file("./" + params.pubDir)

println assembly 

// deduce assembler from assembly name
if (assembly.baseName == "BinPacker") {
    assembler = "binpacker"
} else if (assembly.baseName == "transcripts") {
    assembler = "spades"
} else {
    println "Could not guess assembler: " + assembly.baseName
    exit 1
}

// make publishDir if it doesn't exist
if( !pubDir.exists() ) {
  new File("$pubDir").mkdir()  
}


/*
 * 1) Make BLAST database from assembly
 */
 process makedb {
    publishDir "$dbDir", mode: 'copy'

    input:
    file assembly

    output:
    file db

    """
    makeblastdb -in $assembly -dbtype nucl -out $dbDir/$dbName -parse_seqids -hash_index > db
    """
}


/* 
 * 2) tblastn (query: human Swiss-Prot protein sequences; db: assembly)
 */
process tblastn {
    publishDir "$pubDir", mode: 'copy'
    
    input:
    val db
    file swissprot
    val dbDir
    val dbName
    val assembler

    output:
    file "blast_swissprot2${assembler}_${tissue}.txt" into blast_result
  
    """
    tblastn -query $swissprot \
    -db $dbDir/$dbName \
    -out "blast_swissprot2${assembler}_${tissue}.txt" \
    -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' \
    -evalue 0.0001 \
    -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 \
    -max_target_seqs 3 -num_threads 10 

    """
}

/*
 * 3) get best hits
 */
process getBestHit {
    publishDir "$pubDir", mode: 'copy'

    input:
    val tissue
    val assembler
    file blast_result

 
    output:
    file "bestHits_blast_swissprot2${assembler}_${tissue}.txt" into top_hits
 
    """
    #!/usr/bin/Rscript
      
    # get best hit (contig) for each query (human Swiss-Prot protein)
    currRes <- read.table("${blast_result}", sep="\t", header=F, stringsAsFactors = F)

    colnames(currRes) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart",
                             "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe", "qlen", "slen")
      
    queries <- unique(currRes[,"qaccver"])

    fName <- "bestHits_blast_swissprot2${assembler}_${tissue}.txt"
    file.create(fName)
     
    for (it in seq_along(queries)) {
      currD <- currRes[which(currRes[,"qaccver"]==queries[it]), ]
        
      # sort by evalue and bitscore and choose top one
      currD <- currD[with(currD, order(-bitscore, evalue)),]
        
      write(currD[1, "saccver"], file=fName, append=T)
        
    }
    """
}
 
/*
 * 4) get sequences of best hits
 */
process extract {
    publishDir "$pubDir", mode: 'copy'

    input:
    file top_hits
    val dbDir
    val dbName
 
    output:
    file sequences
 
    """
    blastdbcmd -db $dbDir/$dbName -entry_batch $top_hits > sequences
    """
}


/*
 * 5) blastx (query: best hits; db: human Swiss-Prot protein sequences)
 */
process blastx {
    publishDir "$pubDir", mode: 'copy'

    input:
    file sequences
    val dbDir

    output:
    file "blast_${assembler}2swissprot_${tissue}.txt"
 
    """
    blastx -query $sequences \
    -db $dbDir/$swissprotDB \
    -out "blast_${assembler}2swissprot_${tissue}.txt" \
    -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' \
    -evalue 0.0001 \
    -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 \
    -num_alignments 3 -num_threads 10
    """
}

  workflow.onComplete {
      def subject = 'Reciprocal BLAST execution'
      def recipient = params.email
  
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