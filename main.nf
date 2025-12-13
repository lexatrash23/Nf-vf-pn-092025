#!/usr/bin/env nextflow

// Params output for user 

def printhead() {
    log.info("")
    log.info("════════════════════════════════════════════════════════════")
    log.info("_    _________   ______  __  ___________    ____ _       __")
    log.info("| |  / / ____/ | / / __ \\/  |/  / ____/ /   / __ \\ |     / /")
    log.info("| | / / __/ /  |/ / / / / /|_/ / /_  / /   / / / / | /| / /")
    log.info("| |/ / /___/ /|  / /_/ / /  / / __/ / /___/ /_/ /| |/ |/ /")
    log.info("|___/_____/_/ |_/\\____/_/  /_/_/   /_____/\\____/ |__/|__/")
    log.info("")
    log.info("════════════════════════════════════════════════════════════")
    log.info("")
    log.info("Author:                          ${workflow.manifest.author}")
    log.info("README:                          ${workflow.manifest.homePage}")
    log.info("Description:                     ${workflow.manifest.description}")
    log.info("Version:                         ${workflow.manifest.version}")
    log.info("Start:                           ${workflow.start}")
    log.info("")
    log.info("────────────────────────────────────────────────────────────")
    log.info("")
}

// Process 1: PostTrimFastqc
process PostTrimFastqc {
    errorStrategy 'ignore'

    label 'process_low'

    conda "fastqc=0.12.1"
    container "docker://biocontainers/fastqc:v0.11.9_cv8"

    input:
    tuple val(sample), path(R1), path(R2)

    output:
    tuple val(sample), path("*.zip"), emit: fastqc_zips

    script:
    """
    fastqc ${R1} ${R2}
    """
}

// Process 2: MultiQC
process MultiQC {
    errorStrategy 'ignore'

    label 'process_single'

    conda "multiqc=1.33"
    container "docker://multiqc/multiqc:v1.32"

    publishDir "${sample}/Venomflow/results/Fastqc/posttrim/", mode: 'copy'

    input:
    tuple val(sample), path(fastqc_zips)

    output:
    path 'multiqc_report.html', emit: multiqc

    script:

    """
    multiqc ${fastqc_zips}
    """
}

// Process 3: Bowtie build and index 
process Bowtie {

    label 'process_high'
    label 'process_long'

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "bowtie2=2.5.4 samtools=1.22.1"

    publishDir "${sample}/Venomflow/results/Bowtie/", pattern: "*.bam", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), path(R1), path(R2)

    output:
    path "*.bam"

    script:

    """
    bowtie2-build ${trinity_fasta} ${sample}_transcriptome_index
    bowtie2 -x ${sample}_transcriptome_index -1 ${R1} -2 ${R2} -S ${sample}_mapped_reads.sam -p 8 --no-unal
    samtools view -bS ${sample}_mapped_reads.sam | samtools sort -o ${sample}_mapped_reads.bam



    """
}
// Process 4: TrinityStats
process TrinityStats {

    errorStrategy 'ignore'

    label 'process_single'

    conda "seqkit=2.12.0"

    publishDir "${sample}/Venomflow/results/Stats/", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta)

    output:
    path "*.txt", emit: trinity_stats

    script:

    """
    seqkit stats ${trinity_fasta} > ${sample}_Trinity.stats.txt

    """
}
// Process 5: BUSCO_transcriptome_metazoa
process BUSCO_transcriptome_metazoa {

    errorStrategy 'ignore'

    label 'process_low'

    conda "busco=5.8.3"

    publishDir "${sample}/Venomflow/results/BUSCO/transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), val(metazoa)

    output:
    path "*.txt", emit: busco_transcriptome_met

    script:

    """
    busco -i ${trinity_fasta} -l ${metazoa} -c 10 -o ${sample}_met.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_met.transcriptome/*.txt "."
    """
}

// Process 6: BUSCO_transcriptome_mollusca
process BUSCO_transcriptome_mollusca {

    errorStrategy 'ignore'

    label 'process_low'

    conda "busco=5.8.3"

    publishDir "${sample}/Venomflow/results/BUSCO/transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), val(mollusca)

    output:
    path "*.txt", emit: busco_transcriptome_mollusca

    script:

    """
    busco -i ${trinity_fasta} -l ${mollusca} -c 10 -o ${sample}_mol.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_mol.transcriptome/*.txt "."
    """
}



// Process 7: Kallisto_Trinity
process Kallisto_Trinity {

    conda "kallisto=0.51.1"

    label 'process_low'
    label 'process_long'

    publishDir "${sample}/Venomflow/results/kallisto/trinity/output", mode: 'copy'

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    input:
    tuple val(sample), path(trinity_fasta), path(R1), path(R2), val(Strandedness)

    output:
    path "abundance.tsv", emit: KallistoTrinityAbundance

    script:
    """
    stranded_input="${Strandedness}"

    if [ ${Strandedness} == "fr" ]; then
        stranded_flag="--fr-stranded"
    elif [ ${Strandedness} == "rf" ]; then
        stranded_flag="--rf-stranded"
    else
        stranded_flag=""
    fi
    kallisto index -i index ${trinity_fasta}
    kallisto quant -i index -o ./ -b 100 ${R1} ${R2} \$stranded_flag

    """
}

// Process 8: Blastdatabasecreation
process Blastdatabasecreation {
    errorStrategy 'ignore'

    label 'process_single'

    conda "blast=2.17.0"

    input:
    tuple val(sample), path(database_fasta)

    output:
    tuple val(sample), path("*"), emit: proteindb

    script:
    """
    makeblastdb -in "${database_fasta}" -dbtype prot
    """
}

// Process 9: Blastx
process Blastx {

    label 'process_medium'
    label 'process_long'

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "blast=2.17.0"

    publishDir "${sample}/Venomflow/results/Blast/Blastx/", mode: 'copy'

    input:
    tuple val(sample), val(proteindbdbname), path(trinity_fasta), path(proteindb)

    output:
    path "${sample}.blastx.db.0.txt", emit: blastx0
    path "${sample}.blastx.db.6.txt", emit: blastx6

    script:
    """
   
    blastx -query ${trinity_fasta} -db ${proteindbdbname} -out ${sample}.blastx.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"
    blastx -query ${trinity_fasta} -db ${proteindbdbname} -out ${sample}.blastx.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """
}

// Process 10: Transdecoder
process Transdecoder {

    label 'process_single'

    conda "transdecoder=5.7.1"

    publishDir "${sample}/Venomflow/results/Transdecoder", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta)

    output:
    tuple val(sample), path("*.pep"), emit: pep
    tuple val(sample), path("*.cds"), emit: cds

    script:

    """
    TransDecoder.LongOrfs -t ${trinity_fasta}
    TransDecoder.Predict -t ${trinity_fasta}

    """
}

// Process 11: BUSCO_translatome_metazoa
process BUSCO_translatome_metazoa {

    label 'process_medium'

    errorStrategy 'ignore'

    conda "busco=5.8.3"

    publishDir "${sample}/Venomflow/results/BUSCO/translatome/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), val(metazoa)

    output:
    path "*.txt", emit: busco_translatome_met

    script:

    """
    busco -i ${Transdecoder_pep} -l ${metazoa} -c 10 -o ${sample}_met.protein -m protein -e 1e-5 -f
    mv ${sample}_met.protein/*.txt "."
    
    """
}

// Process 12: BUSCO_translatome_mollusca
process BUSCO_translatome_mollusca {

    label 'process_medium'

    errorStrategy 'ignore'

    conda "busco=5.8.3"

    publishDir "${sample}/Venomflow/results/BUSCO/translatome/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), val(mollusca)

    output:
    path "*.txt", emit: busco_translatome_met

    script:

    """
    busco -i ${Transdecoder_pep} -l ${mollusca} -c 10 -o ${sample}_mol.protein -m protein -e 1e-5 -f
    mv ${sample}_mol.protein/*.txt "."
    
    """
}


// Process 13: Kallisto_Transdecoder
process Kallisto_Transdecoder {

    errorStrategy 'ignore'

    label 'process_medium'

    conda "kallisto=0.51.1"

    publishDir "${sample}/Venomflow/results/kallisto/transdecoder/output", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_cds), path(R1), path(R2), val(Strandedness)

    output:
    path "abundance.tsv", emit: KallistoTransdecoderAbundance

    script:

    """

    if [ ${Strandedness} == "fr" ]; then
        stranded_flag="--fr-stranded"
    elif [ ${Strandedness} == "rf" ]; then
        stranded_flag="--rf-stranded"
    else
        stranded_flag=""
    fi
    kallisto index -i index ${Transdecoder_cds}
    kallisto quant -i index -o ./ -b 100 ${R1} ${R2} \$stranded_flag

    """
}

// Process 14: Blastp
process Blastp {

    errorStrategy 'ignore'

    label 'process_low'

    conda "blast=2.17.0"

    publishDir "${sample}/Venomflow/results/Blast/Blastp/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), path(proteindb), val(proteindbdbname)

    output:
    path "${sample}.blastp.db.0.txt", emit: blastp0
    path "${sample}.blastp.db.6.txt", emit: blastp6

    script:
    """

    blastp -query ${Transdecoder_pep} -db ${proteindbdbname} -out ${sample}.blastp.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"
    blastp -query ${Transdecoder_pep} -db ${proteindbdbname} -out ${sample}.blastp.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """
}

// Process 15: Transdecoder_complete
process Transdecoder_complete {

    errorStrategy 'ignore'

    label 'process_single'

    conda "seqkit=2.12.0"

    publishDir "${sample}/Venomflow/results/Transdecoder", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), path(Transdecoder_cds)

    output:
    tuple val(sample), path("*.pep"), emit: transdecodercomplete_pep
    tuple val(sample), path("*.cds"), emit: transdecodercomplete_cds

    script:

    """
    seqkit grep -n -r -p "ORF type:complete" ${Transdecoder_pep} -o "${sample}.trandescoder.complete.pep"
    seqkit grep -n -r -p "ORF type:complete" ${Transdecoder_cds} -o "${sample}.trandescoder.complete.cds"
    """
}

// Process 16: SignalP
process SignalP {

    errorStrategy 'ignore'

    label 'process_low'

    publishDir "${sample}/Venomflow/results/Signalp", mode: 'copy'

    input:
    tuple val(sample), path(transdecodercomplete_pep)

    output:
    tuple val(sample), path("*_mature.fasta"), emit: maturesequences
    path "*.signalp5", emit: signalpsummary

    script:


    """
    signalp -fasta ${transdecodercomplete_pep} -mature -prefix "${sample}"
    """
}

// Process 17: Filter2
process Filter2 {
    errorStrategy 'ignore'

    label 'process_single'

    conda "seqkit=2.12.0"

    publishDir "${sample}/Venomflow/results/Transdecoder", mode: 'copy'

    input:
    tuple val(sample), path(maturesequences), path(transdecodercomplete_pep)

    output:
    tuple val(sample), path("*.fasta"), emit: transdecoderpep_signalp

    script:

    """
    seqkit grep -f <(seqkit seq -n ${maturesequences}) ${transdecodercomplete_pep} > ${sample}.Transdecoder.complete.signalp.sequences.fasta
    """
}

// Process 18: STATS
process stats {
    errorStrategy 'ignore'

    label 'process_single'

    conda "seqkit=2.12.0"

    publishDir "${sample}/Venomflow/results/Stats", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), path(Transdecoder_cds), path(transdecodercomplete_pep), path(transdecodercomplete_cds), path(maturesequences), path(transdecoderpep_signalp)

    output:
    path "*"

    script:

    """
    seqkit stats ${Transdecoder_pep} > ${sample}_Transdecoder_pep.stats.txt
    seqkit stats ${Transdecoder_cds} > ${sample}_Transdecoder_cds.stats.txt
    seqkit stats ${transdecodercomplete_pep} > ${sample}_complete_pep.stats.txt
    seqkit stats ${transdecodercomplete_cds} > ${sample}_transdecodercomplete_cds.stats.txt
    seqkit stats ${maturesequences} > ${sample}_maturesequences.stats.txt
    seqkit stats ${transdecoderpep_signalp} > ${sample}_transdecoderpep_signalp.stats.txt
    
    """
}

// Process 19: Interproscan
process Interproscan {

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_high'
    label 'process_long'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/Venomflow/results/Interproscan", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep)

    output:
    path "*.tsv", emit: Interproscan

    script:

    """
    awk '{if (\$0 ~ /^>/) print \$0; else {gsub(/\\*/, ""); print \$0}}' ${Transdecoder_pep} > "${sample}.Trinity.fasta.transdecoder.cleaned.pep"

    interproscan.sh -goterms -i "${sample}.Trinity.fasta.transdecoder.cleaned.pep" -pa -t p -d ./ -f TSV 
    """
}

// Process 20: GenomeBlastdatabasecreation
process GenomeBlastdatabasecreation {

    errorStrategy 'ignore'

    conda "blast=2.17.0"

    label 'process_single'

    input:
    tuple val(sample), path(genome_fasta)

    output:
    tuple val(sample), path("*"), emit: genomedbfiles

    script:
    """
    makeblastdb -in "${genome_fasta}" -dbtype nucl
    """
}

// Process 21: GenomeBlasts
process GenomeBlasts {

    label 'process_high'
    label 'process_long'

    conda "blast=2.17.0"

    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    publishDir "${sample}/Venomflow/results/Blast/Blastn/", mode: 'copy'

    input:
    tuple val(sample), val(genomedbname), path(transdecodercompletecds), path(genomedb)

    output:
    path "${sample}.blastn.db.0.txt", emit: blastn0
    path "${sample}.blastn.db.6.txt", emit: blastn6

    script:
    """
    blastn -query ${transdecodercompletecds} -db ${genomedbname} -out ${sample}.blastn.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"
    blastn -query ${transdecodercompletecds} -db ${genomedbname} -out ${sample}.blastn.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """
}







// PARAMETERS DEFINED IN CONFIG AND SAMPLESHEET FILE

//WorkFlow

workflow {

    // Print head 
    printhead()

    // Define CSV channel. Chanel factory creates a channel from a csv file. This csv can be defined in the config or the command line with --input_csv. The CSV is spilt row by row where each row is emitted as a Map (key-value pairs) where column names are news. 
    csv_channel = Channel.fromPath(params.input_csv).splitCsv(header: true, sep: ',')

    // Define Input: Paired Trimmed reads Tuple. Extracts as tuple the sample name and the trimmed reads
    input_R1R2 = csv_channel.map { row -> tuple(row.Sample_name, file(row.R1), file(row.R2)) }

    // Run Process: Fastqc
    input_R1R2 | PostTrimFastqc

    //Define Input: Fastqc htmls. Using the output of Fastqc as an input for MultiQC
    fastqc_output = PostTrimFastqc.out.fastqc_zips
    //Run Process: MultQC
    fastqc_output | MultiQC

    // Define Input: Bowtie 
    input_bowtie = csv_channel.map { row -> tuple(row.Sample_name, file(row.Trinity_fasta), file(row.R1), file(row.R2)) }
    //Run Process: Bowtie
    input_bowtie | Bowtie

    //Define Input: Trinity Fasta 
    input_trinity_fasta = csv_channel.map { row -> tuple(row.Sample_name, file(row.Trinity_fasta)) }
    //Run Process: TrinityStats
    input_trinity_fasta | TrinityStats

    //Define Input: Trinity Fasta + BUSCOlin1 tuple 
    input_BUSCOlin1 = csv_channel.map { row -> tuple(row.Sample_name, file(row.Trinity_fasta), file(row.BUSCO_lin1)) }

    //Run Process: BUSCO_lin1
    input_BUSCOlin1 | BUSCO_transcriptome_metazoa

    //Define Input: Trinity Fasta + BUSCOlin2 tuple 
    input_BUSCOlin2 = csv_channel.map { row -> tuple(row.Sample_name, file(row.Trinity_fasta), file(row.BUSCO_lin2)) }

    //Run Process: BUSCO_lin2
    input_BUSCOlin2 | BUSCO_transcriptome_mollusca

    //Define Input: Trinity fasta + R1 + R2 tuple 
    input_TrinityKallisto = csv_channel.map { row -> tuple(row.Sample_name, file(row.Trinity_fasta), file(row.R1), file(row.R2), row.Strandedness) }

    //Run Process: TrinityKallisto
    input_TrinityKallisto | Kallisto_Trinity

    //Define Input: Database_Fasta. This is set up to allow for multiple different databases if needed.
    input_databasefasta = csv_channel.map { row -> tuple(row.Sample_name, file(row.Protein_fasta_path_for_Blast)) }


    //Run Process: DatabaseCreation
    input_databasefasta | Blastdatabasecreation

    //Define Input: Blastx 
    Blastxinputfasta = csv_channel.map { row ->
        tuple(row.Sample_name, row.Protein_fasta_name, file(row.Trinity_fasta))
    }
    Blastxinputwithdb = Blastxinputfasta.join(Blastdatabasecreation.out.proteindb)


    //Run Process: Blastx
    Blastxinputwithdb | Blastx

    //Run Process: Transdecoder
    input_trinity_fasta | Transdecoder

    //Define Input: Transdecoder pep + BUSCOlin1 tuple 
    Transdecoder_pep = Transdecoder.out.pep
    BUSCOlin1 = csv_channel.map { row -> tuple(row.Sample_name, row.BUSCO_lin1) }
    input_BUSCOlin1_L = Transdecoder_pep.join(BUSCOlin1)

    //Run Process: BUSCO_lin1
    input_BUSCOlin1_L | BUSCO_translatome_metazoa

    //Define Input: Transdecoder pep + BUSCOlin2 tuple 
    BUSCOlin2 = csv_channel.map { row -> tuple(row.Sample_name, row.BUSCO_lin2) }
    input_BUSCOlin2_L = Transdecoder_pep.join(BUSCOlin2)

    //Run Process: BUSCO_lin2
    input_BUSCOlin2_L | BUSCO_translatome_mollusca

    //Define Input: Transdecoder cds + R1 + R2 + Strandedness tuple 
    Transdecoder_cds = Transdecoder.out.cds
    KallistoTransdecoderR1R2S = csv_channel.map { row -> tuple(row.Sample_name, file(row.R1), file(row.R2), row.Strandedness) }
    input_TransKallisto = Transdecoder_cds.join(KallistoTransdecoderR1R2S)

    //Run Process: TransKallisto
    input_TransKallisto | Kallisto_Transdecoder

    //Define Input: Blastp - Match Transdecoder output with databases
    Blastpdb = csv_channel.map { row ->
        tuple(row.Sample_name, row.Protein_fasta_name)
    }
    input_Blastp = Transdecoder_pep.join(Blastdatabasecreation.out.proteindb).join(Blastpdb)

    //Run Process: Blastp
    input_Blastp | Blastp

    //Define Input: Sample name + PEP + CDS tuple
    input_Transdecoder_complete = Transdecoder_pep.join(Transdecoder_cds)

    //Run Process: Transdecoder filter for complete ORFs
    input_Transdecoder_complete | Transdecoder_complete


    //Define Input: Sample name + completepep tuple
    transdecodercomplete_pep = Transdecoder_complete.out.transdecodercomplete_pep

    //Run Process: SignalP
    transdecodercomplete_pep | SignalP

    //Define Input: sample name + mature sequences + completepep tuple 
    maturesequences = SignalP.out.maturesequences
    maturecomplete = maturesequences.join(transdecodercomplete_pep)

    //Run Process: Filter2   
    maturecomplete | Filter2

    //Define Input: stats sample name + pep + cds + completepep + completecds + mature +signalp  tuple 
    stats_join = input_Transdecoder_complete
        .join(Transdecoder_complete.out.transdecodercomplete_pep)
        .join(Transdecoder_complete.out.transdecodercomplete_cds)
        .join(maturesequences)
        .join(Filter2.out.transdecoderpep_signalp)
    //Run Process: STATS   
    stats_join | stats

    //Run Process: Interproscan  (Input was previously defined)
    Transdecoder_pep | Interproscan

    //Define Input: genomefasta (with path tracking)
    Genomefasta = csv_channel
        .map { row -> tuple(row.Genome_fasta_path, file(row.Genome_fasta_path)) }
        .unique { it[0] }
        .filter { genome_path, _genome_file ->
            genome_path != 'NULL' && genome_path.toString() != 'NULL'
        }

    //Run Process: BlastnGenome database creation
    Genomefasta | GenomeBlastdatabasecreation

    //Define Input: Genome BLAST 
    genomedb = GenomeBlastdatabasecreation.out.genomedbfiles
    Blastncds = csv_channel.map { row -> tuple(row.Sample_name, file(row.Genome_fasta_name)) }

    BlastnInput = Blastncds
        .join(Transdecoder_complete.out.transdecodercomplete_cds)
        .join(genomedb)

    //Run Process: BlastnGenome
    BlastnInput | GenomeBlasts
}
