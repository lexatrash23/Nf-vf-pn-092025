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


    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

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
    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_single'
    container "docker://multiqc/multiqc:v1.32"

    conda "multiqc=1.33"
    container "docker://multiqc/multiqc:v1.32"

    publishDir "${sample}/Venomflow/results/Fastqc/posttrim/", mode: 'copy'
    publishDir "${sample}/Analysis/results/htmls/", mode: 'copy'

    input:
    tuple val(sample), path(fastqc_zips)

    output:
    path 'multiqc_report.html', emit: multiqc

    script:

    """
    multiqc ${fastqc_zips}
    """
}

// Process 3: Bowtie build and index # Transcriptome 1 
process Bowtie {

    label 'process_high'
    label 'process_long'

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "bowtie2=2.5.4"
    container "docker://quay.io/biocontainers/bowtie2:2.5.4--he20e202_1"

    publishDir "${sample}/Venomflow/results/Bowtie/Transcriptome1", pattern: "*.log", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), path(R1), path(R2)

    output:
    path "*.log"

    script:

    """
    bowtie2-build ${trinity_fasta} ${sample}_transcriptome_index
    bowtie2 -x ${sample}_transcriptome_index -1 ${R1} -2 ${R2} -S ${sample}_mapped_reads.sam -p 8 --no-unal 2> stats.log



    """
}
// Process 3: Bowtie build and index # Transcriptome 1 

process Bowtie2 {

    label 'process_high'
    label 'process_long'

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "bowtie2=2.5.4"
    container "docker://quay.io/biocontainers/bowtie2:2.5.4--he20e202_1"

    publishDir "${sample}/Venomflow/results/Bowtie/Transcriptome2", pattern: "*.log", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta2), path(R1), path(R2)

    output:
    path "*.log"

    script:

    """
    bowtie2-build ${trinity_fasta2} ${sample}_transcriptome_index
    bowtie2 -x ${sample}_transcriptome_index -1 ${R1} -2 ${R2} -S ${sample}_mapped_reads.sam -p 8 --no-unal 2> stats.log



    """
}


// Process 4: TrinityStats 1 
process TrinityStats {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_single'

    conda "seqkit=2.12.0"
    container "docker://quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0"

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

// Process 4: TrinityStats 2
process TrinityStats2 {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_single'

    conda "seqkit=2.12.0"
    container "docker://quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0"

    publishDir "${sample}/Venomflow/results/Stats/", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta2)

    output:
    path "*.txt", emit: trinity_stats

    script:

    """
    seqkit stats ${trinity_fasta2} > ${sample}_Trinity2.stats.txt

    """
}
// Process 5: BUSCO_transcriptome_metazoa
process BUSCO_transcriptome_metazoa {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_high'
    label 'process_long'

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"


    publishDir "${sample}/Venomflow/results/BUSCO/transcriptome/Transcriptome1", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), val(metazoa), val(transcriptome1_label)

    output:
    path "*.txt"

    script:

    """
    busco -i ${trinity_fasta} -l ${metazoa} -c 10 -o ${sample}_${transcriptome1_label}_met.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_${transcriptome1_label}_met.transcriptome/*.txt "."
    """
}

// Process 6: BUSCO_transcriptome_mollusca
process BUSCO_transcriptome_mollusca {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_high'
    label 'process_long'

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${sample}/Venomflow/results/BUSCO/transcriptome/Transcriptome1", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), val(mollusca), val(transcriptome1_label)

    output:
    path "*.txt"

    script:

    """
    busco -i ${trinity_fasta} -l ${mollusca} -c 10 -o ${sample}_${transcriptome1_label}_mol.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_${transcriptome1_label}_mol.transcriptome/*.txt "."
    """
}

// Process 5: BUSCO_transcriptome_metazoa
process BUSCO_transcriptome_metazoa2 {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_high'
    label 'process_long'

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"


    publishDir "${sample}/Venomflow/results/BUSCO/transcriptome/Transcriptome2", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta2), val(metazoa), val(transcriptome2_label)

    output:
    path "*.txt"

    script:

    """
    busco -i ${trinity_fasta2} -l ${metazoa} -c 10 -o ${sample}_${transcriptome2_label}_met.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_${transcriptome2_label}_met.transcriptome/*.txt "."
    """
}

// Process 6: BUSCO_transcriptome_mollusca
process BUSCO_transcriptome_mollusca2 {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_high'
    label 'process_long'

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${sample}/Venomflow/results/BUSCO/transcriptome/Transcriptome2", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta2), val(mollusca), val(transcriptome2_label)

    output:
    path "*.txt"

    script:

    """
    busco -i ${trinity_fasta2} -l ${mollusca} -c 10 -o ${sample}_${transcriptome2_label}_mol.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_${transcriptome2_label}_mol.transcriptome/*.txt "."
    """
}

// Process 6: Label transcriptomes, combine and remove duplicates 
process Transcriptome_Combined {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_low'

    conda "seqkit=2.12.0"
    container "docker://gfanz/seqkit"

    publishDir "${sample}/Venomflow/results/Transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(transcriptome1), val(transcriptome1_label), path(transcriptome2), val(transcriptome2_label)

    output:
    tuple val(sample), path("*combined.deduplicated.fasta"), emit: transcriptome_combined

    script:

    """
    seqkit replace -p '(.+)' -r '${transcriptome1_label}_\$1' ${transcriptome1} > Transcriptome1_labelled.fasta
    seqkit replace -p '(.+)' -r '${transcriptome2_label}_\$1' ${transcriptome2} > Transcriptome2_labelled.fasta
    cat Transcriptome1_labelled.fasta Transcriptome2_labelled.fasta > Transcriptome_combined.fasta
    seqkit rmdup Transcriptome_combined.fasta -s -o ${sample}_transcriptome_combined.deduplicated.fasta

    """
}


// Process 5: BUSCO_transcriptome_metazoa Combined transcriptome
process BUSCO_transcriptome_metazoa3 {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_high'
    label 'process_long'

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"


    publishDir "${sample}/Venomflow/results/BUSCO/transcriptome/Combined", mode: 'copy'

    input:
    tuple val(sample),val(metazoa),path(combined_trinity)

    output:
    path "*.txt"

    script:

    """
    busco -i ${combined_trinity} -l ${metazoa} -c 10 -o ${sample}_combined.met.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_combined.met.transcriptome/*.txt "."
    """
}

// Process 6: BUSCO_transcriptome_mollusca
process BUSCO_transcriptome_mollusca3 {

    errorStrategy 'retry'
    maxRetries 4

    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_high'
    label 'process_long'

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${sample}/Venomflow/results/BUSCO/transcriptome/Combined", mode: 'copy'

    input:
    tuple val(sample), val(mollusca),path(combined_trinity)

    output:
    path "*.txt"

    script:

    """
    busco -i ${combined_trinity} -l ${mollusca} -c 10 -o ${sample}_combined.mol.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_combined.mol.transcriptome/*.txt "."
    """
}

// Process 7: Kallisto_Trinity
process Kallisto_Trinity {

    conda "kallisto=0.51.1"
    container "docker://quay.io/biocontainers/kallisto:0.51.1--h2b92561_2"


    label 'process_low'
    label 'process_long'

    publishDir "${sample}/Venomflow/results/kallisto/trinity/output", mode: 'copy'

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    input:
    tuple val(sample), path(R1), path(R2), val(Strandedness), path(combined_trinity)

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
    kallisto index -i index ${combined_trinity}
    kallisto quant -i index -o ./ -b 100 ${R1} ${R2} \$stranded_flag

    """
}

// Process 8: Blastdatabasecreation
process Blastdatabasecreation {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_single'

    conda "blast=2.17.0"
    container "docker://ncbi/blast:2.17.0"

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

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "blast=2.17.0"
    container "docker://ncbi/blast:2.17.0"

    publishDir "${sample}/Venomflow/results/Blast/Blastx/", mode: 'copy'

    input:
    tuple val(sample), val(proteindbdbname), path(combined_trinity), path(proteindb)

    output:
    path "${sample}.blastx.db.0.txt", emit: blastx0
    path "${sample}.blastx.db.6.txt", emit: blastx6

    script:
    """
   
    blastx -query ${combined_trinity} -db ${proteindbdbname} -out ${sample}.blastx.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"
    blastx -query ${combined_trinity} -db ${proteindbdbname} -out ${sample}.blastx.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """
}

// Process 10: Transdecoder
process Transdecoder {

    label 'process_medium'

    conda "transdecoder=5.7.1"
    container 'quay.io/biocontainers/transdecoder'


    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    publishDir "${sample}/Venomflow/results/ORFprediction/Transdecoder/", mode: 'copy'

    input:
    tuple val(sample), path(combined_trinity)

    output:
    tuple val(sample), path("*.pep"), emit: transdecoder_pep
    tuple val(sample), path("*.cds"), emit: transdecoder_cds

    script:

    """
    TransDecoder.LongOrfs -t ${combined_trinity} -m 15
    TransDecoder.Predict -t ${combined_trinity}

    """
}

// Process 10: TD2
process TD2 {

    label 'process_medium'

    conda "bioconda::td2"
    container 'quay.io/biocontainers/td2:1.0.7--pyhdfd78af_0'


    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    publishDir "${sample}/Venomflow/results/ORFprediction/TD2/", mode: 'copy'

    input:
    tuple val(sample), path(combined_trinity)

    output:
    tuple val(sample), path("*.pep"), emit: TD2_pep
    tuple val(sample), path("*.cds"), emit: TD2_cds

    script:

    """
    TD2.LongOrfs -t ${combined_trinity} -m 15 
    TD2.Predict -t ${combined_trinity} 

    """
}

// Process 6: Label transcriptomes, combine and remove duplicates 
process ORFs_Combined {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_medium'

    conda "seqkit=2.12.0"
    container "docker://gfanz/seqkit"

    publishDir "${sample}/Venomflow/results/ORFprediction/Combined", mode: 'copy'

    input:
    tuple val(sample), path(transdecoder_pep), path(transdecoder_cds), path(TD2_pep), path(TD2_cds)

    output:
    tuple val(sample), path("*combined.deduplicatedCDS.pep"), emit: combined_pep
    tuple val(sample), path("*combined.deduplicated.cds"), emit: combined_cds
    script:

    """
    seqkit replace -p '(.+)' -r 'TD_\$1' ${transdecoder_cds} > transdecoder_labelled.cds
    seqkit replace -p '(.+)' -r 'TD2_\$1' ${TD2_cds} > TD2_labelled.cds
    cat transdecoder_labelled.cds TD2_labelled.cds > orf_combined.cds
    seqkit rmdup orf_combined.cds -s -o ${sample}_ORF_combined.deduplicated.cds

    seqkit replace -p '(.+)' -r 'TD_\$1' ${transdecoder_pep} > transdecoder_labelled.pep
    seqkit replace -p '(.+)' -r 'TD2_\$1' ${TD2_pep} > TD2_labelled.pep
    cat transdecoder_labelled.pep TD2_labelled.pep > orf_combined.pep
    seqkit seq -n -i ${sample}_ORF_combined.deduplicated.cds > ids_from_cds.txt
    seqkit grep -f ids_from_cds.txt orf_combined.pep -o ${sample}_orf_combined.deduplicatedCDS.pep

    """
}

// Process 11: BUSCO_translatome_metazoa
process BUSCO_translatome_metazoa {

    label 'process_high'
    label 'process_long'

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${sample}/Venomflow/results/BUSCO/translatome/Transdecoder/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), val(metazoa)

    output:
    path "*.txt"

    script:

    """
    busco -i ${Transdecoder_pep} -l ${metazoa} -c 10 -o ${sample}_TD_met.protein -m protein -e 1e-5 -f
    mv ${sample}_TD_met.protein/*.txt "."
    
    """
}

// Process 12: BUSCO_translatome_mollusca
process BUSCO_translatome_mollusca {

    label 'process_high'
    label 'process_long'

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${sample}/Venomflow/results/BUSCO/translatome/Transdecoder/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), val(mollusca)

    output:
    path "*.txt"

    script:

    """
    busco -i ${Transdecoder_pep} -l ${mollusca} -c 10 -o ${sample}_TD.mol.protein -m protein -e 1e-5 -f
    mv ${sample}_TD.mol.protein/*.txt "."
    
    """
}

// Process 11: BUSCO_translatome_metazoa
process BUSCO_translatome_metazoa2 {

    label 'process_high'
    label 'process_long'

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${sample}/Venomflow/results/BUSCO/translatome/TD2/", mode: 'copy'

    input:
    tuple val(sample), path(TD2_pep), val(metazoa)

    output:
    path "*.txt"

    script:

    """
    busco -i ${TD2_pep} -l ${metazoa} -c 10 -o ${sample}_TD2_met.protein -m protein -e 1e-5 -f
    mv ${sample}_TD2_met.protein/*.txt "."
    
    """
}

// Process 12: BUSCO_translatome_mollusca
process BUSCO_translatome_mollusca2 {

    label 'process_high'
    label 'process_long'

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${sample}/Venomflow/results/BUSCO/translatome/Transdecoder/", mode: 'copy'

    input:
    tuple val(sample), path(TD2_pep), val(mollusca)

    output:
    path "*.txt"

    script:

    """
    busco -i ${TD2_pep} -l ${mollusca} -c 10 -o ${sample}_TD2.mol.protein -m protein -e 1e-5 -f
    mv ${sample}_TD2.mol.protein/*.txt "."
    
    """
}

// Process 11: BUSCO_translatome_metazoa
process BUSCO_translatome_metazoa3 {

    label 'process_high'
    label 'process_long'

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${sample}/Venomflow/results/BUSCO/translatome/Combined/", mode: 'copy'

    input:
    tuple val(sample), path(combined_pep), val(metazoa)

    output:
    path "*.txt"

    script:

    """
    busco -i ${combined_pep} -l ${metazoa} -c 10 -o ${sample}_combined_met.protein -m protein -e 1e-5 -f
    mv ${sample}_combined_met.protein/*.txt "."
    
    """
}

// Process 12: BUSCO_translatome_mollusca
process BUSCO_translatome_mollusca3 {

    label 'process_high'
    label 'process_long'

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "busco=5.8.3"
    container "docker://ezlabgva/busco:v5.8.2_cv1"

    publishDir "${sample}/Venomflow/results/BUSCO/translatome/Combined/", mode: 'copy'

    input:
    tuple val(sample), path(combined_pep), val(mollusca)

    output:
    path "*.txt"

    script:

    """
    busco -i ${combined_pep} -l ${mollusca} -c 10 -o ${sample}_combined.mol.protein -m protein -e 1e-5 -f
    mv ${sample}_combined.mol.protein/*.txt "."
    
    """
}


// Process 13: Kallisto_Transdecoder
process Kallisto_Transdecoder {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_medium'

    conda "kallisto=0.51.1"
    container "docker://quay.io/biocontainers/kallisto:0.51.1--h2b92561_2"

    publishDir "${sample}/Venomflow/results/kallisto/transdecoder/output", mode: 'copy'

    input:
    tuple val(sample), path(combined_cds), path(R1), path(R2), val(Strandedness)

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
    kallisto index -i index ${combined_cds}
    kallisto quant -i index -o ./ -b 100 ${R1} ${R2} \$stranded_flag

    """
}

// Process 14: Blastp
process Blastp {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_low'

    conda "blast=2.17.0"
    container "docker://ncbi/blast:2.17.0"

    publishDir "${sample}/Venomflow/results/Blast/Blastp/", mode: 'copy'

    input:
    tuple val(sample), path(combined_pep), path(proteindb), val(proteindbdbname)

    output:
    path "${sample}.blastp.db.0.txt", emit: blastp0
    path "${sample}.blastp.db.6.txt", emit: blastp6

    script:
    """

    blastp -query ${combined_pep} -db ${proteindbdbname} -out ${sample}.blastp.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"
    blastp -query ${combined_pep} -db ${proteindbdbname} -out ${sample}.blastp.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """
}

// Process 15: ORF_complete
process ORF_complete {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_single'

    conda "seqkit=2.12.0"
    container "docker://gfanz/seqkit"

    publishDir "${sample}/Venomflow/results/ORFprediction/Combined/Complete/", mode: 'copy'

    input:
    tuple val(sample), path(combined_pep), path(combined_cds)

    output:
    tuple val(sample), path("*.pep"), emit: complete_pep
    tuple val(sample), path("*.cds"), emit: complete_cds

    script:

    """
    seqkit grep -n -r -p "ORF type:complete" ${combined_pep} -o "${sample}.combine.complete.pep"
    seqkit grep -n -r -p "ORF type:complete" ${combined_cds} -o "${sample}.combine.complete.cds"
    """
}

// Process 16: SignalP
process SignalP {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_low'
    label 'process_long'

    publishDir "${sample}/Venomflow/results/Signalp", mode: 'copy'

    input:
    tuple val(sample), path(complete_pep)

    output:
    tuple val(sample), path("*_mature.fasta"), emit: maturesequences
    path "*.signalp5", emit: signalpsummary

    script:


    """
    signalp -fasta ${complete_pep} -mature -prefix "${sample}"
    """
}

// Process 17: Filter2
process Filter2 {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_single'

    conda "seqkit=2.12.0"
    container "docker://gfanz/seqkit"

    publishDir "${sample}/Venomflow/ORFprediction/Combined/", mode: 'copy'

    input:
    tuple val(sample), path(maturesequences), path(complete_pep)

    output:
    tuple val(sample), path("*.fasta"), emit: transdecoderpep_signalp

    script:

    """
    seqkit grep -f <(seqkit seq -n ${maturesequences}) ${complete_pep} > ${sample}.complete.signalp.sequences.fasta
    """
}

// Process 18: STATS
process stats {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_single'

    conda "seqkit=2.12.0"
    container "docker://gfanz/seqkit"

    publishDir "${sample}/Venomflow/results/Stats", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), path(Transdecoder_cds), path(transdecodercomplete_pep), path(transdecodercomplete_cds), path(maturesequences), path(transdecoderpep_signalp), path(TD2_pep), path(TD2_cds), path(combined_pep), path(combined_cds)

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
    seqkit stats ${TD2_cds} > ${sample}_TD2_cds_pep.stats.txt
    seqkit stats ${TD2_pep} > ${sample}_TD2_pep.stats.txt
    seqkit stats ${combined_cds} > ${sample}_combined_cds.stats.txt
    seqkit stats ${combined_pep} > ${sample}_combined_pep.stats.txt
    
    """
}

// Process 19: Interproscan
process Interproscan {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    label 'process_high'
    label 'process_long'

    publishDir "${sample}/Venomflow/results/Interproscan", mode: 'copy'

    input:
    tuple val(sample), path(complete_pep)

    output:
    path "*.tsv", emit: Interproscan

    script:

    """
    awk '{if (\$0 ~ /^>/) print \$0; else {gsub(/\\*/, ""); print \$0}}' ${complete_pep} > "${sample}.Trinity.fasta.transdecoder.cleaned.pep"

    interproscan.sh -goterms -i "${sample}.Trinity.fasta.transdecoder.cleaned.pep" -pa -t p -d ./ -f TSV 

    """
}

// Process 20: GenomeBlastdatabasecreation
process GenomeBlastdatabasecreation {

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    conda "blast=2.17.0"
    container "docker://ncbi/blast:2.17.0"

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
    container "docker://ncbi/blast:2.17.0"

    errorStrategy 'retry'
    maxRetries 4
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    publishDir "${sample}/Venomflow/results/Blast/Blastn/", mode: 'copy'

    input:
    tuple val(sample), val(genomedbname), path(completecds), path(genomedb)

    output:
    path "${sample}.blastn.db.0.txt", emit: blastn0
    path "${sample}.blastn.db.6.txt", emit: blastn6

    script:
    """
    blastn -query ${completecds} -db ${genomedbname} -out ${sample}.blastn.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"
    blastn -query ${completecds} -db ${genomedbname} -out ${sample}.blastn.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """
}


// PARAMETERS DEFINED IN CONFIG AND SAMPLESHEET FILE


//WorkFlow

workflow {

    // Print head 
    printhead()

    // Define CSV channel
    csv_channel = Channel.fromPath(params.input_csv)
        .splitCsv(header: true, sep: ',')
        .map { row -> row.collectEntries { key, value -> [key.replaceAll('"', ''), value?.toString()?.replaceAll('"', '')] } }

    // Define Input: Paired Trimmed reads Tuple. Extracts as tuple the sample name and the trimmed reads
    input_R1R2 = csv_channel.map { row -> tuple(row.Sample_name, file(row.R1), file(row.R2)) }

    // Run Process: Fastqc
    input_R1R2 | PostTrimFastqc

    //Define Input: Fastqc htmls. Using the output of Fastqc as an input for MultiQC
    fastqc_output = PostTrimFastqc.out.fastqc_zips


    //Run Process: MultQC
    fastqc_output | MultiQC

    // Define Input: Bowtie 
    input_bowtie = csv_channel.map { row -> tuple(row.Sample_name, file(row.Transcriptome1), file(row.R1), file(row.R2)) }

    //Run Process: Bowtie
    input_bowtie | Bowtie

    // Define Input: Bowtie2 
    input_bowtie2 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome2), file(row.R1), file(row.R2)) }

    //Run Process: Bowtie2
    input_bowtie2 | Bowtie2

    //Define Input: Trinity Fasta 
    input_trinity_fasta = csv_channel.map { row -> tuple(row.Sample_name, file(row.Transcriptome1)) }

    //Run Process: TrinityStats
    input_trinity_fasta | TrinityStats


    //Define Input: Trinity Fasta 
    input_trinity_fasta2 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome2)) }

    //Run Process: TrinityStats
    input_trinity_fasta2 | TrinityStats2

    // Transcriptome 1 
    //Define Input: Transcriptome1 Fasta + BUSCOlin1 tuple 
    input_BUSCOlin1 = csv_channel.map { row -> tuple(row.Sample_name, file(row.Transcriptome1), row.BUSCO_lin1, row.Transcriptome1_label) }

    //Run Process: BUSCO_lin1
    input_BUSCOlin1 | BUSCO_transcriptome_metazoa

    //Define Input: Transcriptome1 Fasta + BUSCOlin2 tuple 
    input_BUSCOlin2 = csv_channel.map { row -> tuple(row.Sample_name, file(row.Transcriptome1), row.BUSCO_lin2, row.Transcriptome1_label) }

    //Run Process: BUSCO_lin2
    input_BUSCOlin2 | BUSCO_transcriptome_mollusca

    // Transcriptome 2

    //Define Input: Transcriptome2 Fasta + BUSCOlin1 tuple 
    input_BUSCOlin1_2 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome2), row.BUSCO_lin1, row.Transcriptome2_label) }

    //Run Process: BUSCO_lin1
    input_BUSCOlin1_2 | BUSCO_transcriptome_metazoa2

    //Define Input: Transcriptome2 Fasta + BUSCOlin2 tuple 
    input_BUSCOlin2_2 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome2), row.BUSCO_lin2, row.Transcriptome2_label) }

    //Run Process: BUSCO_lin2
    input_BUSCOlin2_2 | BUSCO_transcriptome_mollusca2

    //Define Input:  Transcriptome_Combined 
    input_Transcriptome_Combined = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome1), row.Transcriptome1_label, file(row.Transcriptome2), row.Transcriptome2_label) }
    // Run process ; 
    input_Transcriptome_Combined | Transcriptome_Combined


    // Transcriptome_Combined
    //Define Input: Combined transcriptome Fasta + BUSCOlin1 tuple 
    input_BUSCOlin1_3 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, row.BUSCO_lin1) }
        .join(Transcriptome_Combined.out.transcriptome_combined)

    //Run Process: BUSCO_lin1
    input_BUSCOlin1_3 | BUSCO_transcriptome_metazoa3

    //Define Input: Combined transcriptome Fasta + BUSCOlin2 tuple 
    input_BUSCOlin2_3 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, row.BUSCO_lin2) }
        .join(Transcriptome_Combined.out.transcriptome_combined)

    //Run Process: BUSCO_lin2
    input_BUSCOlin2_3 | BUSCO_transcriptome_mollusca3






    // Define Input: Trinity fasta + R1 + R2 tuple 
    input_TrinityKallisto_single = csv_channel
        .filter { row -> !row.Transcriptome2?.trim() || row.Transcriptome2.trim().toLowerCase() == 'null' }
        .map { row -> tuple(row.Sample_name, file(row.R1), file(row.R2), row.Strandedness, file(row.Transcriptome1)) }

    input_TrinityKallisto_combined = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.R1), file(row.R2), row.Strandedness) }
        .join(Transcriptome_Combined.out.transcriptome_combined)

    // Combine both channels using mix() operator
    input_TrinityKallisto_all = input_TrinityKallisto_single.mix(input_TrinityKallisto_combined)

    // Run Process: Kallisto_Trinity
    input_TrinityKallisto_all | Kallisto_Trinity


    //Define Input: Database_Fasta. This is set up to allow for multiple different databases if needed.
    input_databasefasta = csv_channel.map { row -> tuple(row.Sample_name, file(row.Protein_fasta_path_for_Blast)) }

    //Run Process: DatabaseCreation
    input_databasefasta | Blastdatabasecreation


    //Define Input: Blastx 
    Blastxinputfasta_single = csv_channel
        .filter { row -> !row.Transcriptome2?.trim() || row.Transcriptome2.trim().toLowerCase() == 'null' }
        .map { row ->
        tuple(row.Sample_name, row.Protein_fasta_name, file(row.Transcriptome1))
    }
        .join(Blastdatabasecreation.out.proteindb)

    Blastxinputfasta_combined = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row ->
            tuple(row.Sample_name, row.Protein_fasta_name)
        }
        .join(Transcriptome_Combined.out.transcriptome_combined)
        .join(Blastdatabasecreation.out.proteindb)

    // Combine both channels using mix() operator
    Blastxinputfasta_all = Blastxinputfasta_single.mix(Blastxinputfasta_combined)

    // Run Process: Blastx 
    Blastxinputfasta_all | Blastx


    //Run Process: Transdecoder
    Transcriptpome1 = csv_channel
        .filter { row -> !row.Transcriptome2?.trim() || row.Transcriptome2.trim().toLowerCase() == 'null' }
        .map { row ->
            tuple(row.Sample_name, file(row.Transcriptome1))
        }

    input_orf = Transcriptpome1.mix(Transcriptome_Combined.out.transcriptome_combined)
    input_orf | Transdecoder


    //Run Process: TD2
    input_orf | TD2

    // Define input: ORFs_Combined
    input_ORFs_Combined = Transdecoder.out.transdecoder_pep.join(Transdecoder.out.transdecoder_cds).join(TD2.out.TD2_pep).join(TD2.out.TD2_cds)
    // Run Process: ORFs_Combined
    input_ORFs_Combined | ORFs_Combined

    // Transdecoder
    //Define Input: Transdecoder pep + BUSCOlin1 tuple 
    Transdecoder_pep = Transdecoder.out.transdecoder_pep
    BUSCOlin1 = csv_channel.map { row -> tuple(row.Sample_name, row.BUSCO_lin1) }
    input_BUSCOlin1_L = Transdecoder_pep.join(BUSCOlin1)

    //Run Process: BUSCO_lin1
    input_BUSCOlin1_L | BUSCO_translatome_metazoa

    //Define Input: Transdecoder pep + BUSCOlin2 tuple 
    BUSCOlin2 = csv_channel.map { row -> tuple(row.Sample_name, row.BUSCO_lin2) }
    input_BUSCOlin2_L = Transdecoder_pep.join(BUSCOlin2)

    //Run Process: BUSCO_lin2
    input_BUSCOlin2_L | BUSCO_translatome_mollusca

    // TD2
    //Define Input: Transdecoder pep + BUSCOlin1 tuple 
    TD2_pep = TD2.out.TD2_pep
    input_BUSCOlin1_L_2 = TD2_pep.join(BUSCOlin1)
    //Run Process: BUSCO_lin1
    input_BUSCOlin1_L_2 | BUSCO_translatome_metazoa2
    //Define Input: Transdecoder pep + BUSCOlin2 tuple 
    BUSCOlin2 = csv_channel.map { row -> tuple(row.Sample_name, row.BUSCO_lin2) }
    input_BUSCOlin2_L_2 = TD2_pep.join(BUSCOlin2)
    //Run Process: BUSCO_lin2
    input_BUSCOlin2_L_2 | BUSCO_translatome_mollusca2


    // Combined
    //Define Input: Transdecoder pep + BUSCOlin1 tuple 
    Combined_pep = ORFs_Combined.out.combined_pep
    input_BUSCOlin1_L_3 = Combined_pep.join(BUSCOlin1)
    //Run Process: BUSCO_lin1
    input_BUSCOlin1_L_3 | BUSCO_translatome_metazoa3
    //Define Input: Transdecoder pep + BUSCOlin2 tuple 
    BUSCOlin2 = csv_channel.map { row -> tuple(row.Sample_name, row.BUSCO_lin2) }
    input_BUSCOlin2_L_3 = Combined_pep.join(BUSCOlin2)
    //Run Process: BUSCO_lin2
    input_BUSCOlin2_L_3 | BUSCO_translatome_mollusca3


    //Define Input: Transdecoder cds + R1 + R2 + Strandedness tuple 
    Combined_cds = ORFs_Combined.out.combined_cds
    KallistoTransdecoderR1R2S = csv_channel.map { row -> tuple(row.Sample_name, file(row.R1), file(row.R2), row.Strandedness) }
    input_TransKallisto = Combined_cds.join(KallistoTransdecoderR1R2S)

    //Run Process: TransKallisto
    input_TransKallisto | Kallisto_Transdecoder

    //Define Input: Blastp - Match Transdecoder output with databases
    Blastpdb = csv_channel.map { row ->
        tuple(row.Sample_name, row.Protein_fasta_name)
    }
    input_Blastp = Combined_pep.join(Blastdatabasecreation.out.proteindb).join(Blastpdb)

    //Run Process: Blastp
    input_Blastp | Blastp

    //Define Input: Sample name + PEP + CDS tuple
    input_ORF_complete = Combined_pep.join(Combined_cds)

    //Run Process: Transdecoder filter for complete ORFs
    input_ORF_complete | ORF_complete


    //Define Input: Sample name + completepep tuple
    input_signalp = ORF_complete.out.complete_pep

    //Run Process: SignalP
    input_signalp | SignalP

    //Define Input: sample name + mature sequences + completepep tuple 
    maturesequences = SignalP.out.maturesequences
    maturecomplete = maturesequences.join(input_signalp)

    //Run Process: Filter2   
    maturecomplete | Filter2

    //Define Input: stats sample name + pep + cds + completepep + completecds + mature +signalp  tuple 
    stats_join = Transdecoder_pep
        .join(Transdecoder.out.transdecoder_cds)
        .join(ORF_complete.out.complete_pep)
        .join(ORF_complete.out.complete_cds)
        .join(maturesequences)
        .join(Filter2.out.transdecoderpep_signalp)
        .join(TD2_pep)
        .join(TD2.out.TD2_cds)
        .join(input_ORF_complete)
    //Run Process: STATS   
    stats_join | stats

    //Run Process: Interproscan  (Complete pep ORFs only)
    input_signalp | Interproscan

    //Define Input: genomefasta 
    Genomefasta = csv_channel
        .filter { row ->
            def path = row.Genome_fasta_path
            path != null && path != 'NULL' && path.toString() != 'NULL' && path.toString().trim() != ''
        }
        .map { row -> tuple(row.Sample_name, file(row.Genome_fasta_path)) }
        .unique { sample_name, genome_file -> sample_name }

    // Run Process: BlastnGenome database creation
    GenomeBlastdatabasecreation(Genomefasta)

    //Define Input: Genome BLAST 
    genomedb = GenomeBlastdatabasecreation.out.genomedbfiles
    Blastncds = csv_channel.map { row -> tuple(row.Sample_name, row.Genome_fasta_name) }

    BlastnInput = Blastncds
        .join(ORF_complete.out.complete_cds)
        .join(genomedb)

    //Run Process: BlastnGenome
    BlastnInput | GenomeBlasts
}
