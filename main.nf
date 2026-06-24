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
    log.info("DeepTMHMM:                       ${params.DeepTMHMM}")
    log.info("Profile:                         ${workflow.profile}")
    log.info("ORFPrediction:                   ${params.ORFPrediction}")
    log.info("Start:                           ${workflow.start}")
    log.info("")
    log.info("────────────────────────────────────────────────────────────")
    log.info("")
}

// Process 1: PostTrimFastqc
process PostTrimFastqc {


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5



    label 'process_low'

    conda "fastqc=0.12.1"
    container 'community.wave.seqera.io/library/fastqc:0.12.1--aa717e1a9d994d74'

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
    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'
    conda "multiqc=1.33"
    container 'community.wave.seqera.io/library/multiqc:1.33--9daaf37cc59ba7dc'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Fastqc/posttrim/", mode: 'copy'
    publishDir "${params.outdir}/${sample}/FinalOutputs/htmls/", mode: 'copy'

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


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    conda "bowtie2=2.5.4"
    container 'community.wave.seqera.io/library/bowtie2:2.5.4--d5022d6316284d3d'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Bowtie/Transcriptome1", pattern: "*.log", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), path(R1), path(R2)

    output:
    path "*.log"

    script:

    """
    bowtie2-build ${trinity_fasta} ${sample}_transcriptome_index -p ${task.cpus}
    bowtie2 -x ${sample}_transcriptome_index -1 ${R1} -2 ${R2} -S ${sample}_mapped_reads.sam -p ${task.cpus} --no-unal 2> stats.log 



    """
}
// Process 3: Bowtie build and index # Transcriptome 1 

process Bowtie2 {

    label 'process_high'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5





    conda "bowtie2=2.5.4"
    container 'community.wave.seqera.io/library/bowtie2:2.5.4--d5022d6316284d3d'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Bowtie/Transcriptome2", pattern: "*.log", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta2), path(R1), path(R2)

    output:
    path "*.log"

    script:

    """
    bowtie2-build ${trinity_fasta2} ${sample}_transcriptome_index -p ${task.cpus}
    bowtie2 -x ${sample}_transcriptome_index -1 ${R1} -2 ${R2} -S ${sample}_mapped_reads.sam -p ${task.cpus} --no-unal 2> stats.log



    """
}


// Process 4: TrinityStats 1 
process TrinityStats {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5



    label 'process_low'

    conda "seqkit=2.12.0"
    container 'community.wave.seqera.io/library/seqkit:2.12.0--430b52150147f163'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Stats/", mode: 'copy'

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

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'


    conda "seqkit=2.12.0"
    container 'community.wave.seqera.io/library/seqkit:2.12.0--430b52150147f163'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Stats/", mode: 'copy'

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

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_medium'


    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/transcriptome/Transcriptome1", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), val(metazoa), val(transcriptome1_label)

    output:
    path "*.txt"

    script:

    """
    busco -i ${trinity_fasta} -l ${metazoa} -c ${task.cpus} -o ${sample}_${transcriptome1_label}_met.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_${transcriptome1_label}_met.transcriptome/*.txt "."
    """
}

// Process 6: BUSCO_transcriptome_mollusca
process BUSCO_transcriptome_mollusca {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_medium'





    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/transcriptome/Transcriptome1", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), val(mollusca), val(transcriptome1_label)

    output:
    path "*.txt"

    script:

    """
    busco -i ${trinity_fasta} -l ${mollusca} -c ${task.cpus} -o ${sample}_${transcriptome1_label}_mol.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_${transcriptome1_label}_mol.transcriptome/*.txt "."
    """
}

// Process 5: BUSCO_transcriptome_metazoa
process BUSCO_transcriptome_metazoa2 {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_medium'


    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'


    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/transcriptome/Transcriptome2", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta2), val(metazoa), val(transcriptome2_label)

    output:
    path "*.txt"

    script:

    """
    busco -i ${trinity_fasta2} -l ${metazoa} -c ${task.cpus} -o ${sample}_${transcriptome2_label}_met.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_${transcriptome2_label}_met.transcriptome/*.txt "."
    """
}

// Process 6: BUSCO_transcriptome_mollusca
process BUSCO_transcriptome_mollusca2 {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_medium'


    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/transcriptome/Transcriptome2", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta2), val(mollusca), val(transcriptome2_label)

    output:
    path "*.txt"

    script:

    """
    busco -i ${trinity_fasta2} -l ${mollusca} -c ${task.cpus} -o ${sample}_${transcriptome2_label}_mol.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_${transcriptome2_label}_mol.transcriptome/*.txt "."
    """
}

// Process 6: Label transcriptomes, combine and remove duplicates 
process Transcriptome_Combined {


    label 'process_low'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5

    conda "seqkit=2.12.0 bioconda::cd-hit=4.8.1"
    container 'community.wave.seqera.io/library/cd-hit_seqkit:27b33ce1ba0d851c'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Combined_Transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(transcriptome1), val(transcriptome1_label), path(transcriptome2), val(transcriptome2_label)

    output:
    tuple val(sample), path("*combined.deduplicated.fasta"), emit: transcriptome_combineddedup
    tuple val(sample), path("*deduplicated.cdhit95.fasta"), emit: transcriptome_combined

    script:

    """
    seqkit replace -p '(.+)' -r '${transcriptome1_label}_\$1' ${transcriptome1} > Transcriptome1_labelled.fasta
    seqkit replace -p '(.+)' -r '${transcriptome2_label}_\$1' ${transcriptome2} > Transcriptome2_labelled.fasta
    cat Transcriptome1_labelled.fasta Transcriptome2_labelled.fasta > Transcriptome_combined.fasta
    seqkit rmdup Transcriptome_combined.fasta -s -o ${sample}_transcriptome_combined.deduplicated.fasta
    cd-hit-est -i ${sample}_transcriptome_combined.deduplicated.fasta  -o ${sample}_transcriptome_combined.deduplicated.cdhit95.fasta -c 0.95 -d 0 -aS 0.9 -aL 0.9 -M 2865

    """
}


// Process 5: BUSCO_transcriptome_metazoa Combined transcriptome
process BUSCO_transcriptome_metazoa3 {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_medium'





    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'


    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/transcriptome/Combined", mode: 'copy'

    input:
    tuple val(sample), val(metazoa), path(combined_trinity)

    output:
    path "*.txt"

    script:

    """
    busco -i ${combined_trinity} -l ${metazoa} -c ${task.cpus} -o ${sample}_combined_met.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_combined_met.transcriptome/*.txt "."
    """
}

// Process 6: BUSCO_transcriptome_mollusca
process BUSCO_transcriptome_mollusca3 {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5



    label 'process_medium'


    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/transcriptome/Combined", mode: 'copy'

    input:
    tuple val(sample), val(mollusca), path(combined_trinity)

    output:
    path "*.txt"

    script:

    """
    busco -i ${combined_trinity} -l ${mollusca} -c ${task.cpus} -o ${sample}_combined_mol.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_combined_mol.transcriptome/*.txt "."
    """
}

// Process 7: Kallisto_Trinity
process Kallisto_Trinity {

    conda "kallisto=0.51.1"
    container 'community.wave.seqera.io/library/kallisto:0.51.1--d7728813dda40c70'


    label 'process_low'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/kallisto/trinity/output", mode: 'copy'

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5

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
    kallisto index -i index ${combined_trinity} -t ${task.cpus}
    kallisto quant -i index -o ./ -b 100 ${R1} ${R2} \$stranded_flag -t ${task.cpus}

    """
}

// Process 8: Blastdatabasecreation
process Blastdatabasecreation {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'

    conda "blast=2.17.0"
    container 'community.wave.seqera.io/library/blast:2.17.0--6279aeee601cb05e'

    input:
    path database_fasta

    output:
    path "ToxProtdb", emit: proteindb

    script:
    """
    makeblastdb -in "${database_fasta}" -dbtype prot -out ToxProtdb/proteindb
    """
}
// Process 9: Blastx
process Blastx {

    label 'process_low'

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5




    conda "blast=2.17.0"
    container 'community.wave.seqera.io/library/blast:2.17.0--6279aeee601cb05e'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Blast/Blastx/", mode: 'copy'

    input:
    tuple val(sample), path(combined_trinity), path(proteindb)

    output:
    path "${sample}.blastx.db.0.txt", emit: blastx0
    path "${sample}.blastx.db.6.txt", emit: blastx6

    script:
    """
   
    blastx -query ${combined_trinity} -db ToxProtdb/proteindb -out ${sample}.blastx.db.6.txt -evalue 1e-5 -num_threads ${task.cpus} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"
    blastx -query ${combined_trinity} -db ToxProtdb/proteindb -out ${sample}.blastx.db.0.txt -evalue 1e-5 -num_threads ${task.cpus} -outfmt '0'

    """
}

// Process 10: Transdecoder
process Transdecoder {

    label 'process_varied'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    conda "transdecoder=5.7.1"
    container 'community.wave.seqera.io/library/transdecoder:5.7.1--bfc613e7081a52d9'



    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/ORFprediction/Transdecoder/", mode: 'copy'

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

    label 'process_varied'



    conda "bioconda::td2"
    container 'community.wave.seqera.io/library/td2:ca3786e862ccfcd7'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5





    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/ORFprediction/TD2/", mode: 'copy'

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

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'

    conda "seqkit=2.12.0"
    container 'community.wave.seqera.io/library/seqkit:2.12.0--430b52150147f163'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/ORFprediction/Combined/All", mode: 'copy'

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

// Process 6: remove duplicates for those with no genomes
process ORFs_Combined_CDHit {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'

    conda "seqkit=2.12.0 bioconda::cd-hit=4.8.1"
    container 'community.wave.seqera.io/library/cd-hit_seqkit:27b33ce1ba0d851c'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/ORFprediction/Combined/All", mode: 'copy'

    input:
    tuple val(sample), path(combined_pep), path(combined_cds)

    output:
    tuple val(sample), path("*combined.deduplicatedcds.cd95.pep"), emit: combined_pep
    tuple val(sample), path("*combined.deduplicated.cd95pep.cds"), emit: combined_cds

    script:

    """
    cd-hit -i ${combined_pep}  -o ${sample}_combined.deduplicatedcds.cd95.pep -c 0.95 -d 0 -aS 0.9 -aL 0.9
    seqkit seq -n -i ${sample}_combined.deduplicatedcds.cd95.pep > ids_from_pep.txt
    seqkit grep -f ids_from_pep.txt ${combined_cds} -o ${sample}_combined.deduplicated.cd95pep.cds

    """
}

// Process 11: BUSCO_translatome_metazoa
process BUSCO_translatome_metazoa {

    label 'process_medium'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5





    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/translatome/Transdecoder/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), val(metazoa)

    output:
    path "*.txt"

    script:

    """
    busco -i ${Transdecoder_pep} -l ${metazoa} -c ${task.cpus} -o ${sample}_TD_met.protein -m protein -e 1e-5 -f
    mv ${sample}_TD_met.protein/*.txt "."
    
    """
}

// Process 12: BUSCO_translatome_mollusca
process BUSCO_translatome_mollusca {

    label 'process_medium'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5





    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/translatome/Transdecoder/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), val(mollusca)

    output:
    path "*.txt"

    script:

    """
    busco -i ${Transdecoder_pep} -l ${mollusca} -c ${task.cpus} -o ${sample}_TD_mol.protein -m protein -e 1e-5 -f
    mv ${sample}_TD_mol.protein/*.txt "."
    
    """
}

// Process 11: BUSCO_translatome_metazoa
process BUSCO_translatome_metazoa2 {

    label 'process_medium'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5





    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/translatome/TD2/", mode: 'copy'

    input:
    tuple val(sample), path(TD2_pep), val(metazoa)

    output:
    path "*.txt"

    script:

    """
    busco -i ${TD2_pep} -l ${metazoa} -c ${task.cpus} -o ${sample}_TD2_met.protein -m protein -e 1e-5 -f
    mv ${sample}_TD2_met.protein/*.txt "."
    
    """
}

// Process 12: BUSCO_translatome_mollusca
process BUSCO_translatome_mollusca2 {

    label 'process_medium'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5





    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/translatome/TD2/", mode: 'copy'

    input:
    tuple val(sample), path(TD2_pep), val(mollusca)

    output:
    path "*.txt"

    script:

    """
    busco -i ${TD2_pep} -l ${mollusca} -c ${task.cpus} -o ${sample}_TD2_mol.protein -m protein -e 1e-5 -f
    mv ${sample}_TD2_mol.protein/*.txt "."
    
    """
}

// Process 11: BUSCO_translatome_metazoa
process BUSCO_translatome_metazoa3 {

    label 'process_medium'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5





    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/translatome/Combined/", mode: 'copy'

    input:
    tuple val(sample), path(combined_pep), val(metazoa)

    output:
    path "*.txt"

    script:

    """
    busco -i ${combined_pep} -l ${metazoa} -c ${task.cpus} -o ${sample}_combined_met.protein -m protein -e 1e-5 -f
    mv ${sample}_combined_met.protein/*.txt "."
    
    """
}

// Process 12: BUSCO_translatome_mollusca
process BUSCO_translatome_mollusca3 {

    label 'process_medium'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5





    conda "busco=5.8.3"
    container 'community.wave.seqera.io/library/busco:5.8.3--dac4836fc2571f70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/BUSCO/translatome/Combined/", mode: 'copy'

    input:
    tuple val(sample), path(combined_pep), val(mollusca)

    output:
    path "*.txt"

    script:

    """
    busco -i ${combined_pep} -l ${mollusca} -c ${task.cpus} -o ${sample}_combined_mol.protein -m protein -e 1e-5 -f
    mv ${sample}_combined_mol.protein/*.txt "."
    
    """
}


// Process 13: Kallisto_Transdecoder
process Kallisto_Transdecoder {



    label 'process_low'


    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5




    conda "kallisto=0.51.1"
    container 'community.wave.seqera.io/library/kallisto:0.51.1--d7728813dda40c70'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/kallisto/transdecoder/output", mode: 'copy'

    input:
    tuple val(sample), path(complete_cds), path(R1), path(R2), val(Strandedness)

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
    kallisto index -i index ${complete_cds} -t ${task.cpus}
    kallisto quant -i index -o ./ -b 100 ${R1} ${R2} \$stranded_flag -t ${task.cpus}

    """
}

// Process 14: Blastp
process Blastp {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'




    conda "blast=2.17.0"
    container 'community.wave.seqera.io/library/blast:2.17.0--6279aeee601cb05e'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Blast/Blastp_Toxin/", mode: 'copy'

    input:
    tuple val(sample), path(combined_pep), path(proteindb)

    output:
    path "${sample}.blastp.db.0.txt", emit: blastp0
    path "${sample}.blastp.db.6.txt", emit: blastp6

    script:
    """

    blastp -query ${combined_pep} -db ToxProtdb/proteindb -out ${sample}.blastp.db.6.txt -evalue 1e-5 -num_threads ${task.cpus} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"
    blastp -query ${combined_pep} -db ToxProtdb/proteindb -out ${sample}.blastp.db.0.txt -evalue 1e-5 -num_threads ${task.cpus} -outfmt '0'

    """
}

// Process 15: ORF_complete
process ORF_complete {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'




    conda "seqkit=2.12.0"
    container 'community.wave.seqera.io/library/seqkit:2.12.0--430b52150147f163'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/ORFprediction/Combined/Complete/", pattern: "*cleaned*", mode: 'copy'

    input:
    tuple val(sample), path(combined_pep), path(combined_cds)

    output:
    tuple val(sample), path("*cleaned.pep"), emit: complete_pep
    tuple val(sample), path("*cleaned.cds"), emit: complete_cds

    script:

    """
    seqkit grep -n -r -p "ORF type:complete" ${combined_pep} -o "${sample}.combine.complete.pep"
    seqkit grep -n -r -p "ORF type:complete" ${combined_cds} -o "${sample}.combine.complete.cds"
    awk '{if (\$0 ~ /^>/) print \$0; else {gsub(/\\*/, ""); print \$0}}' "${sample}.combine.complete.pep" > "${sample}.combine.complete.cleaned.pep"
    awk '{if (\$0 ~ /^>/) print \$0; else {gsub(/\\*/, ""); print \$0}}' "${sample}.combine.complete.cds" > "${sample}.combine.complete.cleaned.cds"

    """
}

// Process 16: SignalP
process SignalP {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_medium'
    conda "seqkit=2.12.0"


    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Secreted/Mature/Signalp", mode: 'copy'

    input:
    tuple val(sample), path(complete_pep)

    output:
    tuple val(sample), path("${sample}_mature.fasta"), emit: maturesequences
    tuple val(sample), path("${sample}_summary.signalp5"), emit: signalpsummary

    script:

    """
    seqkit split ${complete_pep} -s 10000
    for file in  "${complete_pep}".split/*
    do
        signalp -fasta "\$file" -mature
    done

    cat *mature.fasta > ${sample}_mature.fasta
    cat *.signalp5 > ${sample}_summary.signalp5
    """
}

// Process 17: Filter2
process Filter2 {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'

    conda "seqkit=2.12.0"
    container 'community.wave.seqera.io/library/seqkit:2.12.0--430b52150147f163'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Secreted/Full_Secreted/Signalp", mode: 'copy'

    input:
    tuple val(sample), path(maturesequences), path(complete_pep), path(complete_cds)

    output:
    tuple val(sample), path("*.pep.fasta"), emit: complete_pep_signalp
    tuple val(sample), path("*.cds.fasta"), emit: complete_cds_signalp

    script:

    """
    seqkit grep -f <(seqkit seq -n ${maturesequences}) ${complete_pep} > ${sample}.complete.signalp.sequences.pep.fasta
    seqkit grep -f <(seqkit seq -n ${maturesequences}) ${complete_cds} > ${sample}.complete.signalp.sequences.cds.fasta

    """
}

// Process 18: STATS
process stats {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'

    conda "seqkit=2.12.0"
    container 'community.wave.seqera.io/library/seqkit:2.12.0--430b52150147f163'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Stats", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), path(Transdecoder_cds), path(complete_pep), path(complete_cds), path(maturesequences), path(complete_pep_signalp), path(TD2_pep), path(TD2_cds), val(combined_pep), val(combined_cds)

    output:
    path "*"

    script:

    """
    seqkit stats ${Transdecoder_pep} > ${sample}_Transdecoder_pep.stats.txt
    seqkit stats ${Transdecoder_cds} > ${sample}_Transdecoder_cds.stats.txt
    seqkit stats ${complete_pep} > ${sample}_complete_pep.stats.txt
    seqkit stats ${complete_cds} > ${sample}_complete_cds.stats.txt
    seqkit stats ${maturesequences} > ${sample}_maturesequences.stats.txt
    seqkit stats ${complete_pep_signalp} > ${sample}_complete_pep_signalp.stats.txt
    seqkit stats ${TD2_cds} > ${sample}_TD2_cds_pep.stats.txt
    seqkit stats ${TD2_pep} > ${sample}_TD2_pep.stats.txt
    # Check for optional inputs
    if [ "${combined_cds}" != "NULL" ]; then
        seqkit stats ${combined_cds} > ${sample}_combined_cds.stats.txt
        seqkit stats ${combined_pep} > ${sample}_combined_pep.stats.txt

    fi
    
    """
}

// Process 18a:  DeepTMHMM (Optional)

process DeepTMHMM {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5

    label 'process_medium'

    conda "seqkit=2.12.0"
    container 'community.wave.seqera.io/library/seqkit:2.12.0--430b52150147f163'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Secreted/Mature/DeepTMHMM", pattern: "*min5.fasta", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Stats", pattern: "*.txt", mode: 'copy'

    input:
    tuple val(sample), path(complete_pep), path(signalpsummary)

    output:
    tuple val(sample), path('*min5.fasta'), emit: mature
    tuple val(sample), path('*.txt')

    script:
    """
    RUN_DIR="\$(pwd)"
    awk -F'\t' 'NR==2{for(i=1;i<=NF;i++)if(\$i=="SP(Sec/SPI)")col=i} NR>2 && \$col>=0.25 && \$col<=0.5{print \$1}' ${signalpsummary} > labels.txt
    seqkit grep -f labels.txt ${complete_pep} -o deeptmhmmcandidates.pep
    awk '{if (\$0 ~ /^>/) print \$0; else {gsub(/\\*/, ""); print \$0}}' deeptmhmmcandidates.pep > deeptmhmmcandidates.cleaned.pep

    predict --fasta \$RUN_DIR/deeptmhmmcandidates.cleaned.pep --output-dir \$RUN_DIR/${sample}
    cd \$RUN_DIR
    python3 ${workflow.projectDir}/bin/deepout.py ${sample}/predicted_topologies.3line ${sample}_Deep_mature_sequences.fasta 
    seqkit seq -m 5 ${sample}_Deep_mature_sequences.fasta > ${sample}_Deep_mature_sequences.min5.fasta 
    seqkit stats ${sample}_Deep_mature_sequences.min5.fasta > ${sample}_Deep_mature_sequences.stats.txt
    """
}

// Process 18b:  DeepTMHMMfiler (Optional)

process DeepTMHMMFilter {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5

    label 'process_low'



    conda "seqkit=2.12.0"
    container 'community.wave.seqera.io/library/seqkit:2.12.0--430b52150147f163'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Secreted/Full_Secreted/DeepTMHMM", pattern: "*DeepTMHMM*", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Secreted/Full_Secreted/Combined", pattern: "*sequences.deduplicated*", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Stats", pattern: "*.txt", mode: 'copy'
    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Secreted/Mature/Combined/", pattern: "*.combined.mature.deduplicated.pep.fasta", mode: 'copy'

    input:
    tuple val(sample), path(DeepTMHMM_mature), path(complete_pep), path(complete_cds), path(complete_pep_signalp), path(complete_cds_signalp), path(signalpmature)

    output:
    tuple val(sample), path('*DeepTMHMM.sequences.pep.fasta')
    tuple val(sample), path('*DeepTMHMM.sequences.cds.fasta')
    tuple val(sample), path('*.combined.mature.deduplicated.pep.fasta')
    tuple val(sample), path('*sequences.deduplicated.pep.fasta'), emit: complete_pep_secreted
    tuple val(sample), path('*sequences.deduplicated.cds.fasta'), emit: complete_cds_secreted
    tuple val(sample), path('*.txt')

    script:
    """
    seqkit grep -f <(seqkit seq -n ${DeepTMHMM_mature}) ${complete_pep} > ${sample}.complete.DeepTMHMM.sequences.pep.fasta
    seqkit grep -f <(seqkit seq -n ${DeepTMHMM_mature}) ${complete_cds} > ${sample}.complete.DeepTMHMM.sequences.cds.fasta
    cat ${sample}.complete.DeepTMHMM.sequences.pep.fasta ${complete_pep_signalp} > ${sample}.complete.secreted.sequences.pep.fasta
    seqkit rmdup -n ${sample}.complete.secreted.sequences.pep.fasta > ${sample}.complete.secreted.sequences.deduplicated.pep.fasta
    cat ${sample}.complete.DeepTMHMM.sequences.cds.fasta ${complete_cds_signalp} > ${sample}.complete.secreted.sequences.cds.fasta
    seqkit rmdup -n ${sample}.complete.secreted.sequences.cds.fasta > ${sample}.complete.secreted.sequences.deduplicated.cds.fasta
    cat ${DeepTMHMM_mature} ${signalpmature} > ${sample}.combined.mature.pep.fasta
    seqkit rmdup -n ${sample}.combined.mature.pep.fasta > ${sample}.combined.mature.deduplicated.pep.fasta

    seqkit stats ${sample}.complete.DeepTMHMM.sequences.pep.fasta > ${sample}.complete.Deep.sequences.pep.stats.txt
    seqkit stats ${sample}.complete.DeepTMHMM.sequences.cds.fasta > ${sample}.complete.Deep.sequences.cds.stats.txt
    seqkit stats ${sample}.complete.secreted.sequences.deduplicated.pep.fasta > ${sample}.complete.secreted.deduplicated.pep.stats.txt
    seqkit stats ${sample}.complete.secreted.sequences.deduplicated.cds.fasta > ${sample}.complete.secreted.deduplicated.cds.stats.txt
    seqkit stats ${sample}.combined.mature.deduplicated.pep.fasta > ${sample}.combined.mature.deduplicated.pep.stats.txt
    



    """
}
// Process 19: Interproscan
process Interproscan {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5

    label 'process_medium'


    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Interproscan", mode: 'copy'

    input:
    tuple val(sample), path(secreted_pep)

    output:
    path "*.tsv", emit: Interproscan

    script:

    """


    awk '{if (\$0 ~ /^>/) print \$0; else {gsub(/\\*/, ""); print \$0}}' ${secreted_pep} > "${sample}.Trinity.fasta.secreted.cleaned.pep"

    interproscan.sh -goterms -i "${sample}.Trinity.fasta.secreted.cleaned.pep" -pa -t p -d ./ -f TSV 

    """
}

// Process 20: Blastdatabasecreation
process BlastdatabasecreationNonToxin {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'

    conda "blast=2.17.0"
    container 'community.wave.seqera.io/library/blast:2.17.0--6279aeee601cb05e'

    input:
    path database_fasta

    output:
    path "Folder", emit: nontoxinproteindb

    script:
    """
    makeblastdb -in "${database_fasta}" -dbtype prot -out Folder/NonToxinDataBase
    """
}
// Process 21: Blastdatabasecreation

process BlastpNonToxin {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    label 'process_low'



    conda "blast=2.17.0"
    container 'community.wave.seqera.io/library/blast:2.17.0--6279aeee601cb05e'

    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Blast/Blastp_NonToxin/", mode: 'copy'

    input:
    tuple val(sample), path(secreted_pep), path(proteindb)

    output:
    tuple val(sample), path("${sample}.nontoxin.blastp.db.6.txt"), emit: nontoxblastp6

    script:
    """

    blastp -query ${secreted_pep} -db Folder/NonToxinDataBase -out ${sample}.nontoxin.blastp.db.6.txt -evalue 1e-5 -num_threads ${task.cpus} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"

    """
}

// Process 22: GenomeBlastdatabasecreation
process GenomeBlastdatabasecreation {

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    conda "blast=2.17.0"
    container 'community.wave.seqera.io/library/blast:2.17.0--6279aeee601cb05e'

    label 'process_low'

    input:
    tuple val(sample), path(genome_fasta)

    output:
    tuple val(sample), path("*"), emit: genomedbfiles

    script:
    """
    makeblastdb -in "${genome_fasta}" -dbtype nucl -out genomedb
    """
}

// Process 23: GenomeBlasts6
process GenomeBlasts6 {

    label 'process_medium'


    conda "blast=2.17.0"
    container 'community.wave.seqera.io/library/blast:2.17.0--6279aeee601cb05e'

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5


    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Blast/Blastn/", mode: 'copy'

    input:
    tuple val(sample), path(secretedcds), path(genomedb)

    output:
    tuple val(sample), path("${sample}.blastn.db.6.txt"), emit: blastn6


    script:
    """
    blastn -query ${secretedcds} -db genomedb -out ${sample}.blastn.db.6.txt -evalue 1e-5 -num_threads ${task.cpus} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qframe qcovs"

    """
}

// Process 21: GenomeBlasts0
process GenomeBlasts0 {

    label 'process_medium'


    conda "blast=2.17.0"
    container 'community.wave.seqera.io/library/blast:2.17.0--6279aeee601cb05e'

    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5




    publishDir "${params.outdir}/${sample}/Pipelines/Venomflow/Blast/Blastn/", mode: 'copy'

    input:
    tuple val(sample), path(secretedcds), path(genomedb)

    output:
    tuple val(sample), path("${sample}.blastn.db.0.txt"), emit: blastn0

    script:
    """
    blastn -query ${secretedcds} -db genomedb -out ${sample}.blastn.db.0.txt -evalue 1e-5 -num_threads ${task.cpus}  -outfmt '0'

    """
}

// PARAMETERS INED IN CONFIG AND SAMPLESHEET FILE


//WorkFlow

workflow {

    // Print head 
    printhead()

    //Fasta Databases



    // ine CSV channel
    csv_channel = Channel.fromPath(params.input_csv)
        .splitCsv(header: true, sep: ',')
        .map { row -> row.collectEntries { key, value -> [key.replaceAll('"', ''), value?.toString()?.replaceAll('"', '')] } }


    // ine Input: Paired Trimmed reads Tuple. Extracts as tuple the sample name and the trimmed reads
    input_R1R2 = csv_channel.map { row -> tuple(row.Sample_name, file(row.R1), file(row.R2)) }

    // Run Process: Fastqc
    input_R1R2 | PostTrimFastqc

    //ine Input: Fastqc htmls. Using the output of Fastqc as an input for MultiQC
    fastqc_output = PostTrimFastqc.out.fastqc_zips


    //Run Process: MultQC
    fastqc_output | MultiQC

    // ine Input: Bowtie 
    input_bowtie = csv_channel.map { row -> tuple(row.Sample_name, file(row.Transcriptome1), file(row.R1), file(row.R2)) }

    //Run Process: Bowtie
    input_bowtie | Bowtie

    // ine Input: Bowtie2 
    input_bowtie2 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome2), file(row.R1), file(row.R2)) }

    //Run Process: Bowtie2
    input_bowtie2 | Bowtie2

    //ine Input: Trinity Fasta 
    input_trinity_fasta = csv_channel.map { row -> tuple(row.Sample_name, file(row.Transcriptome1)) }

    //Run Process: TrinityStats
    input_trinity_fasta | TrinityStats


    //ine Input: Trinity Fasta 
    input_trinity_fasta2 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome2)) }

    //Run Process: TrinityStats
    input_trinity_fasta2 | TrinityStats2

    // Transcriptome 1 
    //ine Input: Transcriptome1 Fasta + BUSCOlin1 tuple 
    input_BUSCOlin1 = csv_channel.map { row -> tuple(row.Sample_name, file(row.Transcriptome1), row.BUSCO_lin1, row.Transcriptome1_label) }

    //Run Process: BUSCO_lin1
    input_BUSCOlin1 | BUSCO_transcriptome_metazoa

    //ine Input: Transcriptome1 Fasta + BUSCOlin2 tuple 
    input_BUSCOlin2 = csv_channel.map { row -> tuple(row.Sample_name, file(row.Transcriptome1), row.BUSCO_lin2, row.Transcriptome1_label) }

    //Run Process: BUSCO_lin2
    input_BUSCOlin2 | BUSCO_transcriptome_mollusca

    // Transcriptome 2

    //ine Input: Transcriptome2 Fasta + BUSCOlin1 tuple 
    input_BUSCOlin1_2 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome2), row.BUSCO_lin1, row.Transcriptome2_label) }

    //Run Process: BUSCO_lin1
    input_BUSCOlin1_2 | BUSCO_transcriptome_metazoa2

    //ine Input: Transcriptome2 Fasta + BUSCOlin2 tuple 
    input_BUSCOlin2_2 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome2), row.BUSCO_lin2, row.Transcriptome2_label) }

    //Run Process: BUSCO_lin2
    input_BUSCOlin2_2 | BUSCO_transcriptome_mollusca2

    //ine Input:  Transcriptome_Combined 
    input_Transcriptome_Combined = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, file(row.Transcriptome1), row.Transcriptome1_label, file(row.Transcriptome2), row.Transcriptome2_label) }
    // Run process ; 
    input_Transcriptome_Combined | Transcriptome_Combined


    // Transcriptome_Combined
    //ine Input: Combined transcriptome Fasta + BUSCOlin1 tuple 
    input_BUSCOlin1_3 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, row.BUSCO_lin1) }
        .join(Transcriptome_Combined.out.transcriptome_combined)

    //Run Process: BUSCO_lin1
    input_BUSCOlin1_3 | BUSCO_transcriptome_metazoa3

    //ine Input: Combined transcriptome Fasta + BUSCOlin2 tuple 
    input_BUSCOlin2_3 = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row -> tuple(row.Sample_name, row.BUSCO_lin2) }
        .join(Transcriptome_Combined.out.transcriptome_combined)

    //Run Process: BUSCO_lin2
    input_BUSCOlin2_3 | BUSCO_transcriptome_mollusca3






    // ine Input: Trinity fasta + R1 + R2 tuple 
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


    //ine Input: Database_Fasta. This is set up to allow for multiple different databases if needed.
    // ToxinFasta names entry list
    Toxin_fasta_file = file(params.toxprot_fasta, checkIfExists: false).exists()
        ? file(params.toxprot_fasta)
        : file(params.Fallback_toxin_fasta)
    input_databasefasta = Channel.fromPath(Toxin_fasta_file)


    //Run Process: DatabaseCreation
    input_databasefasta | Blastdatabasecreation


    //ine Input: Blastx 
    Blastxinputfasta_single = csv_channel
        .filter { row -> !row.Transcriptome2?.trim() || row.Transcriptome2.trim().toLowerCase() == 'null' }
        .map { row ->
            tuple(row.Sample_name, file(row.Transcriptome1))
        }


    Blastxinputfasta_combined = csv_channel
        .filter { row -> row.Transcriptome2?.trim() && row.Transcriptome2.trim() != '' && row.Transcriptome2.trim().toLowerCase() != 'null' }
        .map { row ->
            tuple(row.Sample_name)
        }
        .join(Transcriptome_Combined.out.transcriptome_combined)

    // Combine both channels using mix() operator
    Blastxinputfasta_all_1 = Blastxinputfasta_single.mix(Blastxinputfasta_combined)
    Blastxinputfasta_all = Blastxinputfasta_all_1.combine(Blastdatabasecreation.out.proteindb)
    // Run Process: Blastx 
    Blastxinputfasta_all | Blastx
    //Run Process: Transdecoder
    Transcriptpome1 = csv_channel
        .filter { row -> !row.Transcriptome2?.trim() || row.Transcriptome2.trim().toLowerCase() == 'null' }
        .map { row ->
            tuple(row.Sample_name, file(row.Transcriptome1))
        }

    input_orf = Transcriptpome1.mix(Transcriptome_Combined.out.transcriptome_combined)


    // Initialize variables that will be used after the conditional blocks
    BUSCOlin1 = csv_channel.map { row -> tuple(row.Sample_name, row.BUSCO_lin1) }
    BUSCOlin2 = csv_channel.map { row -> tuple(row.Sample_name, row.BUSCO_lin2) }

    // Transdecoder
    if (params.ORFPrediction == "TD") {
        // TD only
        input_orf | Transdecoder
        //ine Input: Transdecoder pep + BUSCOlin1 tuple 
        Transdecoder_pep = Transdecoder.out.transdecoder_pep
        input_BUSCOlin1_L = Transdecoder_pep.join(BUSCOlin1)

        //Run Process: BUSCO_lin1
        input_BUSCOlin1_L | BUSCO_translatome_metazoa

        //ine Input: Transdecoder pep + BUSCOlin2 tuple 
        input_BUSCOlin2_L = Transdecoder_pep.join(BUSCOlin2)

        //Run Process: BUSCO_lin2
        input_BUSCOlin2_L | BUSCO_translatome_mollusca
        //ine Input for input_ORF_complete 
        Transdecodercds = Transdecoder.out.transdecoder_cds
        Transdecoderpep = Transdecoder.out.transdecoder_pep

        input_ORFs_Combined_CDHit_Input = Transdecoderpep.join(Transdecodercds)
        input_ORFs_Combined_CDHit_Input | ORFs_Combined_CDHit
        input_ORF_complete = ORFs_Combined_CDHit.out.combined_pep.join(ORFs_Combined_CDHit.out.combined_cds)
        input_Blastp = ORFs_Combined_CDHit.out.combined_pep.combine(Blastdatabasecreation.out.proteindb)

    }
    else if (params.ORFPrediction == "Both") {
        // Both 
        input_orf | Transdecoder
        input_orf | TD2
        // ine input: ORFs_Combined
        input_ORFs_Combined = Transdecoder.out.transdecoder_pep.join(Transdecoder.out.transdecoder_cds).join(TD2.out.TD2_pep).join(TD2.out.TD2_cds)
        // Run Process: ORFs_Combined
        input_ORFs_Combined | ORFs_Combined

        //ine Input: Transdecoder pep + BUSCOlin1 tuple 
        Transdecoder_pep = Transdecoder.out.transdecoder_pep
        input_BUSCOlin1_L = Transdecoder_pep.join(BUSCOlin1)

        //Run Process: BUSCO_lin1
        input_BUSCOlin1_L | BUSCO_translatome_metazoa

        //ine Input: Transdecoder pep + BUSCOlin2 tuple 
        input_BUSCOlin2_L = Transdecoder_pep.join(BUSCOlin2)

        //Run Process: BUSCO_lin2
        input_BUSCOlin2_L | BUSCO_translatome_mollusca

        // TD2
        //ine Input: Transdecoder pep + BUSCOlin1 tuple 
        TD2_pep = TD2.out.TD2_pep
        input_BUSCOlin1_L_2 = TD2_pep.join(BUSCOlin1)
        //Run Process: BUSCO_lin1
        input_BUSCOlin1_L_2 | BUSCO_translatome_metazoa2
        //ine Input: Transdecoder pep + BUSCOlin2 tuple 
        input_BUSCOlin2_L_2 = TD2_pep.join(BUSCOlin2)
        //Run Process: BUSCO_lin2
        input_BUSCOlin2_L_2 | BUSCO_translatome_mollusca2
        // Combined
        //ine Input: Transdecoder pep + BUSCOlin1 tuple 
        Combined_pep = ORFs_Combined.out.combined_pep
        input_BUSCOlin1_L_3 = Combined_pep.join(BUSCOlin1)
        //Run Process: BUSCO_lin1
        input_BUSCOlin1_L_3 | BUSCO_translatome_metazoa3
        //ine Input: Transdecoder pep + BUSCOlin2 tuple 
        input_BUSCOlin2_L_3 = Combined_pep.join(BUSCOlin2)
        //Run Process: BUSCO_lin2
        input_BUSCOlin2_L_3 | BUSCO_translatome_mollusca3

        //ine Input for input_ORF_complete 
        input_ORFs_Combined_CDHit_Input = ORFs_Combined.out.combined_pep.join(ORFs_Combined.out.combined_cds)
        input_ORFs_Combined_CDHit_Input | ORFs_Combined_CDHit

        input_ORF_complete = ORFs_Combined_CDHit.out.combined_pep.join(ORFs_Combined_CDHit.out.combined_cds)
        input_Blastp = ORFs_Combined_CDHit.out.combined_pep.combine(Blastdatabasecreation.out.proteindb)
    }
    else if (params.ORFPrediction == "Done"){
        ORFspep = csv_channel.map { row -> tuple(row.Sample_name, file(row.ORFpep)) }
        ORFscds = csv_channel.map { row -> tuple(row.Sample_name, file(row.ORFcds)) }
        input_ORFs_Combined_CDHit_Input = ORFspep.join(ORFscds)
        input_ORFs_Combined_CDHit_Input | ORFs_Combined_CDHit
        input_ORF_complete = ORFs_Combined_CDHit.out.combined_pep.join(ORFs_Combined_CDHit.out.combined_cds)
        input_Blastp = ORFs_Combined_CDHit.out.combined_pep.combine(Blastdatabasecreation.out.proteindb)
    }
    else {
        // TD2  logic 
        input_orf | TD2
        // TD2
        //ine Input: Transdecoder pep + BUSCOlin1 tuple 
        TD2_pep = TD2.out.TD2_pep
        input_BUSCOlin1_L_2 = TD2_pep.join(BUSCOlin1)
        //Run Process: BUSCO_lin1
        input_BUSCOlin1_L_2 | BUSCO_translatome_metazoa2
        //ine Input: Transdecoder pep + BUSCOlin2 tuple 
        input_BUSCOlin2_L_2 = TD2_pep.join(BUSCOlin2)
        //Run Process: BUSCO_lin2
        input_BUSCOlin2_L_2 | BUSCO_translatome_mollusca2

        //ine Input for input_ORF_complete 
        TD2cds = TD2.out.TD2_cds
        TD2pep = TD2.out.TD2_pep
        input_ORFs_Combined_CDHit_Input = TD2pep.join(TD2cds)
        input_ORFs_Combined_CDHit_Input | ORFs_Combined_CDHit
        input_ORF_complete = ORFs_Combined_CDHit.out.combined_pep.join(ORFs_Combined_CDHit.out.combined_cds)
        input_Blastp = ORFs_Combined_CDHit.out.combined_pep.combine(Blastdatabasecreation.out.proteindb)
    }

    //Run Process: Transdecoder filter for complete ORFs
    input_ORF_complete | ORF_complete


    //ine Input: Transdecoder cds + R1 + R2 + Strandedness tuple 
    Complete_cds = ORF_complete.out.complete_cds
    KallistoTransdecoderR1R2S = csv_channel.map { row -> tuple(row.Sample_name, file(row.R1), file(row.R2), row.Strandedness) }
    input_TransKallisto = Complete_cds.join(KallistoTransdecoderR1R2S)

    //Run Process: TransKallisto
    input_TransKallisto | Kallisto_Transdecoder


    //Run Process: Blastp
    input_Blastp | Blastp


    //ine Input: Sample name + completepep tuple
    input_signalp = ORF_complete.out.complete_pep

    //Run Process: SignalP
    input_signalp | SignalP

    //ine Input: sample name + mature sequences + completepep tuple 
    maturesequences = SignalP.out.maturesequences
    maturecomplete = maturesequences.join(input_signalp).join(ORF_complete.out.complete_cds)

    //Run Process: Filter2   
    maturecomplete | Filter2

    //ine Input: stats sample name + pep + cds + completepep + completecds + mature +signalp  tuple 
    // With the optional NULL's the order here does not match the input tuple. but for this particular proceess it doesnt matter because it is just running seqkit stats on all of them
    if (params.ORFPrediction == "TD") {
        stats_join = Transdecoder_pep
            .join(Transdecoder.out.transdecoder_cds)
            .join(ORF_complete.out.complete_pep)
            .join(ORF_complete.out.complete_cds)
            .join(maturesequences)
            .join(Filter2.out.complete_pep_signalp)
            .join(input_ORF_complete)
            .combine(channel.of("NULL"))
            .combine(channel.of("NULL"))
    }
    else if (params.ORFPrediction == "Both") {
        stats_join = Transdecoder_pep
            .join(Transdecoder.out.transdecoder_cds)
            .join(ORF_complete.out.complete_pep)
            .join(ORF_complete.out.complete_cds)
            .join(maturesequences)
            .join(Filter2.out.complete_pep_signalp)
            .join(TD2_pep)
            .join(TD2.out.TD2_cds)
            .join(input_ORF_complete)
    }
    else if (params.ORFPrediction == "Done") {
        stats_join = ORF_complete.out.complete_pep.join(ORF_complete.out.complete_cds).join(maturesequences).join(Filter2.out.complete_pep_signalp).join(ORFspep).join(ORFscds).join(input_ORF_complete).combine(channel.of("NULL")).combine(channel.of("NULL"))
    }
    else {
        stats_join = ORF_complete.out.complete_pep.join(ORF_complete.out.complete_cds).join(maturesequences).join(Filter2.out.complete_pep_signalp).join(TD2_pep).join(TD2.out.TD2_cds).join(input_ORF_complete).combine(channel.of("NULL")).combine(channel.of("NULL"))
    }



    //Run Process: STATS   
    stats_join | stats

    //ine Input: genomefasta 
    Genomefasta = csv_channel
        .filter { row ->
            def path = row.Genome_fasta_path
            path != null && path != 'NULL' && path.toString() != 'NULL' && path.toString().trim() != ''
        }
        .map { row -> tuple(row.Sample_name, file(row.Genome_fasta_path)) }
        .unique { tuple -> tuple[0] }




    // Run Process: BlastnGenome database creation
    GenomeBlastdatabasecreation(Genomefasta)
    genomedb = GenomeBlastdatabasecreation.out.genomedbfiles


    // Optional Params 
    if (params.DeepTMHMM) {
        signalpsummary = SignalP.out.signalpsummary
        DeepTMHMMinput = input_signalp.join(signalpsummary)
        DeepTMHMMinput | DeepTMHMM
        input_DeepTMHMMFilter = DeepTMHMM.out.mature
            .join(ORF_complete.out.complete_pep)
            .join(ORF_complete.out.complete_cds)
            .join(Filter2.out.complete_pep_signalp)
            .join(Filter2.out.complete_cds_signalp)
            .join(SignalP.out.maturesequences)
        input_DeepTMHMMFilter | DeepTMHMMFilter
        input_Interproscan = DeepTMHMMFilter.out.complete_pep_secreted
        BlastnInput = DeepTMHMMFilter.out.complete_cds_secreted.join(genomedb)
    }
    else {
        input_Interproscan = Filter2.out.complete_pep_signalp
        //ine Input: Genome BLAST 
        BlastnInput = Filter2.out.complete_cds_signalp.join(genomedb)
    }

    //Run Process: Interproscan  (Complete secreted pep ORFs only)
    input_Interproscan | Interproscan


    // ine Input: BlastdatabasecreationNonToxin
    // Panther names entry list
    NonToxin_fasta_file = file(params.nontoxprot_fasta, checkIfExists: false).exists()
        ? file(params.nontoxprot_fasta)
        : file(params.Fallback_nontoxin_fasta)
    input_nontoxindatabasefasta = Channel.fromPath(NonToxin_fasta_file)

    //Run Process: BlastdatabasecreationNonToxin
    input_nontoxindatabasefasta | BlastdatabasecreationNonToxin

    // ine Input: BlastpNonToxin
    input_nontoxinBlastp = input_Interproscan.combine(BlastdatabasecreationNonToxin.out.nontoxinproteindb)

    //Run Process:BlastpNonToxin
    input_nontoxinBlastp | BlastpNonToxin

    //Run Process: BlastnGenome
    BlastnInput | GenomeBlasts6
    //Run Process: BlastnGenome
    BlastnInput | GenomeBlasts0
    


}
