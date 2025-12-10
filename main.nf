#!/usr/bin/env nextflow



// Process 1: PostTrimFastqc
process PostTrimFastqc {

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Fastqc/posttrim/fastqc_zips/", mode: 'copy'

    input:
    tuple path(R1), path(R2)

    output:
    path '*.zip', emit: fastqc_zips

    script:
    """
    fastqc ${R1} ${R2}
    """
}

// Process 2: MultiQC
process MultiQC {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/Fastqc/posttrim/", mode: 'copy'

    input:
    tuple val(sample), path(fastqc_zips)

    output:
    path 'multiqc_report.html', emit: multiqc

    script:

    """
    multiqc ${fastqc_zips}
    """
}

// Process 3: TrinityStats
process TrinityStats {


    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/Stats/", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta)

    output:
    path "*.txt", emit: trinity_stats

    script:

    """
    seqkit stats ${trinity_fasta} > ${sample}_Trinity.stats.txt

    """
}
// Process 4: BUSCO_transcriptome_metazoa
process BUSCO_transcriptome_metazoa {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco.yaml"

    publishDir "${sample}/results/BUSCO/transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), path(metazoa)

    output:
    path "*.txt", emit: busco_transcriptome_met

    script:

    """
    busco -i ${trinity_fasta} -l ${metazoa} -c 10 -o ${sample}_met.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_met.transcriptome/*.txt "."
    """
}

// Process 5: BUSCO_transcriptome_mollusca
process BUSCO_transcriptome_mollusca {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco.yaml"

    publishDir "${sample}/results/BUSCO/transcriptome/", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), path(mollusca)

    output:
    path "*.txt", emit: busco_transcriptome_mollusca

    script:

    """
    busco -i ${trinity_fasta} -l ${mollusca} -c 10 -o ${sample}_mol.transcriptome -m transcriptome -e 1e-5 -f
    mv ${sample}_mol.transcriptome/*.txt "."
    """
}



// Process 6: Kallisto_Trinity
process Kallisto_Trinity {

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/kallisto/trinity/output", mode: 'copy'

    errorStrategy 'retry'
    maxRetries 2
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

// Process 7: Blastdatabasecreation
process Blastdatabasecreation {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    input:
    path database_fasta

    output:
    path "*", emit: proteindb

    script:
    """

    makeblastdb -in "${database_fasta}" -dbtype prot

    """
}

// Process 8: Blastx
process Blastx {

    errorStrategy 'retry'
    maxRetries 2
    cpus { task.attempt * 2 }

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/Blast/Blastx/", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta), path(proteindb), val(proteindbdbname)

    output:
    path "${sample}.blastx.db.0.txt", emit: blastx0
    path "${sample}.blastx.db.6.txt", emit: blastx6

    script:
    """
   
    blastx -query ${trinity_fasta} -db ${proteindbdbname} -out ${sample}.blastx.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch evalue bitscore qframe"
    blastx -query ${trinity_fasta} -db ${proteindbdbname} -out ${sample}.blastx.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """
}

// Process 9: Transdecoder
process Transdecoder {

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    errorStrategy 'ignore'

    publishDir "${sample}/results/Transdecoder", mode: 'copy'

    input:
    tuple val(sample), path(trinity_fasta)

    output:
    path "*.pep", emit: pep
    path "*.cds", emit: cds

    script:

    """
    TransDecoder.LongOrfs -t ${trinity_fasta}
    TransDecoder.Predict -t ${trinity_fasta}

    """
}

// Process 10: BUSCO_translatome_metazoa
process BUSCO_translatome_metazoa {


    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco.yaml"

    publishDir "${sample}/results/BUSCO/translatome/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), path(metazoa)

    output:
    path "*.txt", emit: busco_translatome_met

    script:

    """
    busco -i ${Transdecoder_pep} -l ${metazoa} -c 10 -o ${sample}_met.protein -m protein -e 1e-5 -f
    mv ${sample}_met.protein/*.txt "."
    
    """
}

// Process 11: BUSCO_translatome_mollusca
process BUSCO_translatome_mollusca {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco.yaml"

    publishDir "${sample}/results/BUSCO/translatome/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), path(mollusca)

    output:
    path "*.txt", emit: busco_translatome_met

    script:

    """
    busco -i ${Transdecoder_pep} -l ${mollusca} -c 10 -o ${sample}_mol.protein -m protein -e 1e-5 -f
    mv ${sample}_mol.protein/*.txt "."
    
    """
}


// Process 12: Kallisto_Transdecoder
process Kallisto_Transdecoder {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/kallisto/transdecoder/output", mode: 'copy'

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

// Process 13: Blastp
process Blastp {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/Blast/Blastp/", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), path(proteindb), val(proteindbdbname)

    output:
    path "${sample}.blastp.db.0.txt", emit: blastp0
    path "${sample}.blastp.db.6.txt", emit: blastp6

    script:
    """

    blastp -query ${Transdecoder_pep} -db ${proteindbdbname} -out ${sample}.blastp.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch evalue bitscore qframe"
    blastp -query ${Transdecoder_pep} -db ${proteindbdbname} -out ${sample}.blastp.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """
}

// Process 14: Transdecoder_complete
process Transdecoder_complete {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/Transdecoder", mode: 'copy'

    input:
    tuple val(sample), path(Transdecoder_pep), path(Transdecoder_cds)

    output:
    path "*.pep", emit: transdecodercomplete_pep
    path "*.cds", emit: transdecodercomplete_cds

    script:

    """
    seqkit grep -n -r -p "ORF type:complete" ${Transdecoder_pep} -o "${sample}.trandescoder.complete.pep"
    seqkit grep -n -r -p "ORF type:complete" ${Transdecoder_cds} -o "${sample}.trandescoder.complete.cds"
    """
}

// Process 15: SignalP
process SignalP {

    errorStrategy 'ignore'

    publishDir "${sample}/results/Signalp", mode: 'copy'

    input:
    tuple val(sample), path(transdecodercomplete_pep)

    output:
    path "*_mature.fasta", emit: maturesequences
    path "*.signalp5", emit: signalpsummary

    script:


    """
    signalp -fasta ${transdecodercomplete_pep} -mature -prefix "${sample}"
    """
}

// Process 16: Filter2
process Filter2 {
    errorStrategy 'ignore'
    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/Transdecoder", mode: 'copy'

    input:
    tuple val(sample), path(maturesequences), path(transdecodercomplete_pep)

    output:
    path "*.fasta", emit: transdecoderpep_signalp

    script:

    """
    seqkit grep -f <(seqkit seq -n ${maturesequences}) ${transdecodercomplete_pep} > ${sample}.Transdecoder.complete.signalp.sequences.fasta
    """
}

// Process 17: STATS
process stats {
    errorStrategy 'ignore'
    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/Stats", mode: 'copy'

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

// Process 18: Interproscan
process Interproscan {

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/Interproscan", mode: 'copy'

    errorStrategy 'retry'
    maxRetries 2
    cpus { task.attempt * 2 }
    memory { task.attempt * 2.GB }

    input:
    tuple val(sample), path(Transdecoder_pep)

    output:
    path "*.tsv", emit: Interproscan

    script:

    """
    awk '{if (\$0 ~ /^>/) print \$0; else {gsub(/\\*/, ""); print \$0}}' ${Transdecoder_pep} > "${sample}.Trinity.fasta.transdecoder.cleaned.pep"

    interproscan.sh
      -goterms -i "${sample}.Trinity.fasta.transdecoder.cleaned.pep" 
      -pa -t p -d ./ -f TSV
    """
}

// Process 19: GenomeBlastdatabasecreation
process GenomeBlastdatabasecreation {
    errorStrategy 'ignore'
    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    input:
    path genome_fasta

    output:
    path "*", emit: genomedb

    script:
    """

    makeblastdb -in "${genome_fasta}" -dbtype nucl

    """
}

// Process 20: GenomeBlasts
process GenomeBlasts {
    errorStrategy 'ignore'
    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "${sample}/results/Blast/Blastn/", mode: 'copy'

    input:
    tuple val(sample), path(transdecodercompletecds), path(genomedb), val(genomedbname)

    output:
    path "${sample}.blastn.db.0.txt", emit: blastn0
    path "${sample}.blastn.db.6.txt", emit: blastn6

    script:
    """
    blastn -query ${transdecodercompletecds} -db ${genomedbname} -out ${sample}.blastn.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt '6'
    blastn -query ${transdecodercompletecds} -db ${genomedbname} -out ${sample}.blastn.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """
}







// PARAMETERS DEFINED IN CONFIG FILE

//WorkFlow

workflow {
    // Define CSV channel 
    csv_channel = Channel.fromPath(params.input_csv).splitCsv(header: true, sep: ',')

    //Define Sample_name 
    Sample_name = csv_channel.map { row -> row.Sample_name }

    // Define Input: Paired Trimmed reads Tuple 
    input_R1R2 = csv_channel.map { row -> tuple(file(row.R1), file(row.R2)) }

    // Run Process: Fastqc
    input_R1R2 | PostTrimFastqc

    //Define Input: Fastqc htmls
    def fastqc_output = PostTrimFastqc.out.fastqc_zips
    def fastqc_zips = Sample_name.join(fastqc_output)
    //Run Process: MultQC
    fastqc_zips | MultiQC

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

    //Define Input: Database_Fasta 
    input_databasefasta = csv_channel.map { row -> file(row.Protein_fasta_path_for_Blast) }

    //Run Process: DatabaseCreation
    input_databasefasta | Blastdatabasecreation

    //Define Input: Blastx //FILE OR PATH FOR DATABASE FILES
    input_Blastx = input_trinity_fasta
        .combine(Blastdatabasecreation.out.proteindb)
        .combine(csv_channel.map { row -> tuple(row.Sample_name, row.Protein_fasta_name) })
        .filter { it[0] == it[3] }
        .map { sample, trinity_fasta, proteindb, _sample2, dbname -> tuple(sample, trinity_fasta, proteindb, dbname) }

    //Run Process: Blastx
    input_Blastx | Blastx

    //Run Process: Transdecoder
    input_trinity_fasta | Transdecoder

    //Define Input: Transdecoder pep + BUSCOlin1 tuple 
    def Transdecoder_pep = Transdecoder.out.pep
    def Sample_transdecoder = Sample_name.join(Transdecoder_pep)
    
    input_BUSCOlin1_L = Sample_transdecoder
        .combine(csv_channel.map { row -> tuple(row.Sample_name, file(row.BUSCO_lin1)) })
        .filter { it[0] == it[2] }
        .map { sample, pep, _sample2, busco_lin1 -> tuple(sample, pep, busco_lin1) }

    //Run Process: BUSCO_lin1
    input_BUSCOlin1_L | BUSCO_translatome_metazoa

    //Define Input: Transdecoder pep + BUSCOlin2 tuple 
    input_BUSCOlin2_L = Sample_transdecoder
        .combine(csv_channel.map { row -> tuple(row.Sample_name, file(row.BUSCO_lin2)) })
        .filter { it[0] == it[2] }
        .map { sample, pep, _sample2, busco_lin2 -> tuple(sample, pep, busco_lin2) }

    //Run Process: BUSCO_lin2
    input_BUSCOlin2_L | BUSCO_translatome_mollusca

    //Define Input: Transdecoder cds + R1 + R2 + Strandedness tuple 
    def Transdecoder_cds = Transdecoder.out.cds
    def Sample_transdecoder_cds = Sample_name.join(Transdecoder_cds)
    
    input_TransKallisto = Sample_transdecoder_cds
        .combine(csv_channel.map { row -> tuple(row.Sample_name, file(row.R1), file(row.R2), row.Strandedness) })
        .filter { it[0] == it[2] }
        .map { sample, cds, _sample2, r1, r2, strand -> tuple(sample, cds, r1, r2, strand) }

    //Run Process: TransKallisto
    input_TransKallisto | Kallisto_Transdecoder

    //Define Input: Blastp //FILE OR PATH FOR DATABASE FILES
    input_Blastp = Sample_transdecoder
        .combine(Blastdatabasecreation.out.proteindb)
        .combine(csv_channel.map { row -> tuple(row.Sample_name, row.Protein_fasta_name) })
        .filter { it[0] == it[3] }
        .map { sample, pep, proteindb, _sample2, dbname -> tuple(sample, pep, proteindb, dbname) }

    //Run Process: Blastp
    input_Blastp | Blastp

    //Define Input: Sample name + PEP + CDS tuple
    def Transdecoder_complete = Sample_name.join(Transdecoder_pep).join(Transdecoder_cds)
    //Run Process: Transdecoder filter for complete ORFs
    Transdecoder_complete | Transdecoder_complete


    //Define Input: Sample name + completepep tuple
    def transdecodercomplete_pep = Transdecoder_complete.out.transdecodercomplete_pep
    def transdecodercomplete_pep_sampleid = Sample_name.join(transdecodercomplete_pep)

    //Run Process: SignalP
    transdecodercomplete_pep_sampleid | SignalP

    //Define Input: sample name + mature sequences + completepep tuple 
    def maturesequences = SignalP.out.maturesequences
    maturecomplete = Sample_name.join(maturesequences).join(transdecodercomplete_pep)

    //Run Process: Filter2   
    maturecomplete | Filter2

    //Define Input: stats sample name + pep + cds + completepep + completecds + mature +signalp  tuple 
    def transdecoderpep_signalp = Filter2.out.transdecoderpep_signalp
    def transdecodercomplete_cds = Transdecoder_complete.out.transdecodercomplete_cds
    stats_join = Sample_name
        .join(Transdecoder_pep)
        .join(Transdecoder_cds)
        .join(transdecodercomplete_pep)
        .join(transdecodercomplete_cds)
        .join(maturesequences)
        .join(transdecoderpep_signalp)
    //Run Process: STATS   
    stats_join | stats

    //Define Input: sample name + pep tuple 
    inputInterpro = Sample_name.join(Transdecoder_pep)

    //Run Process: Interproscan   
    inputInterpro | Interproscan

    //Define Input: genomefasta 
    Genomefasta = csv_channel
        .map { row -> file(row.Genome_fasta_path) }
        .unique()
        
    BLASTntuple = Sample_name
        .join(transdecodercomplete_cds)
        .combine(GenomeBlastdatabasecreation.out.genomedb)
        .combine(csv_channel.map { row -> tuple(row.Sample_name, row.Genome_fasta_name) })
        .filter { it[0] == it[3] }
        .map { sample, cds, genomedb, _sample2, genomename -> tuple(sample, cds, genomedb, genomename) }
    
    //Run Process: BlastnGenome   
    Genomefasta
        .filter { it.name != 'NULL' && it.toString() != 'NULL' }
        | GenomeBlastdatabasecreation
    
    BLASTntuple | GenomeBlasts
}
