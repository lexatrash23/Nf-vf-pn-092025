#!/usr/bin/env nextflow

def launchDir = System.getProperty('user.dir')



if ( params.help ) {
    help = """main.nf: Nextflow pipeline that takes as inputs paired trimmed reads and a trinity assembly. Processes inputs through Fastqc, Transdecoder, SignalP, Blast and other analysis to identify putative venom transcripts. For more information please refer to [Insert Final Github link] for README file
             |Required arguments:
             |  --trinity_fasta  Path to trinity fasta file
                --R1 Path to Trimmed Reads 1
                --R2 Path to Trimmed Reads 2
                --TRANSDECODER_PATH Path to Transdecoder directory containing the .longorfs and .predict scripts
                --SIGNALP_PATH Path to SignalP script
                --INTERPROSCAN_PATH Path to Interproscan
                --Sample_name Sample name
                --metazoa Path to BUSCO lineages files 1
                --mollusca Path to BUSCO lineages files 2
                --database_fasta Path to protein fasta file to be used for BLAST
                --database_name Database name (with fasta suffix)
                --stranded_input Library strand orientation (can be 'fr', 'rf', or '')
                --genomefasta Genome fasta (NULL if not available)
                --genomefastaname Genome fasta database name (with fasta suffix, NULL if not available)""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}


// Process 1: PostTrimFastqc
process PostTrimFastqc {

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

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

    publishDir "results/Fastqc/posttrim/", mode: 'copy'

    input:
    path fastqc_zips

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

    publishDir "results/Stats/", mode: 'copy'

    input:
    path trinity_fasta

    output:
    path "*.txt", emit: trinity_stats


    script:

    """
    seqkit stats ${trinity_fasta} > ${params.Sample_name}_Trinity.stats.txt

    """

}
// Process 4: BUSCO_transcriptome_metazoa
process BUSCO_transcriptome_metazoa {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco.yaml"

    publishDir "results/BUSCO/transcriptome/", mode: 'copy'

    input:
    tuple path(trinity_fasta), path(metazoa)

    output:
    path "*.txt", emit: busco_transcriptome_met


    script:

    """
    busco -i ${trinity_fasta} -l ${metazoa} -c 10 -o ${params.Sample_name}_met.transcriptome -m transcriptome -e 1e-5 -f
    mv ${params.Sample_name}_met.transcriptome/*.txt "."
    """

}

// Process 5: BUSCO_transcriptome_mollusca
process BUSCO_transcriptome_mollusca {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco.yaml"

    publishDir "results/BUSCO/transcriptome/", mode: 'copy'

    input:
    tuple path(trinity_fasta), path(mollusca)

    output:
    path "*.txt", emit: busco_transcriptome_mollusca


    script:

    """
    busco -i ${trinity_fasta} -l ${mollusca} -c 10 -o ${params.Sample_name}_mol.transcriptome -m transcriptome -e 1e-5 -f
    mv ${params.Sample_name}_mol.transcriptome/*.txt "."
    """

}



// Process 7: Kallisto_Trinity
process Kallisto_Trinity {

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/kallisto/trinity/output", mode: 'copy'
    
    errorStrategy 'retry'
    maxRetries 2
    cpus { task.attempt * 2 }
    memory { (task.attempt * 2) * 1.9.GB }

    input:
    tuple path(trinity_fasta), path(R1), path(R2)

    output:
    path "abundance.tsv", emit: KallistoTrinityAbundance


    script:
    """
    stranded_input="${params.stranded_input}"

    if [ "\$stranded_input" == "fr" ]; then
        stranded_flag="--fr-stranded"
    elif [ "\$stranded_input" == "rf" ]; then
        stranded_flag="--rf-stranded"
    else
        stranded_flag=""
    fi
    kallisto index -i index ${trinity_fasta}
    kallisto quant -i index -o ./ -b 100 ${R1} ${R2} \$stranded_flag

    """

}

// Process 8a: Blastdatabasecreation
process Blastdatabasecreation {
    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    input:
    path(database_fasta)

    output:
    path "*", emit: proteindb

    script:
    """

    makeblastdb -in "${database_fasta}" -dbtype prot

    """

}

// Process 8: Blastx
process Blastx {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Blast/Blastx/", mode: 'copy'

    input:
    path(trinity_fasta)
    path(proteindb)

    output:
    path "${params.Sample_name}.blastx.db.0.txt", emit: blastx0
    path "${params.Sample_name}.blastx.db.6.txt", emit: blastx6



    script:
    """
   
    blastx -query ${trinity_fasta} -db ${params.database_name} -out ${params.Sample_name}.blastx.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch evalue bitscore qframe"
    blastx -query ${trinity_fasta} -db ${params.database_name} -out ${params.Sample_name}.blastx.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """

}

// Process 9: Transdecoder
process Transdecoder {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Transdecoder", mode: 'copy'

    input:
    path trinity_fasta

    output:
    path "*.pep", emit: pep
    path "*.cds", emit: cds


    script:

    """
    ${params.TRANSDECODER_PATH}/TransDecoder.LongOrfs -t ${trinity_fasta}
    ${params.TRANSDECODER_PATH}/TransDecoder.Predict -t ${trinity_fasta}

    """

}

// Process 10: BUSCO_translatome_metazoa
process BUSCO_translatome_metazoa {


    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco.yaml"

    publishDir "results/BUSCO/translatome/", mode: 'copy'


    input:
    tuple path(Transdecoder_pep), path(metazoa)

    output:
    path "*.txt", emit: busco_translatome_met


    script:

    """
    busco -i ${Transdecoder_pep} -l ${metazoa} -c 10 -o ${params.Sample_name}_met.protein -m protein -e 1e-5 -f
    mv ${params.Sample_name}_met.protein/*.txt "."
    
    """

}

// Process 11: BUSCO_translatome_mollusca
process BUSCO_translatome_mollusca {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/busco.yaml"

    publishDir "results/BUSCO/translatome/", mode: 'copy'


    input:
    tuple path(Transdecoder_pep), path(mollusca)

    output:
    path "*.txt", emit: busco_translatome_met


    script:

    """
    busco -i ${Transdecoder_pep} -l ${mollusca} -c 10 -o ${params.Sample_name}_mol.protein -m protein -e 1e-5 -f
    mv ${params.Sample_name}_mol.protein/*.txt "."
    
    """

}


// Process 12: Kallisto_Transdecoder
process Kallisto_Transdecoder {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/kallisto/transdecoder/output", mode: 'copy'

    input:
    tuple path(Transdecoder_cds), path(R1), path(R2)

    output:
    path "abundance.tsv", emit: KallistoTransdecoderAbundance


    script:
    if (params.stranded_input == 'fr') {
        stranded_flag = '--fr-stranded'
    } else if (params.stranded_input == 'rf') {
        stranded_flag = '--rf-stranded'
    } else {
        stranded_flag = ''
    }
    """
    kallisto index -i index ${Transdecoder_cds}
    kallisto quant -i index -o ./ -b 100 ${stranded_flag} ${R1} ${R2}

    """

}

// Process 13: Blastp
process Blastp {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Blast/Blastp/", mode: 'copy'

    input:
    path(Transdecoder_pep)
    path(proteindb)

    output:
    path "${params.Sample_name}.blastp.db.0.txt", emit: blastp0
    path "${params.Sample_name}.blastp.db.6.txt", emit: blastp6


    script:
    """

    blastp -query ${Transdecoder_pep} -db ${params.database_name} -out ${params.Sample_name}.blastp.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch evalue bitscore qframe"
    blastp -query ${Transdecoder_pep} -db ${params.database_name} -out ${params.Sample_name}.blastp.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """

}

// Process 14: Transdecoder_complete
process Transdecoder_complete {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Transdecoder", mode: 'copy'

    input:
    tuple path(Transdecoder_pep), path(Transdecoder_cds)

    output:
    path "*.pep", emit: transdecodercomplete_pep
    path "*.cds", emit: transdecodercomplete_cds

    script:

    """
    seqkit grep -n -r -p "ORF type:complete" ${Transdecoder_pep} -o "${params.Sample_name}.trandescoder.complete.pep"
    seqkit grep -n -r -p "ORF type:complete" ${Transdecoder_cds} -o "${params.Sample_name}.trandescoder.complete.cds"
    """

}

// Process 15: SignalP
process SignalP {

    errorStrategy 'ignore'

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Signalp", mode: 'copy'

    input:
    path(transdecodercomplete_pep)

    output:
    path "*_mature.fasta", emit: maturesequences
    path "*.signalp5", emit:signalpsummary


    script:
    

    """
    CURRENT_WORK_DIR="\$PWD"
    echo "\$CURRENT_WORK_DIR"
    input_abs_path=\$(readlink -f "${transdecodercomplete_pep}")
    cd ${params.SIGNALP_PATH}
    ./signalp -fasta "\$input_abs_path" -mature -prefix "${params.Sample_name}"
    echo "\$CURRENT_WORK_DIR"
    mv "${params.Sample_name}_mature.fasta" "\$CURRENT_WORK_DIR"
    mv "${params.Sample_name}_summary.signalp5" "\$CURRENT_WORK_DIR"
    echo "\$CURRENT_WORK_DIR"
    """

}

// Process 16: Filter2
process Filter2 {
    errorStrategy 'ignore'
    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Transdecoder", mode: 'copy'

    input:
    tuple path(maturesequences), path(transdecodercomplete_pep)

    output:
    path "*.fasta", emit: transdecoderpep_signalp


    script:

    """
    seqkit grep -f <(seqkit seq -n ${maturesequences}) ${transdecodercomplete_pep} > ${params.Sample_name}.Transdecoder.complete.signalp.sequences.fasta
    """

}

// Process 17: STATS
process stats {
    errorStrategy 'ignore'
    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Stats", mode: 'copy'

    input:
    tuple path(Transdecoder_pep), path(Transdecoder_cds), path(transdecodercomplete_pep), path(transdecodercomplete_cds), path(maturesequences), path(transdecoderpep_signalp)

    output:
    path "*"


    script:

    """
    seqkit stats ${Transdecoder_pep} > ${params.Sample_name}_Transdecoder_pep.stats.txt
    seqkit stats ${Transdecoder_cds} > ${params.Sample_name}_Transdecoder_cds.stats.txt
    seqkit stats ${transdecodercomplete_pep} > ${params.Sample_name}_complete_pep.stats.txt
    seqkit stats ${transdecodercomplete_cds} > ${params.Sample_name}_transdecodercomplete_cds.stats.txt
    seqkit stats ${maturesequences} > ${params.Sample_name}_maturesequences.stats.txt
        seqkit stats ${transdecoderpep_signalp} > ${params.Sample_name}_transdecoderpep_signalp.stats.txt
    
    """

}

// Process 18: Interproscan
process Interproscan {

    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Interproscan", mode: 'copy'
    
    errorStrategy 'retry'
    maxRetries 2
    cpus { task.attempt * 2 }
    memory { (task.attempt * 2) * 1.9.GB }

    input:
    path (Transdecoder_pep)

    output:
    path "*.tsv", emit: Interproscan


    script:

    """
    awk '{if (\$0 ~ /^>/) print \$0; else {gsub(/\\*/, ""); print \$0}}' ${Transdecoder_pep} > "${params.Sample_name}.Trinity.fasta.transdecoder.cleaned.pep"

    ${params.INTERPROSCAN_PATH} \
      -goterms \
      -i "${params.Sample_name}.Trinity.fasta.transdecoder.cleaned.pep" \
      -pa -t p -d ./ -f TSV
    """

}

// Process 19a: GenomeBlastdatabasecreation
process GenomeBlastdatabasecreation {
    errorStrategy 'ignore'
    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    input:
    path(genome_fasta)

    output:
    path "*", emit: genomedb

    script:
    """

    makeblastdb -in "${genome_fasta}" -dbtype nucl

    """

}

// Process 19: GenomeBlasts
process GenomeBlasts {
    errorStrategy 'ignore'
    conda "${workflow.projectDir}/bin/Setup/VF.yaml"

    publishDir "results/Blast/Blastn/", mode: 'copy'

    input:
    path(transdecodercomplete_cds)
    path(genomedb)

    output:
    path "${params.Sample_name}.blastn.db.0.txt", emit: blastn0
    path "${params.Sample_name}.blastn.db.6.txt", emit: blastn6


    script:
    """
    blastn -query ${transdecodercomplete_cds} -db ${params.genomefastaname} -out ${params.Sample_name}.blastn.db.6.txt -evalue 1e-5 -num_threads 10 -outfmt '6'
    blastn -query ${transdecodercomplete_cds} -db ${params.genomefastaname} -out ${params.Sample_name}.blastn.db.0.txt -evalue 1e-5 -num_threads 10 -outfmt '0'

    """

}


process final {
  
  output:
  path("metadata.txt")

  script:
  """
  cat <<EOF > metadata.txt
  ${params.manifest.author}
  ${params.manifest.version}
  ${workflow.workDir}
  ${workflow.userName}
  ${workflow.start}
  """
}






// PARAMETERS DEFINED IN CONFIG FILE
params.manifest=manifest

//WorkFlow

workflow {
log.info """\
         ${params.manifest.name} v${params.manifest.version}
         ==========================
         input from   : ${params.input_file}
         output to    : ${params.output_dir}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}:${workflow.container}
         """
         .stripIndent()
         
def R1 = Channel.fromPath(params.R1)
def R2 = Channel.fromPath(params.R2)
def R1R2 = R1.combine(R2)
def trinity_fasta = Channel.fromPath(params.trinity_fasta)
R1R2 | PostTrimFastqc

def fastqc_zips = PostTrimFastqc.out.fastqc_zips
def metazoa = Channel.fromPath(params.metazoa)
def mollusca = Channel.fromPath(params.mollusca)
def buscoTranscriptome_met = trinity_fasta.combine(metazoa)
def buscoTranscriptome_mol = trinity_fasta.combine(mollusca)
def TrinityR1R2 = trinity_fasta.combine(R1).combine(R2)
def database_fasta = Channel.fromPath(params.database_fasta)

database_fasta | Blastdatabasecreation

Blastx ( trinity_fasta, Blastdatabasecreation.out.proteindb )

trinity_fasta | Transdecoder
def Transdecoder_pep = Transdecoder.out.pep
def buscoTranslatome_met = Transdecoder_pep.combine(metazoa)
def buscoTranslatome_mol = Transdecoder_pep.combine(mollusca)
def Transdecoder_cds = Transdecoder.out.cds
def TransdecoderR1R2 = Transdecoder_cds.combine(R1).combine(R2)


Blastp ( Transdecoder_pep, Blastdatabasecreation.out.proteindb )


def Transdecoder_pep_cds = Transdecoder_pep.combine(Transdecoder_cds)
Transdecoder_pep_cds | Transdecoder_complete
def transdecodercomplete_pep = Transdecoder_complete.out.transdecodercomplete_pep
transdecodercomplete_pep | SignalP
def maturesequences = SignalP.out.maturesequences
def mature_complete = maturesequences.combine(transdecodercomplete_pep)
def transdecodercomplete_cds = Transdecoder_complete.out.transdecodercomplete_cds
mature_complete | Filter2
def transdecoderpep_signalp = Filter2.out.transdecoderpep_signalp

def stats_combine = Transdecoder_pep
                     .combine(Transdecoder_cds).
                     combine(transdecodercomplete_pep).
                     combine(transdecodercomplete_cds).
                     combine(maturesequences).
                     combine(transdecoderpep_signalp)



fastqc_zips | MultiQC

trinity_fasta | TrinityStats

buscoTranscriptome_met | BUSCO_transcriptome_metazoa

buscoTranscriptome_mol | BUSCO_transcriptome_mollusca


TrinityR1R2 | Kallisto_Trinity



buscoTranslatome_met | BUSCO_translatome_metazoa

buscoTranslatome_mol | BUSCO_translatome_mollusca

TransdecoderR1R2 | Kallisto_Transdecoder



stats_combine | stats

Transdecoder_pep | Interproscan

if (params.genomefasta != 'NULL') {
   def genomefasta = Channel.fromPath(params.genomefasta)
   genomefasta | GenomeBlastdatabasecreation
   GenomeBlasts ( transdecodercomplete_cds, GenomeBlastdatabasecreation.out.genomedb )
}

final()

}
