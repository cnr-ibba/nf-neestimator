#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// check parameters
if (!params.prefix) { exit 1, "Error: 'prefix' parameter not specified" }
if (params.steps) {
    def range = (1..params.steps)
    steps_ch = Channel.of( range )
} else {
    exit 1, "Error: 'steps' parameter not specified"
}
if (params.individuals) { individuals_ch = Channel.fromList( params.individuals ) } else { exit 1, "Error: 'individuals' parameter not specified" }


process PLINK_SUBSET {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink=1.90b6.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), val(seed)
    tuple path(bed), path(bim), path(fam)
    val(species_opts)

    output:
    tuple val(meta), path("*.map"), emit: map
    tuple val(meta), path("*.ped"), emit: ped
    path "versions.yml"           , emit: versions

    script:
    def step = "${meta.step}"
    def individuals = "${meta.individuals}"
    def prefix = "${bed.getBaseName()}"
    def outfile = "${bed.getBaseName()}_${individuals}_individuals_${step}_step"
    """
    plink \\
        $species_opts \\
        --seed $seed \\
        --bfile $prefix \\
        --thin-indiv-count ${meta.individuals} \\
        --threads $task.cpus \\
        --recode \\
        --out ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """

    stub:
    def step = "${meta.step}"
    def individuals = "${meta.individuals}"
    def prefix = "${bed.getBaseName()}"
    def outfile = "${bed.getBaseName()}_${individuals}_individuals_${step}_step"
    """
    touch ${outfile}.map
    touch ${outfile}.ped
    touch versions.yml
    """
}


process PED2GENEPOP {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pgdspider=2.1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pgdspider:2.1.1.5--hdfd78af_1' :
        'quay.io/biocontainers/pgdspider:2.1.1.5--hdfd78af_1' }"

    input:
    tuple val(meta), path(ped)
    tuple val(meta), path(map)

    output:
    tuple val(meta), path("*.txt"), emit: genepop
    tuple val(meta), path("*.spid"), emit: spid
    path "versions.yml"           , emit: versions

    script:
    def prefix = "${ped.getBaseName()}"
    def xmx = (task.memory.mega*0.95).intValue()
    def xms = task.memory.mega / 2
    """
    cat <<EOF > PED_GENEPOP.spid
    # PED Parser questions
    PARSER_FORMAT=PED

    # Is the Phenotype absent in the input file?
    PED_PARSER_PHENOTYPE_QUESTION=false
    # Is the Family ID column absent in the input file?
    PED_PARSER_FAMILY_ID_QUESTION=false
    # Do you want to include a MAP file with loci information?
    PED_PARSER_INCLUDE_MAP_QUESTION=true
    # Is the Sex column absent in the input file?
    PED_PARSER_SEX_QUESTION=false
    # Are the Paternal ID and the Maternal ID columns absent in the input file?
    PED_PARSER_PATERNAL_MATERNAL_ID_QUESTION=false
    # Open MAP file
    PED_PARSER_MAP_FILE_QUESTION=${map}
    # Group individuals into populations according to:
    PED_PARSER_POPULATION_QUESTION=FAMILY
    # Is the Individual ID column absent in the input file?
    PED_PARSER_IND_ID_QUESTION=false

    # GENEPOP Writer questions
    WRITER_FORMAT=GENEPOP

    # Specify the locus/locus combination you want to write to the GENEPOP file:
    GENEPOP_WRITER_LOCUS_COMBINATION_QUESTION=
    # Specify which data type should be included in the GENEPOP file  (GENEPOP can only analyze one data type per file):
    GENEPOP_WRITER_DATA_TYPE_QUESTION=SNP
    EOF

    PGDSpider2-cli \\
        -Xmx${xmx}M \\
        -Xms${xms}M \\
        -inputfile ${ped} \\
        -inputformat PED \\
        -outputfile ${prefix}_genepop.txt \\
        -outputformat GENEPOP \\
        -spid PED_GENEPOP.spid

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(PGDSpider2-cli -h) | grep -i Copyright | sed 's/ PGDSpider //;s/Copyright.*//')
    END_VERSIONS
    """

    stub:
    def prefix = "${ped.getBaseName()}"
    """
    touch ${prefix}_genepop.txt
    touch PED_GENEPOP.spid
    touch versions.yml
    """
}


process LDNE {
    tag "$meta.id"
    label 'process_single'
    label 'error_retry'

    container "bunop/neestimator2x:0.2"

    input:
    tuple val(meta), path(genepop)

    output:
    tuple val(meta), path("*_info.txt"), emit: info
    tuple val(meta), path("*_option.txt"), emit: option
    tuple val(meta), path("*_Ne.txt"), emit: output
    tuple val(meta), path("*_NexLD.txt"), emit: tabular
    tuple val(meta), path("*NoDat.txt"), emit: missing

    script:
    def prefix = "${genepop.getBaseName()}"
    def info_file = "${prefix}_info.txt"
    def option_file = "${prefix}_option.txt"
    def output_file = "${prefix}_Ne.txt"
    def missing_file = "${prefix}NoDat.txt"
    """
    cat <<EOF > ${info_file}
    1           * A number = sum of method(s) to run: LD(=1), Het(=2), Coan(=4), Temporal(=8).
    \$PWD/      * Input Directory
    ${genepop}  * Input file name
    2           * 1 = FSTAT format, 2 = GENEPOP format
    \$PWD/      * Output Directory
    ${output_file}*     * Output file name (put asterisk adjacent to the name to append)
    1           * Number of critical values, added 1 if a run by rejecting only singleton alleles is included
    0.02 -1     * Critical values, a special value '1' is for rejecting only singleton alleles
    1           * 0: Random mating, 1: Monogamy (LD method)
    EOF

    cat <<EOF > ${option_file}
    1  0  1  1  * First number = sum of method(s) to have extra output: LD(=1), Het(=2), Coan(=4), Temporal(=8)
    0           * Maximum individuals/pop. If 0: no limit
    0           * First entry n1 = 0: No Freq output. If n1 = -1: Freq. output up to population 50. Two entries n1, n2 with n1 <= n2: Freq output for populations from n1 to n2. Max. populations to have Freq output is set at 50
    0           * For Burrow output file (up to 50 populations can have output). See remark below
    1           * Parameter CI: 1 for Yes, 0 for No
    1           * Jackknife CI: 1 for Yes, 0 for No
    0           * Up to population, or range of populations to run (if 2 entries). If first entry = 0: no restriction
    0           * All loci are accepted
    1           * Enter 1: A file is created to document missing data if there are any in input file. Enter 0: no file created
    0           * Line for chromosomes/loci option and file
    EOF

    Ne2-1L i:${info_file} o:${option_file}
    """

    stub:
    def prefix = "${genepop.getBaseName()}"
    """
    touch ${prefix}_info.txt
    touch ${prefix}_option.txt
    touch ${prefix}_Ne.txt
    touch ${prefix}_NexLD.txt
    touch ${prefix}NoDat.txt
    """
}


process SUMMARIZE {
    tag "$meta individuals"
    label 'process_low'

    conda "biopython:1.78"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'quay.io/biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(ne_output)
    tuple path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*_individuals.csv"), emit: csv

    script:
    def prefix = "${bed.getBaseName()}"
    def summary = "${prefix}_Ne_${meta}_individuals.csv"
    """
    parse_ne_results.py > ${summary}
    """

    stub:
    def summary = "Ne_${meta}_individuals.csv"
    """
    touch ${summary}
    """
}


workflow LDNE_PIPELINE {
    iterations_ch = individuals_ch.combine(steps_ch)//.view()
        .map{ iteration -> [[
            id:"${iteration[0]}_individual_${iteration[1]}_step",
            individuals: "${iteration[0]}",
            step: "${iteration[1]}"], iteration[0] * iteration[1] ]}//.view()

    plink_input_ch = Channel.fromPath( "${params.prefix}.{bim,bed,fam}")
        .collect()//.view()

    // subsetting plink files
    PLINK_SUBSET(iterations_ch, plink_input_ch, params.species_opts)

    // create GENEPOP file
    PED2GENEPOP(PLINK_SUBSET.out.ped, PLINK_SUBSET.out.map)

    // launch LDNE
    LDNE(PED2GENEPOP.out.genepop)

    neestimator_ch = LDNE.out.tabular
        .map { meta, path -> [
            meta["individuals"],
            path
        ]}
        .groupTuple( by: 0, size: params.steps )
        //.view()

    // collect and parse results
    SUMMARIZE(neestimator_ch, plink_input_ch)
}


workflow {
    LDNE_PIPELINE()
}
