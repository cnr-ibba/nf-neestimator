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
    tuple val(meta), val(ped), path(spid)

    output:
    tuple val(meta), path("*.txt"), emit: genepop
    path "versions.yml"           , emit: versions

    script:
    def prefix = "${ped.getBaseName()}"
    def xmx = (task.memory.mega*0.95).intValue()
    def xms = task.memory.mega / 2
    """
    PGDSpider2-cli \\
        -Xmx${xmx}M \\
        -Xms${xms}M \\
        -inputfile $ped \\
        -inputformat PED \\
        -outputfile ${prefix}_genepop.txt \\
        -outputformat GENEPOP \\
        -spid $spid

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(PGDSpider2-cli -h) | grep -i Copyright | sed 's/ PGDSpider //;s/Copyright.*//')
    END_VERSIONS
    """

    stub:
    def prefix = "${ped.getBaseName()}"
    """
    touch ${prefix}_genepop.txt
    touch versions.yml
    """
}


process RLDNE {
    tag "$meta.id"
    label 'process_single'
    label 'error_retry'

    container "bunop/rldne:0.2"

    input:
    tuple val(meta), path(genepop)

    output:
    tuple val(meta), path("*_Ne_params.txt"), emit: params
    tuple val(meta), path("*_Ne_out.txt"), emit: file
    tuple val(meta), path("*_Ne_outxLD.txt"), emit: tab

    script:
    def prefix = "${genepop.getBaseName()}"
    def ne_out_tab = "${prefix}_Ne_outxLD.txt"
    def param_file = "${prefix}_Ne_params.txt"
    def ne_out_file = "${prefix}_Ne_out.txt"
    """
    #!/usr/bin/env Rscript

    library(RLDNe)

    param_files <- NeV2_LDNe_create(
        input_file="${genepop}",
        param_file="${param_file}",
        NE_out_file="${ne_out_file}",
        matingsystem = 1,
        crit_vals = 0.02
    )

    run_LDNe(LDNe_params = param_files\$param_file)
    """

    stub:
    def prefix = "${genepop.getBaseName()}"
    """
    touch ${prefix}_Ne_out.txt
    touch ${prefix}_Ne_params.txt
    touch ${prefix}_Ne_outxLD.txt
    """

}


workflow RLDNE_PIPELINE {
    iterations_ch = individuals_ch.combine(steps_ch)//.view()
        .map{ iteration -> [[
            id:"${iteration[0]}_individual_${iteration[1]}_step",
            individuals: "${iteration[0]}",
            step: "${iteration[1]}"], iteration[0] * iteration[1] ]}//.view()

    plink_input_ch = Channel.fromPath( "${params.prefix}.{bim,bed,fam}")
        .collect()//.view()

    // subsetting plink files
    PLINK_SUBSET(iterations_ch, plink_input_ch, params.species_opts)

    // get SPID path
    pgdspider_spid_ch = Channel.fromPath(params.pgdspider_spid)

    // a new channel for data conversion
    pgdspider_input_ch = PLINK_SUBSET.out.ped.combine(pgdspider_spid_ch)//.view()

    // create GENEPOP file
    PED2GENEPOP(pgdspider_input_ch)

    // launch RLDNE
    RLDNE(PED2GENEPOP.out.genepop)
}


workflow {
    RLDNE_PIPELINE()
}
