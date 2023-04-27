#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.prefix = 'testdata/murciano_test'
params.species_opts = '--allow-extra-chr --chr-set 29'

steps_ch = Channel.from( 1..5 )
individuals_ch = Channel.of( 20, 50, 100 )


process PLINK_SUBSET {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink=1.90b6.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), val(iteration)
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
    def outfile = "${bed.getBaseName()}_${individuals}_individuals_step_${step}"
    """
    plink \\
        $species_opts \\
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
    def outfile = "${bed.getBaseName()}_${individuals}_individuals_step_${step}"
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
    tuple val(meta), path("*.txt"), emit: map
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


workflow RNELD {
    iterations_ch = individuals_ch.combine(steps_ch)//.view()
        .map{ iteration -> [[
            id:"individuals_${iteration[0]}_step_${iteration[1]}",
            individuals: "${iteration[0]}",
            step: "${iteration[1]}"], iteration ]}//.view()

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
}


workflow {
    RNELD()
}
