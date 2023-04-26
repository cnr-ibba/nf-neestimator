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
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.bim"), emit: bim
    tuple val(meta), path("*.fam"), emit: fam
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
        --make-bed \\
        --threads $task.cpus \\
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
    touch ${outfile}.bed
    touch ${outfile}.bim
    touch ${outfile}.fam
    touch version.yml
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

    PLINK_SUBSET(iterations_ch, plink_input_ch, params.species_opts)
}

workflow {
    RNELD()
}
