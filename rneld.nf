#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.prefix = 'murciano_gmmphk'
params.species = '--sheep'

steps_ch = Channel.from( 1..5 )
individuals_ch = Channel.of( 20, 50, 100 )


process PLINK_SUBSET {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), val(iteration)
    val(species)
    val(prefix)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.bim"), emit: bim
    tuple val(meta), path("*.fam"), emit: fam

    script:
    def step = "${meta.step}"
    def individuals = "${meta.individuals}"
    def outfile = "${prefix}_${individuals}_${step}"
    """
    plink \\
        $species \\
        --bfile $prefix \\
        --thin-indiv-count ${meta.individuals} \\
        --out ${outfile} \\
        --make-bed >> ${meta.id}.txt
    """

    stub:
    def step = "${meta.step}"
    def individuals = "${meta.individuals}"
    def outfile = "${prefix}_${individuals}_${step}"
    """
    touch ${outfile}.bed
    touch ${outfile}.bim
    touch ${outfile}.fam
    """
}

workflow RNELD {
    iterations_ch = individuals_ch.combine(steps_ch)//.view()
        .map{ iteration -> [[
            id:"individuals_${iteration[0]}_step_${iteration[1]}",
            individuals: "${iteration[0]}",
            step: "${iteration[1]}"], iteration ]}//.view()

    PLINK_SUBSET(iterations_ch, params.species, params.prefix)
}

workflow {
    RNELD()
}
