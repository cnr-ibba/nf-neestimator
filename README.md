# nf-rldne
A nextflow pipeline for RLDNe calculations with bootstrapping

## Background

This pipeline is an attempt to call [RLDNe](https://github.com/zakrobinson/RLDNe)
by bootstrapping individuals from a plink binary file. You require both
[nextflow](https://www.nextflow.io/) and [singularity](https://apptainer.org/)
(or [docker](https://www.docker.com/)) configured in order to call this pipeline
properly
