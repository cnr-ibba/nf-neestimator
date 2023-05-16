# nf-neestimator

A nextflow pipeline for Ne estimates with bootstrapping

## Background

This pipeline is an attempt to call [NeEstimator2](http://www.molecularfisherieslaboratory.com.au/neestimator-software/)
by bootstrapping individuals from a plink binary file. You require both
[nextflow](https://www.nextflow.io/) and [singularity](https://apptainer.org/)
(or [docker](https://www.docker.com/)) configured in order to call this pipeline
properly

## Calling pipeline

Simply type:

```bash
nextflow run cnr-ibba/nf-neestimator -resume -profile singularity --prefix <plink file prefix>
```

where the ``--prefix`` parameter is the same prefix you would pass as plink ``--bfile``

### Customize configuration

This pipeline is configured for bootstrapping 10, 20 and 100 individuals 1000
times. If you need to customize the bootstrapping, you can create a custom
configuration file like this:

```conf
params {
    prefix      = '<plink file prefix>'
    steps       = 50
    individuals = [10, 20, 30]
}
```

Add then provide this custom configuration file with the nextflow `-config`
parameter. See [configuration file](https://www.nextflow.io/docs/latest/config.html#configuration-file)
in the nexflow documentation for more information
