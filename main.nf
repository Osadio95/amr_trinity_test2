#!/usr/bin/env nextflow

// Inclure vos modules EXISTANTS - NE PAS MODIFIER
include { amrfinderplus } from './modules/amrfinderplus.nf'
include { resfinder } from './modules/resfinder.nf'
include { rgi } from './modules/rgi.nf'

// Processus hamronize - MINIMAL
process hamronize {
    publishDir "${params.out_dir}/harmonized", mode: 'copy'
    
    input:
    tuple val(id), val(tool), path(metadata), path(results)
    
    output:
    path "${id}-${tool}.tsv"
    
    script:
    """
    METADATA=\$(cat $metadata)
    hamronize \$tool "\$METADATA" $results > "${id}-${tool}.tsv"
    """
}

// Processus summarise - MINIMAL
process summarise {
    publishDir "${params.out_dir}/summary", mode: 'copy'
    
    input:
    path inputs
    
    output:
    path('report.tsv')
    path('report.json')
    path('report.html')
    
    script:
    """
    hamronize summarize -t tsv -o report.tsv $inputs
    hamronize summarize -t json -o report.json $inputs
    hamronize summarize -t interactive -o report.html $inputs
    """
}

// Workflow principal - MINIMAL
workflow {
    // NOTE IMPORTANTE : Vos modules attendent 'contigs' mais le schema EPI2ME utilise 'assembly'
    // On adapte le nom ici seulement
    Channel.fromPath(params.assemblies_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.id, row.species, file(row.assembly)) }
        | (amrfinderplus & resfinder & rgi)
    
    Channel.of().mix(amrfinderplus.out, resfinder.out, rgi.out)
        | hamronize 
        | collect
        | summarise
}
