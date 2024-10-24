process MULTIQC {

    if (workflow.containerEngine == 'singularity') {
        container = params.multiqc_singularity
    } else {
        container = params.multiqc_docker
    }

    publishDir "${params.output}/multiqc", mode: 'copy'

    input:
    path multiqc_files, stageAs: "?/*"
    path assets_dir

    output:
    path "*.html", emit: report
    path "*_data", emit: data

    """
    multiqc -c ${assets_dir}/multiqc_config.yml .
    mv *.html multiqc_${params.exp}.html
    mv *_data multiqc_${params.exp}_data
    """
}
