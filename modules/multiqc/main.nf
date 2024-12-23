process MULTIQC {

    container workflow.containerEngine == 'singularity' ?
        params.multiqc_singularity : params.multiqc_docker

    publishDir "${params.output}/multiqc", mode: 'copy'

    input:
    path multiqc_files, stageAs: "?/*"

    output:
    path "*.html", emit: report
    path "*_data", emit: data

    """
    cp ${projectDir}/assets/* .
    multiqc -c ./multiqc_config.yml .
    mv *.html multiqc_${params.exp}.html
    mv *_data multiqc_${params.exp}_data
    """
}
