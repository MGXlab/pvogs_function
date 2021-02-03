rule download_archive:
    output:
        tar_gz = "pvogs_function.data.tar.gz"
    log:
        "results/.logs/get_data/download_archive.log"
    conda:
        "envs/zenodo.yml"
    params:
        sandbox = "--sandbox" if config["zenodo"]["use_sandbox"] is True else "",
        zenodo_id = config["zenodo"]["doi"]
    shell:
        "zenodo_get {params.sandbox} --record={params.zenodo_id}"

rule extract_data:
    input:
        tar_gz = rules.download_archive.output.tar_gz
    output:
        genomes_fasta = "data/genomes/phages_refseq.fasta",
        intact_txt= "data/interactions/intact.txt",
        pvogs_profiles = "data/pvogs/all.hmm",
        pvogs_annotations = "data/pvogs/VOGProteinTable.txt",
        taxonomy_db = "data/taxonomy_db/taxa.sqlite",
        taxonomy_pkl = "data/taxonomy_db/taxa.sqlite.traverse.pkl"
    log:
        "results/.logs/get_data/extract_data.log"
    shell:
        "tar -xzvf {input.tar_gz} &>{log}"

