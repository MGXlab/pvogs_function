rule filter_intact:
	input:
		intact_raw = ancient("data/interactions/intact.txt"),
		tax_db = ancient("data/taxonomy_db/taxa.sqlite")
	output:
		intact_phages = "results/interaction_datasets/01_filter_intact/intact_phages.txt"
	conda:
		"../envs/pvogs.yml"
	log:
		"results/logs/construct_datasets/filter_intact.log"
	shell:
		"python workflow/scripts/filter_intact.py "
		"--tax-db {input.tax_db} "
		"--phages-only "
		"-i {input.intact_raw} "
		"-o {output.intact_phages} &>{log}"

rule summarize_intact:
	input:
		intact_phages = rules.filter_intact.output.intact_phages,
		tax_db = ancient("data/taxonomy_db/taxa.sqlite")
	output:
		metadata_tsv = "results/interaction_datasets/02_summarize_intact/metadata.tsv"
	params:
		summary_dir = "results/interaction_datasets/02_summarize_intact"
	conda:
		"../envs/pvogs.yml"
	log:
		"results/logs/construct_datasets/summarize_intact.log"
	shell:
		"python workflow/scripts/summarize_intact.py "
		"--tax-db {input.tax_db} "
		"-i {input.intact_phages} "
		"-o {params.summary_dir} &>{log}"

rule get_uniprot_ids:
	input:
		metadata_tsv = rules.summarize_intact.output.metadata_tsv
	output:
		uniprot_ids_txt = "results/interaction_datasets/03_uniprot/uniprot_ids.txt",
		uniprot_only = "results/interaction_datasets/03_uniprot/uniprot.phages_only.tsv"
	conda:
		"../envs/pvogs.yml"
	log:
		"results/logs/construct_datasets/get_uniprot_ids.log"
	shell:
		"python workflow/scripts/get_uniprot_ids.py "
		"-i {input.metadata_tsv} "
		"-l {output.uniprot_ids_txt} "
		"-o {output.uniprot_only} &>{log}"
	
rule download_uniprots:
	input:
		uniprot_ids_txt = rules.get_uniprot_ids.output.uniprot_ids_txt
	output:
		swissprot_txt = "results/interaction_datasets/03_uniprot/uniprot.sprot.txt"
	log:
		"results/logs/construct_datasets/query_uniprot.log"
	conda:
		"../envs/pvogs.yml"
	shell:
		"python workflow/scripts/query_uniprot.py "
		"--filter-list "
		"-i {input.uniprot_ids_txt} "
		"-o {output.swissprot_txt} &>{log}"

rule process_uniprot:
	input:
		uniprot_ids_txt = rules.get_uniprot_ids.output.uniprot_ids_txt,
		uniprot_only = rules.get_uniprot_ids.output.uniprot_only,
		swissprot_txt = rules.download_uniprots.output.swissprot_txt
	output:
		proteins_fasta = "results/interaction_datasets/04_process_uniprot/proteins.faa",
		skipped = "results/interaction_datasets/04_process_uniprot/skipped.txt",
		duplicates = "results/interaction_datasets/04_process_uniprot/duplicates.tsv",
		interactions_filtered = "results/interaction_datasets/04_process_uniprot/interactions_filtered.tsv",
		ncbi_interactions = "results/interaction_datasets/04_process_uniprot/ncbi_interactions.tsv",
		genomes_accessions = "results/interaction_datasets/04_process_uniprot/genome_accessions.txt",
		ncbi_to_uniprot = "results/interaction_datasets/04_process_uniprot/ncbi2uniprot.mapping.txt",
		uniprot_to_ncbi = "results/interaction_datasets/04_process_uniprot/uniprot2ncbi.mapping.txt"
	conda:
		"../envs/pvogs.yml"
	log:
		"results/logs/construct_datasets/process_uniprot.log"
	params:
		output_dir = "results/interaction_datasets/04_process_uniprot"
	shell:
		"python workflow/scripts/process_uniprot.py "
		"-i {input.uniprot_only} "
		"-l {input.uniprot_ids_txt} "
		"-s {input.swissprot_txt} "
		"-p {params.output_dir} &>{log}"

rule download_genomes:
	input:
		genomes_accessions = rules.process_uniprot.output.genomes_accessions
	output:
		genomes_gb = "results/interaction_datasets/05_interaction_datasets/genomes.gb"
	log:
		"results/logs/construct_datasets/download_genomes.log"
	params:
		email = config.get('email')
	conda:
		"../envs/pvogs.yml"
	shell:
		"python workflow/scripts/download_genomes.py "
		"-i {input.genomes_accessions} "
		"-o {output.genomes_gb} "
		"-e {params.email} &>{log}"

rule extract_proteins_from_genomes:
	input:
		genomes_gb = rules.download_genomes.output.genomes_gb
	output:
		all_proteins_faa = "results/interaction_datasets/05_genomes/all_proteins.faa"
	conda:
		"../envs/pvogs.yml"
	log:
		"results/logs/interaction_datasets/05_genomes/.extract_proteins.log"
	shell:
		"python workflow/scripts/extract_proteins_from_gb.py "
		"-i {input.genomes_gb} "
		"-o {output.all_proteins_faa} &>{log}"

## Copy the proteins file and create a 2-column tsv in a positives directory
## The same naming structure helps with expansion for the rules.
rule copy_positives_proteins:
	input:
		proteins_faa = rules.process_uniprot.output.proteins_fasta
	output:
		positives_faa = "results/interaction_datasets/positives/positives.proteins.faa"
	log:
		"results/logs/construct_datasets/copy_positives_scores.log"
	shell:
		"cp {input.proteins_faa} {output.positives_faa} 2>{log}"


rule ncbi_positives:
	input:
		ncbi_interactions = rules.process_uniprot.output.ncbi_interactions,
		genomes_gb = rules.download_genomes.output.genomes_gb
	output:
		ncbi_positives_tsv = "results/interaction_datasets/positives/positives.interactions.tsv"
	log:
		"results/logs/construct_datasets/ncbi_positives.log"
	shell:
		"cut -f3,4 {input.ncbi_interactions} > {output.ncbi_positives_tsv} 2>{log}"


rule make_negatives:
	input:
		ncbi_positives_tsv = rules.ncbi_positives.output.ncbi_positives_tsv,
		all_proteins_faa = rules.extract_proteins_from_genomes.output.all_proteins_faa,
		genomes_gb = rules.download_genomes.output.genomes_gb

	output:
		expand(["results/interaction_datasets/N{I}/N{I}.interactions.tsv",
				"results/interaction_datasets/N{I}/N{I}.proteins.faa"],
				I = [i+1 for i in range(0, NEGATIVES)])
	conda:
		"../envs/pvogs.yml"
	params:
		no_negatives = NEGATIVES
	log:
		"results/logs/construct_datasets/make_negatives.log"
	shell:
		"for i in `seq 1 {params.no_negatives}`;do "
		"	python workflow/scripts/make_protein_combos.py "
		"	-i {input.genomes_gb} "
		"	-a {input.all_proteins_faa} "
		"	-o results/interaction_datasets/N${{i}}/N${{i}}.proteins.faa "
		"	-x results/interaction_datasets/N${{i}}/N${{i}}.interactions.tsv "
		"	--exclude {input.ncbi_positives_tsv} "
		"	--sample-size 2 "
		"	--random-seed ${{i}}; "
		"done &>{log}"	

rule hmmsearch:
	input:
		proteins_fasta = "results/interaction_datasets/{dataset}/{dataset}.proteins.faa",
		all_pvogs_profiles = "data/pvogs/all.hmm"
	output:
		hmm_out_txt = "results/interaction_datasets/06_map_proteins_to_pvogs/{dataset}/{dataset}.hmmout.txt",
		hmm_tblout_tsv = "results/interaction_datasets/06_map_proteins_to_pvogs/{dataset}/{dataset}.hmmtblout.tsv" 
	log:
		"results/logs/construct_datasets/{dataset}/hmmsearch.log"
	threads: 8
	conda:
		"../envs/pvogs.yml"
	shell:
		"hmmsearch --cpu {threads} "
		"-o {output.hmm_out_txt} "
		"--tblout {output.hmm_tblout_tsv} "
		"{input.all_pvogs_profiles} "
		"{input.proteins_fasta}"

rule refseqs_to_pvogs:
	input:
		interactions_tsv = "results/interaction_datasets/{dataset}/{dataset}.interactions.tsv",
		proteins_faa = "results/interaction_datasets/{dataset}/{dataset}.proteins.faa",
		hmm_tblout = "results/interaction_datasets/06_map_proteins_to_pvogs/{dataset}/{dataset}.hmmtblout.tsv"
	output:
		pvogs_interactions = "results/interaction_datasets/{dataset}/{dataset}.pvogs_interactions.tsv"
	conda:
		"../envs/pvogs.yml"
	log:
		"results/logs/construct_datasets/{dataset}/refseqs_to_pvogs.log"
	shell:
		"python workflow/scripts/refseqs_to_pvogs.py "
		"-i {input.interactions_tsv} "
		"-f {input.proteins_faa} "
		"-hmm {input.hmm_tblout} "
		"-o {output.pvogs_interactions} "
		"&> {log}"

rule create_features_tables:
	input:
		filtered_master_tsv = rules.filter_scores_table.output.filtered_master_tsv,
		interactions_tsv = "results/interaction_datasets/{dataset}/{dataset}.pvogs_interactions.tsv"
	output:
		features_tsv = "results/interaction_datasets/{dataset}/{dataset}.features.tsv"
	log:
		"results/logs/construct_datasets/{dataset}/create_features_tables.log"
	conda:
		"../envs/pvogs.yml"
	shell:
		"python workflow/scripts/subset_scores.py "
		"-s {input.filtered_master_tsv} "
		"-i {input.interactions_tsv} "
		"-o {output.features_tsv}"

