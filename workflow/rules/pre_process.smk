## pre_process
## PREPARE NECESSARY FILES FOR FEATURE CALCULATIONS

rule translate_genomes:
	input:
		refseq_genomes = ancient("data/genomes/phages_refseq.fasta")
	output:
		transeq_genomes = "results/pre_process/transeq/transeq.genomes.fasta"
	log: "results/logs/pre_process/transeq/transeq_genomes.log"
	conda:
		"../envs/pvogs.yml"
	shell:
		"transeq -frame 6 "
		"-table 11 -clean "
		"-sequence {input.refseq_genomes} "
		"-outseq {output.transeq_genomes} 2>{log}"

rule hmmsearch_transeq:
	input:
		transeq_genomes = rules.translate_genomes.output.transeq_genomes,
		pvogs_all_profiles = ancient("data/pvogs/all.hmm")
	output:
		hmmout_txt = "results/pre_process/hmmsearch/transeq.hmmout.txt",
		hmmtblout_tsv = "results/pre_process/hmmsearch/transeq.hmmtblout.tsv"
	log:
		"results/logs/pre_process/hmmsearch/hmmsearch_transeq.log"
	conda:
		"../envs/pvogs.yml"
	threads:
		config['hmmsearch_transeq'].get('threads', 10)
	shell:
		"hmmsearch --cpu {threads} "
		"-o {output.hmmout_txt} "
		"--tblout {output.hmmtblout_tsv} "
		"{input.pvogs_all_profiles} "
		"{input.transeq_genomes} "
		"&>{log}"


rule split_genomes:
	input:
		refseq_genomes = ancient("data/genomes/phages_refseq.fasta")
	output:
		reflist_txt = "results/pre_process/reflist.txt"
	log:
		"results/logs/pre_process/split_genomes.log"
	conda:
		"../envs/pvogs.yml"
	params:
		genomes_dir = "results/pre_process/all_genomes"
	shell:
		"mkdir -p {params.genomes_dir} && "
		"python workflow/scripts/split_multifasta.py "
		"--write-reflist "
		"-i {input.refseq_genomes} "
		"-o {params.genomes_dir} 2>{log}"


## FASTANI
rule fastani:
	input:
		reflist_txt = rules.split_genomes.output.reflist_txt
	output:
		fastani_raw = "results/pre_process/fastani/fastani.out",
		fastani_matrix = "results/pre_process/fastani/fastani.out.matrix"
	log:
		"results/logs/pre_process/fastani/fastani.log"
	conda:
		"../envs/pvogs.yml"
	threads: 
		config['fastani'].get('threads', 8)
	params:
		fragLen = 300,
		minFraction = 0.1
	shell:
		"fastANI -t {threads} "
		"--ql {input.reflist_txt} "
		"--rl {input.reflist_txt} "
		"--fragLen {params.fragLen} "
		"--minFraction {params.minFraction} "
		"-o {output.fastani_raw} "
		"--matrix 2>{log}" # Output the matrix

rule fastani_matrix_to_square:
	input:
		fastani_matrix = rules.fastani.output.fastani_matrix
	output:
		fastani_square_mat = "results/pre_process/fastani/fastani.square.matrix"
	conda:
		"../envs/pvogs.yml"
	log:
		"results/logs/pre_process/fastani_matrix_to_square.log"
	shell:
		"python workflow/scripts/fastani_mat_to_square.py "
		"-i {input.fastani_matrix} "
		"-o {output.fastani_square_mat} "
		"--process-names &>{log}" # Chop full paths to last node, which is the accession number



## COMPAREM
rule comparem_call_genes:
	input:
		rules.split_genomes.output.reflist_txt
	output:
		done_file = touch("results/pre_process/comparem/aai_wf/genes/.done.txt")
	conda:
		"../envs/pvogs.yml"
	log:
		# This is produced by comparem by default
		"results/pre_process/comparem/aai_wf/genes/comparem.log"
	params:
		genes_dir = "results/pre_process/comparem/aai_wf/genes",
		genomes_dir = "results/pre_process/all_genomes"
	threads:
		config['comparem_call_genes'].get('threads', 10)
	shell:
		"comparem call_genes -c {threads} --silent -x fasta "
		"{params.genomes_dir} "
		"{params.genes_dir}"

rule remove_empty_files:
	input:
		rules.comparem_call_genes.output.done_file
	output:
		info_file = "results/pre_process/comparem/aai_wf/genes/genomes_skipped.txt"
	params:
		genes_dir = "results/pre_process/comparem/aai_wf/genes"
	log:
		"results/logs/pre_process/remove_empty_files.log"
	shell:
		"python workflow/scripts/remove_empty_files.py "
		"-i {params.genes_dir} "
		"-o {output.info_file}"

rule comparem_similarity:
	input:
		info_file = rules.remove_empty_files.output.info_file
	output:
		hits_sorted_tsv = "results/pre_process/comparem/aai_wf/similarity/hits_sorted.tsv",
		query_genes_dmnd = "results/pre_process/comparem/aai_wf/similarity/query_genes.dmnd",
		query_genes_faa = "results/pre_process/comparem/aai_wf/similarity/query_genes.faa"
	params:
		similarity_dir = "results/pre_process/comparem/aai_wf/similarity",
		genes_dir = "results/pre_process/comparem/aai_wf/genes"
	log:
		# Produced by comparem by default
		"results/pre_process/comparem/aai_wf/similarity/comparem.log"
	threads:
		config['comparem_similarity'].get('threads', 10)
	conda:
		"../envs/pvogs.yml"
	shell:
		"comparem similarity "
		"-c {threads} -x faa "
		"--silent "
		"{params.genes_dir} {params.genes_dir} {params.similarity_dir}"


rule comparem_aai:
	input:
		query_genes_faa = rules.comparem_similarity.output.query_genes_faa,
		hits_sorted_tsv = rules.comparem_similarity.output.hits_sorted_tsv
	output:
		aai_summary_tsv = "results/pre_process/comparem/aai_wf/aai/aai_summary.tsv"
	conda:
		"../envs/pvogs.yml"
	log:
		# Produced by comparem by default
		"results/pre_process/comparem/aai_wf/aai/comparem.log"
	params:
		aai_dir = "results/pre_process/comparem/aai_wf/aai"
	threads:
		config['comparem_aai'].get('threads', 10)
	shell:
		"comparem aai --silent "
		"-c {threads} "
		"{input.query_genes_faa} {input.hits_sorted_tsv} "
		"{params.aai_dir}"

rule process_comparem_matrix:
	input:
		aai_summary_tsv = rules.comparem_aai.output.aai_summary_tsv
	output:
		aai_summary_square = "results/pre_process/comparem/aai_wf/aai/aai_summary.square.tsv"
	conda:
		"../envs/pvogs.yml"
	log:
		"results/logs/pre_process/proecess_comparem_matrix.log"
	shell:
		"python workflow/scripts/process_comparem.py "
		"-i {input.aai_summary_tsv} "
		"-o {output.aai_summary_square} &>{log}"

rule calculate_all_scores:
	input:
		ani_square_matrix = rules.fastani_matrix_to_square.output.fastani_square_mat,
		aai_square_matrix = rules.process_comparem_matrix.output.aai_summary_square,
		phages_genomes_fasta = ancient("data/genomes/phages_refseq.fasta"),
		pvogs_all_profiles = ancient("data/pvogs/all.hmm"),
		hmmout_txt = rules.hmmsearch_transeq.output.hmmout_txt
	output:
		master_table = "results/scores.tsv"
	conda: "../envs/pvogs.yml"
	log: 
		"results/logs/pre_process/calculate_scores.log"
	shell:
		"python workflow/scripts/calculate_all_scores.py "
		"--profiles-file {input.pvogs_all_profiles} "
		"--genomes {input.phages_genomes_fasta} "
		"--input-hmm {input.hmmout_txt} "
		"--ani-matrix {input.ani_square_matrix} "
		"--aai-matrix {input.aai_square_matrix} "
		"-o {output.master_table} &>{log}"

rule filter_scores_table:
	input:
		master_table = rules.calculate_all_scores.output.master_table
	output:
		filtered_master_tsv = "results/filtered_scores.tsv"
	log:
		"results/pre_process/filter_scores_table.log"
	shell:
		"""awk '{{if ( $10 != "1000000") print $0}}' {input.master_table} > {output.filtered_master_tsv}"""

rule parse_annotations:
	input:
		raw_annotations_fp = ancient("data/pvogs/VOGProteinTable.txt")
	output:
		processed_annotations_fp = "results/annotations.tsv"
	conda:
		"../envs/pvogs.yml"
	log:
		"results/pre_process/process_annotations.log"
	shell:
		"python workflow/scripts/process_annotations.py "
		"-i {input.raw_annotations_fp} "
		"-o {output.processed_annotations_fp} "
		"&>{log}"


