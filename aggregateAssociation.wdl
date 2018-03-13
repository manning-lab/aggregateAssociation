task fitNull {
	File genotype_file
	File? phenotype_file
	String? outcome_name
	String? outcome_type
	String? covariates_string
	File? sample_file
	String label
	File? kinship_matrix
	String? id_col

	Int memory
	Int disk

	command {
		echo "Input files" > fitNull.log
		echo "genotype_file: ${genotype_file}" >> fitNull.log
		echo "phenotype_file: ${phenotype_file}" >> fitNull.log
		echo "outcome_name: ${outcome_name}" >> fitNull.log
		echo "outcome_type: ${outcome_type}" >> fitNull.log
		echo "covariates_string: ${covariates_string}" >> fitNull.log
		echo "sample_file: ${sample_file}" >> fitNull.log
		echo "label: ${label}" >> fitNull.log
		echo "kinship_matrix: ${kinship_matrix}" >> fitNull.log
		echo "id_col: ${id_col}" >> fitNull.log
		echo "memory: ${memory}" >> fitNull.log
		echo "disk: ${disk}" >> fitNull.log
		echo "" >> fitNull.log
		dstat -c -d -m --nocolor 10 1>>fitNull.log &
		R --vanilla --args ${genotype_file} ${phenotype_file} ${outcome_name} ${outcome_type} ${default="NA" covariates_string} "NA" "NA" ${sample_file} ${label} ${kinship_matrix} ${id_col} < /aggregateAssociation/genesis_nullmodel.R
	}

	runtime {
		docker: "manninglab/aggregateassociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File model = "${label}_null.RDa"
		File log = "fitNull.log"
	}
}

task aggAssocTest {
	File gds_file
	File null_file
	File group_file
	String label
	String? test
	String? pval
	String? weights

	Int memory
	Int disk

	command {
		echo "Input files" > aggAssocTest.log
		echo "gds_file: ${gds_file}" >> aggAssocTest.log
		echo "null_file: ${null_file}" >> aggAssocTest.log
		echo "group_file: ${group_file}" >> aggAssocTest.log
		echo "label: ${label}" >> aggAssocTest.log
		echo "test: ${test}" >> aggAssocTest.log
		echo "pval: ${pval}" >> aggAssocTest.log
		echo "weights: ${weights}" >> aggAssocTest.log
		echo "memory: ${memory}" >> aggAssocTest.log
		echo "disk: ${disk}" >> aggAssocTest.log
		echo "" >> aggAssocTest.log
		dstat -c -d -m --nocolor 10 1>>aggAssocTest.log &
		R --vanilla --args ${gds_file} ${null_file} ${group_file} ${label} ${default="SKAT" test} ${default="kuonen" pval} ${default="1,25" weights} < /aggregateAssociation/aggregateAssociation.R
	}

	meta {
		author: "Tim Majarian"
		email: "tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "manninglab/aggregateassociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File assoc = "${label}.assoc.RData"
		File log = "aggAssocTest.log"
	}
}

task summary {
	String label
	Array[File] assoc

	Int memory
	Int disk

	command {
		echo "Input files" > summary.log
		echo "label: ${label}" >> summary.log
		echo "assoc: ${sep = ',' assoc}" >> summary.log
		echo "memory: ${memory}" >> summary.log
		echo "disk: ${disk}" >> summary.log
		echo "" >> summary.log
		dstat -c -d -m --nocolor 10 1>>summary.log &
		R --vanilla --args ${label} ${sep = ',' assoc} < /aggregateAssociation/aggregateSummary.R
	}
	
	runtime {
		docker: "manninglab/aggregateassociation:latest"
  	    disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File mhplot = "${label}_association_plots.png"
		File assoc_res = "${label}.groupAssoc.csv"
		File log = "summary.log"
	}
}

workflow group_assoc_wf {
	# fitNull inputs
	File? this_phenotype_file
	String? this_outcome_name
	String? this_outcome_type
	String? this_covariates_string
	File? this_sample_file
	String this_label
	File? this_kinship_matrix
	String? this_id_col

	# aggAssocTest inputs
	Array[File] these_gds_files
	File? this_null_file
	Array[File] these_group_files
	String? this_test
	String? this_pval
	String? this_weights
	
	# other inputs
	Int this_memory
	Int this_disk
	
	# gds and group files must be in the same order, one group file per gds file
	Array[Pair[File,File]] these_gds_groups = zip(these_gds_files, these_group_files)

	File null_genotype_file = these_gds_files[0]
	Boolean have_null = defined(this_null_file)

	if (!have_null) {

		call fitNull {
				input: genotype_file = null_genotype_file, phenotype_file = this_phenotype_file, outcome_name = this_outcome_name, outcome_type = this_outcome_type, covariates_string = this_covariates_string, sample_file = this_sample_file, label = this_label, kinship_matrix = this_kinship_matrix, id_col = this_id_col, memory = this_memory, disk = this_disk
			}

		scatter(this_gds_group in these_gds_groups) {
			
			call aggAssocTest {
				input: gds_file = this_gds_group.left, null_file = fitNull.model, group_file = this_gds_group.right, label = this_label, test = this_test, pval = this_pval, weights = this_weights, memory = this_memory, disk = this_disk
			}
		}
	

		call summary {
			input: assoc = aggAssocTest.assoc, label = this_label, memory = this_memory, disk = this_disk
		}
	}

	if (have_null) {

		scatter(this_gds_group in these_gds_groups) {
			
			call aggAssocTest as aggAssocTest_null_in {
				input: gds_file = this_gds_group.left, null_file = this_null_file, group_file = this_gds_group.right, label = this_label, test = this_test, pval = this_pval, weights = this_weights, memory = this_memory, disk = this_disk
			}
		}
	

		call summary as summary_null_in {
			input: assoc = aggAssocTest.assoc, label = this_label, memory = this_memory, disk = this_disk
		}
	}
}