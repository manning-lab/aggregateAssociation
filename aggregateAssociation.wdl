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
		echo "Input files" > fitNull_out.log
		echo "genotype_file: ${genotype_file}" >> fitNull_out.log
		echo "phenotype_file: ${phenotype_file}" >> fitNull_out.log
		echo "outcome_name: ${outcome_name}" >> fitNull_out.log
		echo "outcome_type: ${outcome_type}" >> fitNull_out.log
		echo "covariates_string: ${covariates_string}" >> fitNull_out.log
		echo "sample_file: ${sample_file}" >> fitNull_out.log
		echo "label: ${label}" >> fitNull_out.log
		echo "kinship_matrix: ${kinship_matrix}" >> fitNull_out.log
		echo "id_col: ${id_col}" >> fitNull_out.log
		echo "memory: ${memory}" >> fitNull_out.log
		echo "disk: ${disk}" >> fitNull_out.log
		echo "" >> fitNull_out.log
		dstat -c -d -m --nocolor 10 1>>fitNull_out.log &
		R --vanilla --args ${genotype_file} ${phenotype_file} ${outcome_name} ${outcome_type} ${default="NA" covariates_string} "NA" "NA" "NA" ${sample_file} ${label} ${kinship_matrix} ${id_col} < /aggregateAssociation/genesis_nullmodel.R
	}

	runtime {
		docker: "manninglab/aggregateassociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File model = "${label}_null.RDa"
		File log = "fitNull_out.log"
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
		echo "Input files" > aggAssocTest_out.log
		echo "gds_file: ${gds_file}" >> aggAssocTest_out.log
		echo "null_file: ${null_file}" >> aggAssocTest_out.log
		echo "group_file: ${group_file}" >> aggAssocTest_out.log
		echo "label: ${label}" >> aggAssocTest_out.log
		echo "test: ${test}" >> aggAssocTest_out.log
		echo "pval: ${pval}" >> aggAssocTest_out.log
		echo "weights: ${weights}" >> aggAssocTest_out.log
		echo "memory: ${memory}" >> aggAssocTest_out.log
		echo "disk: ${disk}" >> aggAssocTest_out.log
		echo "" >> aggAssocTest_out.log
		dstat -c -d -m --nocolor 10 1>>aggAssocTest_out.log &
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
		File log = "aggAssocTest_out.log"
	}
}

task summary {
	String label
	Array[File] assoc

	Int memory
	Int disk

	command {
		echo "Input files" > summary_out.log
		echo "label: ${label}" >> summary_out.log
		echo "assoc: ${sep = ',' assoc}" >> summary_out.log
		echo "memory: ${memory}" >> summary_out.log
		echo "disk: ${disk}" >> summary_out.log
		echo "" >> summary_out.log
		dstat -c -d -m --nocolor 10 1>>summary_out.log &
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
		File log = "summary_out.log"
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
	File this_group_file
	String? this_test
	String? this_pval
	String? this_weights
	
	# other inputs
	Int this_memory
	Int this_disk
	
	File null_genotype_file = these_gds_files[0]
	Boolean have_null = defined(this_null_file)

	if (!have_null) {

		call fitNull {
				input: genotype_file = null_genotype_file, phenotype_file = this_phenotype_file, outcome_name = this_outcome_name, outcome_type = this_outcome_type, covariates_string = this_covariates_string, sample_file = this_sample_file, label = this_label, kinship_matrix = this_kinship_matrix, id_col = this_id_col, memory = this_memory, disk = this_disk
			}

		scatter(this_gds_file in these_gds_files) {
			
			call aggAssocTest {
				input: gds_file = this_gds_file, null_file = fitNull.model, group_file = this_group_file, label = this_label, test = this_test, pval = this_pval, weights = this_weights, memory = this_memory, disk = this_disk
			}
		}
	

		call summary {
			input: assoc = aggAssocTest.assoc, label = this_label, memory = this_memory, disk = this_disk
		}
	}

	if (have_null) {

		scatter(this_gds_file in these_gds_files) {
			
			call aggAssocTest as aggAssocTest_null_in {
				input: gds_file = this_gds_file, null_file = this_null_file, group_file = this_group_file, label = this_label, test = this_test, pval = this_pval, weights = this_weights, memory = this_memory, disk = this_disk
			}
		}
	

		call summary as summary_null_in {
			input: assoc = aggAssocTest_null_in.assoc, label = this_label, memory = this_memory, disk = this_disk
		}
	}
}