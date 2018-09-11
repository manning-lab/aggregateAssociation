# makeGroupFile.wdl
task makeGroups {
	File genotype_file
	File region_file
	File? variant_file
	Float? max_maf
	Float? min_maf
	String out_pref
	Int index
	
	Int memory
	Int disk

	command {
		R --vanilla --args ${genotype_file} ${region_file} ${default="NA" variant_file} ${default="NA" max_maf} ${default="NA" min_maf} ${out_pref} ${index} < /aggregateAssociation/makeGroupFile.R
	}

	runtime {
		docker: "manninglab/aggregateassociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File group_file = "${out_pref}.${index}.csv"
	}
}

task combineGroups {
	Array[File] these_groups
	String out_pref

	command <<<
		echo "args <- commandArgs(trailingOnly=T)" > script.R && \
		echo "g.files <- unlist(strsplit(args[1],','))" >> script.R && \
		echo "out.pref <- args[2]" >> script.R && \
		echo "library(data.table)" >> script.R && \
		echo "all.g <- list()" >> script.R && \
		echo "for (g in g.files){" >> script.R && \
		echo "all.g[[length(all.g)+1]] <- fread(g, data.table = F, stringsAsFactors = F)}" >> script.R && \
		echo "all.g <- do.call(rbind,all.g)" >> script.R && \
		echo "write.table(all.g, file=paste0(out.pref,'.csv'), quote = F, sep = ',', row.names = F)" >> script.R && \
		R --vanilla --args ${sep="," these_groups} ${out_pref} < script.R
	>>>

	runtime {
		docker: "manninglab/aggregateassociation:latest"
		disks: "local-disk 20 SSD"
		memory: "7 GB"
	}

	output {
		File final_groups = "${out_pref}.csv"
	}
}

workflow makeGroupFile {
	# makeGroups inputs
	Array[File] these_genotype_files
	File this_region_file
	File? this_variant_file
	Float? this_max_maf
	Float? this_min_maf
	String this_out_pref
	
	# other inputs
	Int this_memory
	Int this_disk

	Array[Int] these_genotype_indices = range(length(these_genotype_files))
	
	scatter(this_genotype_index in these_genotype_indices) {
		call makeGroups {
			input: genotype_file = these_genotype_files[this_genotype_index], region_file = this_region_file, variant_file = this_variant_file, max_maf = this_max_maf, min_maf = this_min_maf, out_pref = this_out_pref, index = this_genotype_index, memory = this_memory, disk = this_disk
		}
	}

	call combineGroups {
		input: these_groups = makeGroups.group_file, out_pref = this_out_pref
	}

	output {
		File groups = combineGroups.final_groups
	}
}