#module load R
#module load BCFTOOLS
#module load HTSLIB

# Read the list of individuals (i.e. columns in vcf file, which we stored in a text file separately)
individuals <- readLines("individuals_ids.txt")

results <- data.frame(Individual = character(), TotalDiscrepancies = numeric(), stringsAsFactors = FALSE)

common_data_path = "Common/Path/For/VCF/Files"

First_file_name = paste("First_vcf_file",sep="")
other_file_name = paste("Second_VCF_File",sep="")
name = "First_vs_Second_File"
output_file_name = paste("First_vs_Second",sep="")
joint_stats_file = "joined_stats_First_vs_Second"
sample_1_file_name = paste(common_data_path,"sample1.vcf",sep="")
sample_2_file_name = paste(common_data_path,"sample2.vcf",sep="")


system(paste("bgzip -c ",common_data_path,"first_vcf/", First_file_name, ".vcf > ", common_data_path, "first_vcf/", First_file_name, ".vcf.gz",sep=""))

system(paste("bgzip -c ",common_data_path,"second_vcf/", other_file_name, ".vcf > ",common_data_path,"second_vcf/", other_file_name, ".vcf.gz",sep=""))

system(paste("tabix -p vcf ",common_data_path,"first_vcf/", First_file_name, ".vcf.gz",sep=""))

system(paste("tabix -p vcf ",common_data_path,"second_vcf/", other_file_name, ".vcf.gz",sep=""))

counter = 1
# Loop through each individual
for (individual in individuals) {
	print(paste("Running ",name,", Individual: ",individual, ", Running: ",counter,"/",length(individuals),sep=""))
	counter = counter + 1

	temp_command_1 = paste("bcftools view -Oz -c1 -s ",individual," ",common_data_path,"first_vcf/", First_file_name, ".vcf.gz > ", sample_1_file_name, ".gz",sep="")
	system(temp_command_1)

	temp_command_2 = paste("bcftools view -Oz -c1 -s ",individual," ",common_data_path,"second_vcf/", other_file_name, ".vcf.gz > ", sample_2_file_name, ".gz",sep="")
	system(temp_command_2)

	system(paste("tabix -p vcf ", sample_1_file_name, ".gz",sep=""))

	system(paste("tabix -p vcf ", sample_2_file_name, ".gz",sep=""))

	system(paste("bcftools stats ", sample_1_file_name, ".gz ", sample_2_file_name, ".gz > ", common_data_path, joint_stats_file, ".txt",sep=""))

	system(paste("rm ", sample_1_file_name, ".gz ", sample_2_file_name, ".gz ", sample_1_file_name, ".gz.tbi ", sample_2_file_name, ".gz.tbi",sep=""))

	system(paste("awk '
		/^SN.*number of SNPs/ {snp[++s]=$6}
		/^SN.*number of indels/ {indel[++i]=$6}
		END {
			snp_disc = snp[3] - (snp[1] + snp[2]);
			print \"SNP discrepancies:\", snp[3], \" - (\", snp[1], \"+\", snp[2],\") = \",snp_disc;
			indel_disc = indel[3] - (indel[1] + indel[2]);
			print \"Indel discrepancies:\", indel[3], \" - (\", indel[1], \"+\", indel[2],\") = \",indel_disc;
			total_disc = snp_disc + (indel_disc < 0 ? -indel_disc : indel_disc);
			print \"Total discrepancies:\", snp_disc, \"+\",  indel_disc, \"=\", total_disc
		}
	' ", common_data_path, joint_stats_file, ".txt",sep=""))
	
	total_disc <- as.numeric(system(paste("awk '
    /^SN.*number of SNPs/ {snp[++s]=$6}
    /^SN.*number of indels/ {indel[++i]=$6}
    END {
        snp_disc = snp[3] - (snp[1] + snp[2]);
        indel_disc = indel[3] - (indel[1] + indel[2]);
        total_disc = snp_disc + (indel_disc < 0 ? -indel_disc : indel_disc);
        print total_disc
    }
' ", common_data_path, joint_stats_file, ".txt",sep=""), intern = TRUE))


	results <- rbind(results, data.frame(Individual = individual, TotalDiscrepancies = total_disc))

}

write.table(results, file = paste(common_data_path,output_file_name,".txt",sep=""), row.names = FALSE, quote = FALSE, sep = "\t")