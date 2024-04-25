#install.packages("ieugwasr")
#install.packages("dplyr")
library(dplyr)
library(ieugwasr)

for (chrom in 1:22) {
    file_path <- paste("../data/snp/by_chromosome/chromosome", chrom,".csv", sep = "")
    
    df <- read.csv(file_path, header = FALSE)

    ld_matrix_data_frame <- as.data.frame(ld_matrix(df[,1]))

    ld_snps <- colnames(ld_matrix_data_frame)
    ld_snps <- substr(ld_snps, 1, nchar(ld_snps) - 4)
    filtered_df <- subset(df, df[,1] %in% ld_snps)

    for(i in 1:nrow(filtered_df)) {
        my_rsid <- filtered_df[i, 1]
        my_bp <- filtered_df[i,2]
        rsids <- system(paste("python3 find_near_rsids.py", my_rsid, chrom, my_bp), intern = TRUE)
        if(length(rsids) > 1) {
            rsid_vector <- readLines(textConnection(rsids))
            my_snp <- rsid_vector[1]
            
            check_ld_ref <- function(rsid, pop = "EUR", opengwas_jwt = get_opengwas_jwt()) {
                ld_ref <- ld_reflookup(rsid, pop, opengwas_jwt)
                return(!is.null(ld_ref))
            }

            print(my_snp)
            
            if (check_ld_ref(my_snp)) {    
                ld_matrix_data_frame <- as.data.frame(ld_matrix(unlist(strsplit(rsid_vector, split = "\n"))))
                column_index <- grep(paste0("^", my_snp), colnames(ld_matrix_data_frame))
                final_matrix_for_my_snp <- ld_matrix_data_frame[column_index, ]
                file <- paste("../data/matrixes/chrom", chrom, "/", my_snp,".csv", sep = "")
                write.csv(final_matrix_for_my_snp, file)
                
                if (ncol(final_matrix_for_my_snp) > 0) {
                    final_matrix_for_my_snp <- final_matrix_for_my_snp[, which(final_matrix_for_my_snp[1, ] >= 0.75)]
                    print(final_matrix_for_my_snp)
                } else {
                    print("Final matrix has no columns.")
                
                }
            }
        }
    }
}

