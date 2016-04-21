### mixcr analysis

# Function to compile and format IMGT results using R
compile_imgt_output <- function(summary_file, output_dir, project) {
    # Build and use Unix commands
    imgt_file <- file.path(output_dir, 
                           paste(project, "compiled_imgt_output.txt", 
                                 sep = "_"))
    
    if (!file.exists(imgt_file)) {
        imgt_summary <- read.delim(summary_file, stringsAsFactors = FALSE)
        
        imgt_compiled <- imgt_summary %>% 
            filter(Functionality == "productive") %>% 
            select(Sequence.ID, V.GENE.and.allele, J.GENE.and.allele, 
                   AA.JUNCTION) %>% 
            transmute(lib_id = str_extract(Sequence.ID, "lib|SRR[0-9]+"),
                      v_gene = str_extract(V.GENE.and.allele, 
                                           "TR.*?(?=(\\*))"),
                      j_gene = str_extract(J.GENE.and.allele, 
                                           "TR.*?(?=(\\*))"),
                      junction = AA.JUNCTION)
        
        write.table(imgt_compiled, imgt_file, 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    return(imgt_file)
}

# Function to format IMGT results using R
format_imgt_jxns <- function(imgt_file) {
    # Read in IMGT junctions
    imgt_jxns <- read_delim(imgt_file, delim = "\t")
    
    if ("libID" %in% names(imgt_jxns)) {
        imgt_jxns <- imgt_jxns %>% 
            rename(lib_id = libID,
                   v_gene = `V-GENE`,
                   j_gene = `J-GENE`,
                   junction = JUNCTION)
    }
    
    return(imgt_jxns)
}

# Function to filter IMGT junctions
filter_imgt_jxns <- function(imgt_jxns, min_length = 6) {
    
    imgt_jxns <- imgt_jxns %>% 
        filter(str_detect(v_gene, "^((?![C-G]).)*$"),
               str_detect(j_gene, "^((?![C-G]).)*$"),
               str_detect(junction, "^C"),
               str_detect(junction, "^((?!(\\*|_)).)*$"),
               !duplicated(.[, c(1:4)]),
               str_length(junction) > min_length)
    
    return(imgt_jxns)
}

# Function to combine individual * mixcr_clns.txt files using unix: 1) add header
# from one file to newfile. 2) grep contents of all files and 3) append them to
# newfile
combine_mixcr_outputs <- function(mixcr_dir, output_dir, project) {

    # Build and use Unix commands
    mixcr_combined_file <- file.path(output_dir, 
                                     paste(project, "compiled_mixcr_output.txt", 
                                           sep = "_"))
    
    if (!file.exists(mixcr_combined_file)) {
        
        # Select one file to grab the column headers
        mixcr_tmp_file <- data_frame(
            file = list.files(mixcr_dir, full.names = TRUE)) %>% 
            mutate(size = file.size(file)) %>% 
            filter(str_detect(str_to_lower(file), "clns.txt"),
                   size > 0) %>% # make sure not to select empty files
            slice(1) %>% 
            select(file) %>% 
            unlist()
        
        header_cmd <- sprintf("head -1 %s > %s", mixcr_tmp_file, mixcr_combined_file)
        system(header_cmd)
        
        compile_cmd <- sprintf("grep '*' %s/*lns.txt >> %s",
                               mixcr_dir, mixcr_combined_file)
        system(compile_cmd)
    }
    
    return(mixcr_combined_file)
}

# Function to read in MiXCR data then extract & format variables
format_mixcr_jxns_old <- function(mixcr_combined_file) {
    # Read in MiXCR clones
    mixcr_jxns <- read.delim(mixcr_combined_file)
    print(nrow(mixcr_jxns))
    
    print(head(mixcr_jxns))
    # Extract & format key variables
    mixcr_jxns <- mixcr_jxns %>% 
        transmute(lib_id = as.character(str_match(Clone.count, "(lib|SRR)[0-9]+")),
                  cln_count = as.numeric(str_match(Clone.count, "(?<=:)[0-9]+")),
                  v_gene = str_extract(All.V.hits,
                                      "TR[A-Z]+[0-9]*(\\-[0-9])*(DV[0-9]+)*"),
                  v_gene_score = str_extract(All.V.hits, "(?<=\\()[0-9]+") %>% 
                      as.numeric(),
                  j_gene = str_extract(All.J.hits, 
                                      "TR[A-Z]+[0-9]*(\\-[0-9][A-Z]*)*"),
                  j_gene_score = str_extract(All.J.hits, "(?<=\\()[0-9]+") %>% 
                      as.numeric(),
                  junction = as.character(AA..seq..CDR3))
    
    return(mixcr_jxns)
}

format_mixcr_jxns <- function(mixcr_combined_file) {
    # Read in MiXCR clones
    mixcr_jxns <- read_delim(mixcr_combined_file, delim = "\t")

    # Extract & format key variables
    mixcr_jxns <- mixcr_jxns %>% 
        transmute(lib_id = as.character(str_extract(`Clone count`, "(lib|SRR)[0-9]+")),
                  cln_count = as.numeric(str_extract(`Clone count`, "(?<=:)[0-9]+")),
                  v_gene = str_extract(`All V hits`,
                                       "TR[A-Z]+[0-9]*(\\-[0-9])*(DV[0-9]+)*"),
                  v_gene_score = str_extract(`All V hits`, "(?<=\\()[0-9]+") %>% 
                      as.numeric(),
                  v_align = `All V alignment`,
                  j_gene = str_extract(`All J hits`, 
                                       "TR[A-Z]+[0-9]*(\\-[0-9][A-Z]*)*"),
                  j_align = `All J alignment`,
                  j_gene_score = str_extract(`All J hits`, "(?<=\\()[0-9]+") %>% 
                      as.numeric(),
                  junction = as.character(`AA. seq. CDR3`)) %>% 
        rowwise() %>% 
        mutate(v_region_nt_overlap = str_split(v_align, pattern = "\\|") %>% 
                   unlist() %>% 
                   .[length(.) - 2],
               v_region_align_score = str_split(v_align, pattern = "\\|") %>% 
                   unlist() %>% 
                   .[length(.)],
               j_region_nt_overlap = str_split(j_align, pattern = "\\|") %>% 
                   unlist() %>% 
                   .[length(.) - 2],
               j_region_align_score = str_split(j_align, pattern = "\\|") %>% 
                   unlist() %>% 
                   .[length(.)]) %>% 
        select(one_of(c("lib_id", "cln_count", "v_gene", "v_gene_score",
                        "v_align", "v_region_nt_overlap",
                        "v_region_align_score",
                        "j_gene", "j_gene_score",
                        "j_align", "j_region_nt_overlap",
                        "j_region_align_score",
                        "junction")))
    
    return(mixcr_jxns)
}

# Function to filter MiXCR junctions
filter_mixcr_jxns <- function(mixcr_jxns, min_count = 1, min_length = 6) {
    
    mixcr_jxns <- mixcr_jxns %>% 
        filter(str_detect(v_gene, "^((?![C-G]).)*$"), # TRA/B only
               str_detect(j_gene, "^((?![C-G]).)*$"), # TRA/B only
               str_extract(v_gene, "(?<=TR)[A-Z]") == 
                   str_extract(j_gene, "(?<=TR)[A-Z]"), # segments must match
               str_detect(junction, "^C"), # jxn AA cannot start with C
               str_detect(junction, "^((?!(\\*|_)).)*$"), # functional only
               !duplicated(.[, c(1, 3:5)]), # remove dups
               cln_count >= min_count,
               str_length(junction) > min_length)
    
    return(mixcr_jxns)
}

# There should be only one alpha and one beta junction for each lib; select
# the top hit for each lib, sorting first by clone count and second by
# alignment _score
select_top_jxns <- function(mixcr_jxns, allow_multi = FALSE, detail = FALSE) {

   mixcr_jxns <- mixcr_jxns %>% 
        mutate(a_b = ifelse(str_detect(v_gene, "TRA"), "TRA", "TRB")) %>% 
        group_by(lib_id, a_b) %>% 
        arrange(desc(cln_count), desc(v_gene_score)) 
   
   if (allow_multi) {
       mixcr_jxns <- mixcr_jxns %>% 
           slice(1:2) %>% 
           ungroup() %>% 
           group_by(lib_id, a_b) %>% 
           mutate(count_diff = c(0, diff(cln_count)) / lag(cln_count),
                  score_diff = c(0, diff(v_gene_score)) / lag(v_gene_score),
                  best = ifelse(is.na(count_diff), TRUE, FALSE),
                  second = ifelse(a_b == "TRA" & !best &
                                      (count_diff > -0.5 | score_diff > -0.1),
                                  TRUE, FALSE)) %>% 
           ungroup() %>% 
           filter(best | second)
   } else {
       mixcr_jxns <- mixcr_jxns %>% 
           slice(1) %>% 
           ungroup()
   }

   if (detail) {
       mixcr_jxns %>% 
           select(-a_b, -count_diff, -score_diff, -best, -second)
   } else {
       mixcr_jxns %>% 
           select(lib_id, v_gene, j_gene, junction)
   }

}
