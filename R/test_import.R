
# load packages -----------------------------------------------------------

library(readr)
library(stringr)
library(dplyr)


# load functions & set params ---------------------------------------------

source("R/import_data.R")
sample_regex = "(lib|SRR)[0-9]+"


# test IMGT read functions -------------------------------------------------

# for single file...
imgt_file <- "/Users/jaeddy/data/analyses/ilc_tcr_analysis/imgt/ilc_tcr_1.txz"

extract_cmd <- sprintf("tar -xf '%s' -O '1_Summary.txt'", imgt_file)
imgt_summary <- pipe(extract_cmd) %>% 
    read_file() %>% 
    read_imgt_summary()

# test IMGT parse functions -----------------------------------------------

imgt_df <- parse_imgt_summary(imgt_summary)

# test MiXCR read functions ------------------------------------------------

mixcr_file <- "/Users/jaeddy/data/analyses/ilc_tcr_analysis/mixcr/SRR2088075_mixcrClns.txt"

mixcr_clones <- read_mixcr_clones(mixcr_file)


# test MiXCR parse functions ----------------------------------------------

mixcr_df <- parse_mixcr_clones(mixcr_clones)
