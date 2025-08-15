# author: Libby Natola
# date: 30 May 2025
# script to test converting NCBI taxonomy to WoRMS taxonomy for OBIS upload. Test files are from False Creek bioblitz
# Minor updates made by Tim van der Stap, July 4, 2025, primarily to add Darwin Core terms. 

# Install requirements for microViz and phyloseq:
install.packages("microViz", repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos")))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

# load libraries ----
library(tidyverse)
library(taxize) # c 0.9.100
library(microViz)
library(data.table)
library(phyloseq)
library(here)

# data processing ----
## 12S should convert more easily as it's mostly fish and WoRMS is designed for marine taxa, also fish taxonomy has greater consensus than many taxa we sequence with COI 
taxtab.12S <- fread(here("data", "metabarcoding", "12S", "FC_taxonomy_12s.csv"), sep=",", header=T, data.table=FALSE)
taxtab.12S <- taxtab.12S %>% 
  column_to_rownames("ASV")

# convert to taxonomyTable class
tax_tab_12S <- tax_table(as.matrix(taxtab.12S))

# clean up "unknown" values and propagate lowest taxonomic rank out
# can use this for the shiny app if you want the interactive version
# tax_fix_interactive(tax_tab_12S, app_options = list(launch.browser = TRUE))
tax_tab_12S_fixed <- tax_tab_12S %>% 
  tax_fix(
    min_length = 4,
    unknowns = c("unknown"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified") %>%
  as.data.frame()

# Store the verbatimIdentification separately that we'll append to the 
# taxonomic data later:
verbatimIdentification <- tax_tab_12S_fixed$species
names(verbatimIdentification) <- rownames(tax_tab_12S_fixed)
verbatimIdentification <- as.data.frame(verbatimIdentification)

# get list of taxa, remove the taxon ranks (eg "Pleuronectidae family") from the species names
ncbi_taxa_12S <- as.data.frame(tax_tab_12S_fixed) %>%
  transmute(species_clean = str_remove(species, " sp\\.$| family$")) %>%
  pull(species_clean)

# keep the ASV names 
names(ncbi_taxa_12S) <- rownames(tax_tab_12S_fixed)

### convert to WoRMS ----
# get the unique taxa so we don't have to look up the same species in taxize multiple times
tax_lookup_tbl_12S <- tibble(
  ASV = names(ncbi_taxa_12S),
  Species = ncbi_taxa_12S)

unique_species_12S <- unique(tax_lookup_tbl_12S$Species)

# create lookup table (cache results by species)
get_worms_lineage <- function(tax) {
  result <- tryCatch(classification(tax, db = "worms", ), error = function(e) NULL)
  worms_data <- result[[1]]
  
  if (!is.data.frame(worms_data)) {
    return(tibble(Species = tax))  # fallback: just return species name
  } else {
    # convert worms_data into wide form
    lineage <- as_tibble(t(worms_data$name), .name_repair = "minimal")
    colnames(lineage) <- worms_data$rank
    lineage$Species <- tax  # retain for joining
    
    # get id number from worms, format it for AphiaID/LSID
    aphia_num <- last(worms_data$id)
    AphiaID <- paste0("urn:lsid:marinespecies.org:taxname:", aphia_num)
    lineage$AphiaID <- AphiaID
    
    # get the rank of the lowest taxon
    rank <- last(worms_data$rank)
    lineage$Rank <- rank
    
    return(lineage)
  }
}

# map over unique species names
# query worms function - babysit this and manually input prompt responses - I generally choose the first one with status == "accepted"
worms_taxonomy_tbl_12S <- map_dfr(unique_species_12S, get_worms_lineage)

# join lineage info back to ASVs
worms_df_12S <- tax_lookup_tbl_12S %>%
  left_join(worms_taxonomy_tbl_12S, by = "Species") %>% 
  column_to_rownames("ASV")

# specify the order for the columns so the taxonomic ranks are in biological order
cols_order_standard <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Add in the verbatimIdentification, merging based on row names (ASV)
worms_tax_tab_12S_full <- merge(worms_df_12S, verbatimIdentification, by = "row.names", all = TRUE)
rownames(worms_tax_tab_12S_full) <- worms_tax_tab_12S_full$Row.names
worms_tax_tab_12S_full$Row.names <- NULL

worms_tax_tab_12S_full <- worms_tax_tab_12S_full %>%
  rename(scientificNameID = AphiaID,
         taxonRank = Rank)

# For one taxa, Ptychocheilus, taxize could not find taxonomic information. 
# Use the worrms package to look up this information.
ptychocheilus <- worrms::wm_records_name("Ptychocheilus", marine_only = F) %>% filter(AphiaID == "270659")
aphia_id <- ptychocheilus$AphiaID[1]
classification <- worrms::wm_classification(aphia_id)

taxonomy <- classification[classification$rank %in% cols_order_standard, c("rank", "scientificname")] %>%
  pivot_wider(names_from = rank, values_from = scientificname) %>%
  mutate(scientificName = Genus, 
         scientificNameID = paste0("urn:lsid:marinespecies.org:taxname:", ptychocheilus$AphiaID),
         taxonRank = ptychocheilus$rank)

# Add into data frame, however the number of cols differs between data frames, and we don't want
# to duplicate the taxonomy df cols across the row. Find the common cols:
common_cols <- intersect(names(worms_tax_tab_12S_full), names(taxonomy))

# Find matching row index:
row_idx <- which(worms_tax_tab_12S_full$Species == "Ptychocheilus")

# In that row, replace the values where columns match:
worms_tax_tab_12S_full[row_idx, common_cols] <- taxonomy[1, common_cols]

# Select the right columns, renaming slightly:
worms_tax_tab_12S_full <- worms_tax_tab_12S_full %>%
  select(kingdom = Kingdom,
         phylum = Phylum,
         class = Class,
         order = Order,
         family = Family,
         genus = Genus,
         scientificName = Species,
         scientificNameID,
         taxonRank,
         verbatimIdentification)

# Grab authority data from the AphiaID
worms_tax_tab_12S_full$aphiaID <- as.numeric(sub("urn:lsid:marinespecies.org:taxname:", "", worms_tax_tab_12S_full$scientificNameID))

get_authority <- function(id) {
  tryCatch({
    rec <- worrms::wm_record(id = id)
    rec$authority
  }, error = function(e) NA)
}

worms_tax_tab_12S_full <- worms_tax_tab_12S_full %>%
  mutate(scientificNameAuthorship = sapply(worms_tax_tab_12S_full$aphiaID, get_authority)) %>%
  select(-aphiaID)

# Change rownames to column, remove NAs"
worms_tax_tab_12S_full <- rownames_to_column(worms_tax_tab_12S_full, var = "ASV")
worms_tax_tab_12S_full <- worms_tax_tab_12S_full %>%
  mutate(across(everything(), ~ replace(.x, is.na(.x), "")))

# Save the file:
write.csv(worms_tax_tab_12S_full, here("data", "metabarcoding", "12S", "FC_taxonomy_12s_worms.csv"), quote = FALSE, row.names = TRUE)
