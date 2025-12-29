# author: Libby Natola
# date: 30 May 2025
# script to test converting NCBI taxonomy to WoRMS taxonomy for OBIS upload. Test files are from False Creek bioblitz
# Minor updates made by Tim van der Stap, July 4, 2025, primarly to add some Darwin Core terms. 

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


taxtab.16S <- fread(here("data", "metabarcoding", "16S", "FC_taxonomy_16s.csv"), sep=",", header=T, data.table=FALSE)
tax_tab_16S <- as.data.frame(tax_table(as.matrix(taxtab.16S)))

# Where the lowest taxonomic information contains uncultured or unidentified, make these blank cells
# Save the lowest taxonomic resolution in a column 'verbatimIdentification'
tax_tab_16S <- tax_tab_16S %>%
  mutate(across(everything(), ~ ifelse(grepl("^uncultured_|^unidentified_", .), NA, .)),
         across(everything(), ~ ifelse(grepl("marine_metagenome|bacterium_enrichment|soda_lake|metagenome|hydrothermal_vent|seawater_metagenome", .), NA, .)))

# Convert to taxonomyTable class
tax_tab_16S <- tax_tab_16S %>%
  mutate(across(everything(), ~na_if(., "uncultured")),
         across(everything(), ~na_if(., "Unknown_Family")))

# Record verbatimIdentification
tax_tab_16S$verbatimIdentification <- apply(tax_tab_16S, 1, function(row) {
  # Remove NA and blank values
  valid_values <- row[!is.na(row) & row != ""]
  if (length(valid_values) > 0) {
    tail(valid_values, 1)
  } else {
    NA
  }
})

# Lots of 'wonky' nomenclature, let's clean up from phylum moving down:
tax_tab_16S <- tax_tab_16S %>%
  mutate(across(.cols = -verbatimIdentification, .fns = ~ ifelse(grepl("NB1-j|Sva0485|Marinimicrobia_\\(SAR406_clade\\)|SAR324_clade\\(Marine_group_B\\)", .), "", .)))

# Class
tax_tab_16S <- tax_tab_16S %>%
  mutate(across(.cols = -verbatimIdentification, .fns = ~ ifelse(grepl("BD2-11_terrestrial_group|OM190|Pla3_lineage", .), "", .)))

# Order
tax_tab_16S <- tax_tab_16S %>%
  mutate(across(.cols = -verbatimIdentification, .fns = ~ ifelse(grepl("B2M28|Arctic97B-4_marine_group|PeM15|JGI_0000069-P22|Bacteroidetes_VC2\\.1_Bac22|PB19|Blfdi19|SAR11_clade|SAR202_clade|Ga0077536|UBA10353_marine_group|Milano-WF1B-44|Marine_Group_II|HOC36", .), "", .)))

# Family
family_pattern <- paste("SAR116_clade", "AEGEAN-169_marine_group", "Clade_II", "NS9_marine_group", 
                        "DEV007", "KI89A_clade", "S25-593", "Sva1033", "Bacteroidetes_BD2-2", "PS1_clade", 
                        "Clade_III", "SAR86_clade", "OM182_clade", "Clade_IV", "OCS116_clade", "NS11-12_marine_group", "NS7_marine_group", "Clade_I",
                        sep = "|")
tax_tab_16S <- tax_tab_16S %>%
  mutate(across(.cols = -verbatimIdentification, .fns = ~ ifelse(grepl(family_pattern, .), "", .)))

tax_tab_16S$family <- gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax_tab_16S$family)

# Genus
genus_pattern <- paste("OM27_clade", "OM60\\(NOR5\\)_clade", "SAR92_clade", "C1-B045", "Roseobacter_clade_NAC11-7_lineage", "NS10_marine_group", "Cyanobium_PCC-6307",
                       "OM75_clade", "\\[Desulfobacterium\\]_catecholicum_group", "MD3-55", "RS62_marine_group", "SEEP-SRB1", "JL-ETNP-F27", "endosymbionts", 
                       "\\[Caedibacter\\]_taeniospiralis_group", "CL500-3", "BD1-7_clade", "P3OB-42", "HIMB11", "SUP05_cluster", "LS-NOB", "Sva0996_marine_group",
                       "NS4_marine_group", "NS5_marine_group", "Synechococcus_CC9902", "IS-44", "SCGC_AAA164-E04", "Defluviitaleaceae_UCG-011", "LCP-80", "NS2b_marine_group",
                       "Sva0081_sediment_group", "NS3a_marine_group", "OM43_clade", "Atelocyanobacterium_\\(UCYN-A\\)", "ML602J-51", "MB11C04_marine_group", "Prevotella_7", 
                       "hgcI_clade", "Prevotella_9", sep = "|")

tax_tab_16S <- tax_tab_16S %>%
  mutate(across(.cols = -verbatimIdentification, .fns = ~ ifelse(grepl(genus_pattern, .), "", .)))

# species
species_pattern <- paste("bacterium_RCC", "alpha_proteobacterium", "SAR86_cluster", "bacterium_WHC3-3", "bacterium_WHC4-9", sep = "|")

tax_tab_16S <- tax_tab_16S %>%
  mutate(across(.cols = -verbatimIdentification, .fns = ~ ifelse(grepl(species_pattern, .), "", .)))

# Replace _ and - with blank space, and make sure that only the first word is capitalized
tax_tab_16S <- tax_tab_16S %>%
  mutate(across(where(is.character) & !any_of("verbatimIdentification"), ~ {
    str_replace_all(., "[-_]", " ") %>%
      str_to_sentence()
  }))

tax_tab_16S[tax_tab_16S == ""] <- NA

verbatimIdentification <- tax_tab_16S %>% select(ASV, verbatimIdentification)

# Save the cleaned, but non-WoRMS-matched, version to GitHub:
write.csv(tax_tab_16S, here("data", "metabarcoding", "16S", "FC_taxonomy_16S_cleaned.csv"))

# Grab the most granular taxonomic information, except verbatimIdentification:
tax_tab_16S <- tax_tab_16S %>%
  rowwise() %>%
  mutate(scientificName = {
    row_data <- c_across(-verbatimIdentification)
    last_value <- tail(na.omit(row_data), 1)
    if (length(last_value) == 0) NA_character_ else last_value
  }) %>%
  ungroup()

unique_species_16S <- unique(tax_tab_16S$scientificName)

# map over unique species names
# query worms function - babysit this and manually input prompt responses - I generally choose the first one with status == "accepted"
worms_taxonomy_tbl_16S <- map_dfr(unique_species_16S, get_worms_lineage)

worms_taxonomy_tbl_16S <- worms_taxonomy_tbl_16S %>%
  rename(scientificName = Species,
         scientificNameID = AphiaID,
         taxonRank = Rank)

# Grab authority data from the AphiaID
worms_taxonomy_tbl_16S <- worms_taxonomy_tbl_16S %>%
  mutate(aphiaID = ifelse(is.na(scientificNameID), NA_character_, gsub("urn:lsid:marinespecies.org:taxname:", "", scientificNameID)))

get_authority <- function(id) {
  tryCatch({
    if (is.na(id) || id == "") return(NA_character_)
    rec <- worrms::wm_record(as.integer(id))
    rec$authority
  }, error = function(e) NA_character_)
}

worms_taxonomy_tbl_16S <- worms_taxonomy_tbl_16S %>%
  mutate(scientificNameAuthorship = sapply(worms_taxonomy_tbl_16S$aphiaID, get_authority))

worms_taxonomy_tbl_16S <- worms_taxonomy_tbl_16S %>%
  select(kingdom = Kingdom, 
         phylum = Phylum, 
         class = Class, 
         order = Order, 
         family = Family, 
         genus = Genus,
         scientificName,
         scientificNameID,
         taxonRank,
         scientificNameAuthorship)

# Merge this in with the original tax_tab_16S df:
occ_16S <- tax_tab_16S %>% select(scientificName, verbatimIdentification, ASV) %>%
  left_join(worms_taxonomy_tbl_16S, by = "scientificName")

occ_16S[is.na(occ_16S)] <- ""

# Save: 
write.csv(occ_16S, here("data", "metabarcoding", "16S", "FC_taxonomy_16S_worms.csv"))
