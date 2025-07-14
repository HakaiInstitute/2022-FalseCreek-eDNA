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

# data processing ----
## 12S ----
## 12S should convert more easily as it's mostly fish and WoRMS is designed for marine taxa, also fish taxonomy has greater consensus than many taxa we sequence with COI 
### load/clean data ----
# sequence table
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

# Store the verbatimIdentification separately that we'll cbind to the 
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
  Species = ncbi_taxa_12S
)

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

# did we lose anyone?
all(rownames(worms_tax_tab_12S_full) %in% rownames(tax_tab_12S)) & all(rownames(tax_tab_12S) %in% rownames(worms_tax_tab_12S_full))

# Change rownames to column, remove NAs"
worms_tax_tab_12S_full <- rownames_to_column(worms_tax_tab_12S_full, var = "ASV")
worms_tax_tab_12S_full <- worms_tax_tab_12S_full %>%
  mutate(across(everything(), ~ replace(.x, is.na(.x), "")))

# Save the file:
write.csv(worms_tax_tab_12S_full, here("data", "metabarcoding", "12S", "FC_taxonomy_12s_worms.csv"), quote = FALSE, row.names = TRUE)

## COI ----
# COI is a marker that amplifies a broader taxonomic group, many of which won't be well annotated within WoRMS

taxtab.COI <- fread(here("data", "metabarcoding", "COI", "FC_taxonomy_COI_Revised_June2025.csv"), sep=",", header=T, data.table=FALSE)
colnames(taxtab.COI)[-1] <- str_to_title(colnames(taxtab.COI)[-1])
taxtab.COI <- taxtab.COI %>% column_to_rownames("ASV")

# convert to taxonomyTable class
tax_tab_COI <- as.data.frame(tax_table(as.matrix(taxtab.COI)))

# make uncultureds unknowns
tax_tab_COI$Species <- gsub("^uncultured.*", "unknown", tax_tab_COI$Species)

# Convert to taxonomyTable class
tax_tab_COI <- tax_tab_COI %>%
  mutate(across(everything(), ~na_if(., "unknown")),
         across(everything(), ~na_if(., "no identification")),
         across(everything(), ~na_if(., "unknown class")),
         across(everything(), ~na_if(., "unknown order")),
         across(everything(), ~na_if(., "unknown family")),
         across(everything(), ~na_if(., "unknown genus")),
         across(everything(), ~na_if(., "uncultured *")))

# Grab the most granular taxonomic data (i.e., the last column populated before
# 'NA' information is provided):
tax_tab_COI$scientificName <- apply(tax_tab_COI, 1, function(row) {
  last_value <- tail(na.omit(row), 1) 
  if (length(last_value) == 0) NA else last_value
})

verbatimIdentification_COI <- tax_tab_COI$scientificName
names(verbatimIdentification_COI) <- rownames(tax_tab_COI)
verbatimIdentification_COI <- as.data.frame(verbatimIdentification_COI)

# Capture identificationQualifiers and remove them from the scientificName:
tax_tab_COI_fixed <- tax_tab_COI %>%
  mutate(identificationQualifier = case_when(
    grepl("cf. promare", Species) ~ "cf. promare",
    grepl("cf. depressum", Species) ~ "cf. depressum",
    grepl("cf. urtica", Species) ~ "cf. urtica",
    grepl("aff. granulata", Species) ~ "cf. aff. granulata"))
tax_tab_COI_fixed$scientificName <- gsub("cf. promare|cf. depressum|cf. urtica|aff. granulata", "", tax_tab_COI_fixed$scientificName) %>% trimws()

# remove these strings
tax_tab_COI_fixed$scientificName <- gsub("sp\\..*", "", tax_tab_COI_fixed$scientificName)
tax_tab_COI_fixed$scientificName <- gsub("CCMP1545|CMC01|EAC-2010|KSL-2016", "", tax_tab_COI_fixed$scientificName) %>% trimws()

# get list of taxa, remove the taxon ranks (eg "Pleuronectidae family") from the species names
ncbi_taxa_COI <- as.data.frame(tax_tab_COI_fixed) %>%
  pull(scientificName)

# keep the ASV names 
names(ncbi_taxa_COI) <- rownames(tax_tab_COI_fixed)

### convert to WoRMS ----
# get the unique taxa so we don't have to look up the same species in taxize multiple times
tax_lookup_tbl_COI <- tibble(
  ASV = names(ncbi_taxa_COI),
  Species = ncbi_taxa_COI
)

unique_species_COI <- unique(tax_lookup_tbl_COI$Species)

# map over unique species names
# query worms function - babysit this and manually input prompt responses - I generally choose the first one with status == "accepted"
## 1   149045  Nitzschia      A.H. Hassall, 1845    accepted
## 1   111022  Amathia        Lamouroux, 1812       accepted
## 2   160576  Chlorella      M.Beijerinck, 1890    accepted
worms_taxonomy_tbl_COI <- map_dfr(unique_species_COI, get_worms_lineage) 

# join lineage info back to ASVs
worms_df_COI <- tax_lookup_tbl_COI %>%
  left_join(worms_taxonomy_tbl_COI, by = "Species") %>% 
  column_to_rownames("ASV") %>%
  rename(scientificNameID = AphiaID,
         taxonRank = Rank)

# Add in the verbatimIdentification, merging based on row names (ASV)
worms_tax_tab_COI_full <- merge(worms_df_COI, verbatimIdentification_COI, by = "row.names", all = TRUE)
rownames(worms_tax_tab_COI_full) <- worms_tax_tab_COI_full$Row.names
worms_tax_tab_COI_full$Row.names <- NULL
worms_tax_tab_COI_full <- worms_tax_tab_COI_full %>% rename(verbatimIdentification = verbatimIdentification_COI)

# Grab authority data from the AphiaID
worms_tax_tab_COI_full$aphiaID <- as.numeric(sub("urn:lsid:marinespecies.org:taxname:", "", worms_tax_tab_COI_full$scientificNameID))

get_authority <- function(id) {
  tryCatch({
    rec <- worrms::wm_record(id = id)
    rec$authority
  }, error = function(e) NA)
}

worms_tax_tab_COI_full <- worms_tax_tab_COI_full %>%
  mutate(scientificNameAuthorship = sapply(worms_tax_tab_COI_full$aphiaID, get_authority))

worms_tax_tab_COI_full <- worms_tax_tab_COI_full %>%
  select(kingdom = Kingdom, 
         phylum = Phylum, 
         class = Class, 
         order = Order, 
         family = Family, 
         genus = Genus,
         scientificName = Species,
         scientificNameID,
         taxonRank,
         verbatimIdentification, aphiaID, scientificNameAuthorship)

# Upon inspection, there are 8 unique taxa for which there is no associated WoRMS LSID. 
# Three of these taxa can be manually found on WoRMS, and their AphiaIDs have been verified with Libby Natola
# Their taxonomic information, lsid, rank and scientificNameAuthorship are created below and then merged
# into the worms_tax_tab_COI_full data. 

names <- c("Nitzschia inconspicua", "Prorocentrum pervagatum", "Synchaeta kitina")
records <- worrms::wm_records_names(name = names, marine_only = F)

# Extract AphiaID, authority, etc.
flat_info <- lapply(seq_along(records), function(i) {
  rec <- records[[i]]
  name <- names[i]
  
  if (!is.null(rec) && nrow(rec) > 0) {
    tibble(
      scientificName = name,
      kingdom = rec$kingdom[1],
      phylum = rec$phylum[1],
      class = rec$class[1],
      order = rec$order[1],
      family = rec$family[1],
      genus = rec$genus[1],
      aphiaID = rec$AphiaID[1],
      scientificNameAuthorship = rec$authority[1],
      taxonRank = rec$rank[1],
    )
  } else {
    tibble(
      scientificName = name,
      scientificNameID = NA,
      kingdom = NA,
      phylum = NA,
      class = NA,
      order = NA,
      family = NA,
      genus = NA,
      aphiaID = NA,
      scientificNameAuthorship = NA,
      taxonRank = NA
    )
  }
})

df <- bind_rows(flat_info) %>%
  mutate(scientificNameID = paste0("urn:lsid:marinespecies.org:taxname:", aphiaID),
         verbatimIdentification = scientificName)
df <- df[c("kingdom", "phylum", "class", "order", "family", "genus", "scientificName",
           "scientificNameID", "taxonRank", "verbatimIdentification", "aphiaID", "scientificNameAuthorship")]

tax_Nitzschia <- df %>% filter(scientificName == "Nitzschia inconspicua")
tax_Proro <- df %>% filter(scientificName == "Prorocentrum pervagatum")
tax_kitina <- df %>% filter(scientificName == "Synchaeta kitina")

# Add into data frame:
worms_tax_tab_COI_full[worms_tax_tab_COI_full$scientificName == "Nitzschia inconspicua", ] <- tax_Nitzschia
worms_tax_tab_COI_full[worms_tax_tab_COI_full$scientificName == "Prorocentrum pervagatum", ] <- tax_Proro
worms_tax_tab_COI_full[worms_tax_tab_COI_full$scientificName == "Synchaeta kitina", ] <- tax_kitina

# Add identificationQualifier column from the tax_tab_COI_fixed df:
worms_tax_tab_COI_full$identificationQualifier <- tax_tab_COI_fixed[rownames(worms_tax_tab_COI_full), "identificationQualifier"]

# Still a few AphiaIDs missing - Verified with Matt Lemay that these taxa can be omitted from submission to OBIS:
worms_tax_tab_COI_full <- worms_tax_tab_COI_full %>%
  filter(!scientificName %in% c("Eukaryota", "Picobiliphyte", "Lagenidium caudatum", "Austrarchaea", "Trieres chinensis")) %>%
  select(-aphiaID)

# Change rownames to column, remove NAs"
worms_tax_tab_COI_full <- rownames_to_column(worms_tax_tab_COI_full, var = "ASV")
worms_tax_tab_COI_full <- worms_tax_tab_COI_full %>%
  mutate(across(everything(), ~ replace(.x, is.na(.x), "")))

write.csv(worms_tax_tab_COI_full, here("data", "metabarcoding", "COI", "FC_taxonomy_COI_worms.csv"))
