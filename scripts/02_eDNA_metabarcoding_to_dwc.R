
# Install packages:
install.packages("devtools")
devtools::install_github("iobis/obistools")

# Load packages
library(tidyr)
library(dplyr)
library(obistools)
library(janitor)
library(here)
library(readxl)
library(readr)
library(lubridate)
library(hms)

# Read in the metabarcoding eDNA metadata tables:
sample_meta <- read_csv(here("data", "metabarcoding", "FalseCreek_eDNA_sample_metadata.csv")) %>% janitor::clean_names()

# Create project level information:
FCBB <- data.frame(
  eventID = "Hakai-FCBB",
  language = "en",
  license = "http://creativecommons.org/licenses/by/4.0/legalcode",
  bibliographicCitation = "Lemay, M., & Kellogg, C. (2025). Environmental DNA metabarcoding data from the False Creek Bioblitz, 2022 (v1.0). Hakai Institute. https://doi.org/10.21966/8sgn-jg38",
  rightsHolder = "Hakai Institute",
  modified = lubridate::today(),
  country = "Canada",
  countryCode = "CA",
  geodeticDatum = "WGS84",
  stateProvince = "British Columbia",
  county = "Burrard Inlet", 
  municipality = "False Creek")

# Create table for sampling event metadata: 
event <- sample_meta %>% 
  mutate(parentEventID = "Hakai-FCBB",
         eventType = "Niskin bottle sampling",
         eventID = paste(parentEventID, sample_id, sep = "-"),
         month = case_when(grepl("Sept", sample_meta$collection_date) ~ 9),
         day = as.numeric(stringr::str_extract(sample_meta$collection_date, "[0-9]+")),
         year = as.numeric(regmatches(sample_meta$collection_date, gregexpr("\\d{4}", sample_meta$collection_date))),
         eventDate = (paste(year, month, day, sep = "-")),
         minimumDepthInMeters = depth_m,
         maximumDepthInMeters = depth_m,
         institutionCode = "Hakai Institute",
         institutionID = "https://edmo.seadatanet.org/report/5148") %>%
  dplyr::rename(
    materialSampleID = sample_id,
    verbatimLocality = location,
    decimalLatitude = lat,
    decimalLongitude = lon,
    verbatimEventDate = collection_date) %>%
  distinct() %>%
  select(-c(station_id, treatment, depth_m))

# Combine for event table eDNA metabarcoding, and flatten:
FC2022_event <- bind_rows(FCBB, event)
FC2022_event <- obistools::flatten_event(FC2022_event)
FC2022_event[is.na(FC2022_event)] <- ""

# Save this data table in the obis folder
write_csv(FC2022_event, here("obis", "interim_obis", "edna", "FC2022_event.csv"))

#################################################################################################################################

## Occurrence extension

# Create separate data tables for the 12S, COI and 16S data, and join those.
# Start by reading in the 12S data:
FC_ASV_12S <- read_csv(here("data", "metabarcoding", "12S", "FC_ASVs_12s.csv"))
fasta_12S <- seqinr::read.fasta(file = here("data", "metabarcoding", "12S", "FC_12S_ASV_sequences.fasta"), as.string = TRUE, forceDNAtolower = FALSE)
taxonomy_12S <- read_csv(here("data", "metabarcoding", "12S", "FC_taxonomy_12s_worms.csv"))
sequence_meta_12S <- read_csv(here("data", "metabarcoding", "12S", "FC12S_Metabarcoding_OBIS_sample_metadata.csv"))

#  Pivot the ASV table
FC_ASV_12S <- FC_ASV_12S %>%
  dplyr::rename(ASV = row_name) %>%
  pivot_longer(FC_002:FC_074,
               names_to = "materialSampleID",
               values_to = "organismQuantity") %>%
  filter(organismQuantity != 0) %>%
  mutate(organismQuantityType = "DNA sequence reads",
         sampleSizeUnit = "DNA sequence reads",
         basisOfRecord = "MaterialSample") %>%
  mutate(materialSampleID = gsub("_", "-", materialSampleID)) %>%
  mutate(eventID = paste("Hakai", materialSampleID, sep = "-")) %>%
  group_by(materialSampleID, ASV) %>%
  mutate(occurrenceID = paste(eventID, ASV, "12S", "occ", sep = "-"),
         sampleSizeValue = sum(organismQuantity)) %>%
  ungroup()

# Match ASV to the sequence and taxa information:
ASV <- names(fasta_12S)
DNA_sequence <- unname(unlist(fasta_12S))
fasta_12S <- cbind(ASV, DNA_sequence) %>% as.data.frame()

# Upon manual inspection, for some ASVs information is not provided at species level, 
# for these taxa information is e.g. provided up to the genus level and species
# information is recorded as "unknown". Change 'unknown' to NA
tax_table_12S <- taxonomy_12S %>%
  mutate(across(everything(), ~na_if(., "unknown"))) 
colnames(tax_table_12S)[1] <- "ASV"

# Grab the most granular taxonomic data (i.e., the last column populated before
# 'NA' information is provided):
tax_table_12S$scientificName <- apply(tax_table_12S, 1, function(row) {
  last_value <- tail(na.omit(row), 1) 
  if (length(last_value) == 0) NA else last_value
})

# Join ASV to samples and taxonomic information:
DNAtable_12S <- dplyr::left_join(tax_table_12S, FC_ASV_12S, by = "ASV")
DNAtable_12S <- dplyr::left_join(DNAtable_12S, event, by = c("materialSampleID", "eventID"))
DNAtable_12S <- dplyr::left_join(DNAtable_12S, fasta_12S, by = "ASV")

# Capture original recording of species under verbatimIdentifation:
DNAtable_12S$verbatimIdentification <- DNAtable_12S$Species

# Occurrence table: each unique sequence by sample combination is considered one occurrence.
sequence_meta_12S <- sequence_meta_12S %>% mutate(materialSampleID = gsub("_", "-", sequence_meta_12S$materialSampleID))

# Join DNA table with sequence metadata:
DNA_occ_12S <- left_join(DNAtable_12S, sequence_meta_12S, by = "materialSampleID") %>%
  select(eventID,
         occurrenceID,
         basisOfRecord,
         scientificName,
         organismQuantity,
         organismQuantityType,
         sampleSizeValue,
         sampleSizeUnit,
         verbatimIdentification,
         samplingProtocol,
         associatedSequences,
         identificationRemarks,
         materialSampleID)

# To match with the WoRMS taxonomic database we'll need to remove the ' sp.' from the scientificName
DNA_occ_12S$scientificName <- trimws(gsub(" sp\\..*$", "", DNA_occ_12S$scientificName))

DNA_worms_12S <- worrms::wm_records_names(unique(DNA_occ_12S$scientificName), marine_only = F) %>% 
  dplyr::bind_rows() %>%
  filter(status == "accepted") %>%
  dplyr::rename(scientificName = scientificname) %>%
  select(scientificName, 
         scientificNameAuthorship = authority, 
         taxonRank = rank, 
         kingdom, phylum, class, order, family, genus,
         scientificNameID = lsid)

# Join the occurrence data with the WoRMS taxa information:
FC_occ_S12 <- left_join(DNA_occ_12S, DNA_worms_12S, by = "scientificName") %>%
  mutate(occurrenceStatus = "present")

# Save this data table in the obis folder
write_csv(FC_occ_S12, here("obis", "interim_obis", "edna", "FC2022_occ_S12.csv"))

####################################################################################################################################

# Read in the data files for COI:
FC_ASV_COI <- read_csv(here("data", "metabarcoding", "COI", "FC_ASVs_COI.csv"))
fasta_COI <- seqinr::read.fasta(file = here("data", "metabarcoding", "COI", "FC_CO1_ASV_sequences.fasta"),
                                as.string = TRUE,
                                forceDNAtolower = FALSE)

# Taxa data was claned by Libby Natolia following script 01_FCBB_COI_taxonomy_cleaning.R
taxonomy_COI_worms <- read_csv(here("data", "metabarcoding", "COI", "FC_taxonomy_COI_worms.csv"))
sequence_meta_COI <- readxl::read_xlsx(here("data", "metabarcoding", "COI", "FC-COI_Metabarcoding_OBIS_sample_metadata.xlsx"), sheet = "Sheet1")

# Occurrence table: each unique sequence by sample combination is considered one occurrence.
sequence_meta_COI <- sequence_meta_COI %>% mutate(materialSampleID = gsub("_", "-", sequence_meta_COI$materialSampleID))

colnames(taxonomy_COI_worms)[1] <- "ASV"
taxonomy_COI_worms <- taxonomy_COI_worms %>%
  select(ASV,
         kingdom = Kingdom,
         phylum = Phylum,
         class = Class,
         order = Order,
         family = Family,
         genus = Genus,
         species = Species) %>%
  mutate(scientificName = species)

#  Pivot the ASV table
FC_ASV_COI <- FC_ASV_COI %>%
  dplyr::rename(ASV = row_names) %>%
  pivot_longer(FC_002:FC_074,
               names_to = "materialSampleID",
               values_to = "organismQuantity") %>%
  filter(organismQuantity != 0) %>%
  mutate(organismQuantityType = "DNA sequence reads",
         sampleSizeUnit = "DNA sequence reads",
         basisOfRecord = "MaterialSample") %>%
  mutate(materialSampleID = gsub("_", "-", materialSampleID)) %>%
  mutate(eventID = paste("Hakai", materialSampleID, sep = "-")) %>%
  group_by(materialSampleID, ASV) %>%
  mutate(occurrenceID = paste(eventID, ASV, "COI", "occ", sep = "-"),
         sampleSizeValue = sum(organismQuantity)) %>%
  ungroup()

# Match ASV to the sequence and taxa information:
ASV_COI <- names(fasta_COI)
DNA_sequence_COI <- unname(unlist(fasta_COI))
fasta_COI <- cbind(ASV_COI, DNA_sequence_COI) %>% as.data.frame() %>%
  dplyr::rename(ASV = ASV_COI)

# As per conversation with Matt Lemay, for species names with sp. followed by a string
# we can remove this, but keep original recording under verbatimIdentification:
taxonomy_COI_worms$verbatimIdentification <- taxonomy_COI_worms$species
taxonomy_COI_worms$scientificName <- gsub("sp\\..*", "", taxonomy_COI_worms$scientificName)
taxonomy_COI_worms$scientificName <- gsub("CCMP1545|CMC01", "", taxonomy_COI_worms$scientificName) %>% trimws()

# Capture identificationQualifiers and remove them from the scientificName:
taxonomy_COI_worms <- taxonomy_COI_worms %>%
  mutate(identificationQualifier = case_when(
    grepl("cf. promare", taxonomy_COI_worms$scientificName) ~ "cf. promare",
    grepl("cf. depressum", taxonomy_COI_worms$scientificName) ~ "cf. depressum"))
taxonomy_COI_worms$scientificName <- gsub("cf. promare|cf. depressum", "", taxonomy_COI_worms$scientificName) %>% trimws()

# Join ASV to samples and taxonomic information:
DNAtable_COI <- dplyr::left_join(taxonomy_COI_worms, FC_ASV_COI, by = "ASV")
DNAtable_COI <- dplyr::left_join(DNAtable_COI, event, by = c("materialSampleID", "eventID"))
DNAtable_COI <- dplyr::left_join(DNAtable_COI, fasta_COI, by = "ASV")

# Join DNA table with sequence metadata:
DNA_occ_COI <- left_join(DNAtable_COI, sequence_meta_COI, by = "materialSampleID") %>%
  select(eventID,
         occurrenceID,
         basisOfRecord,
         scientificName,
         organismQuantity,
         organismQuantityType,
         sampleSizeValue,
         sampleSizeUnit,
         verbatimIdentification,
         samplingProtocol,
         associatedSequences,
         identificationRemarks,
         materialSampleID)

# For each row, grab the following taxonomic hierarchical information: scientificNameAuthorship, taxonRank and 
# scientificNameID:
DNA_worms_COI <- worrms::wm_records_names(unique(taxonomy_COI_worms$scientificName), marine_only = F) %>%
  dplyr::bind_rows() %>%
  filter(status == "accepted") %>%
  filter(AphiaID != 119270) %>%
  select(scientificName = scientificname,
         taxonRank = rank,
         scientificNameAuthorship = authority,
         kingdom, phylum, class, order, family, genus,
         scientificNameID = lsid)

FC_occ_COI <- left_join(DNA_occ_COI, DNA_worms_COI, by = "scientificName") %>%
  mutate(occurrenceStatus = "present") %>%
  mutate(ASV = stringr::str_extract(occurrenceID, "ASV\\d+"))

# Save this data table in the obis folder
write_csv(FC_occ_COI, here("obis", "interim_obis", "edna", "FC2022_occ_COI.csv"))

#####################################################################################################################################

## 16S
# Read in the 16S data files:




#####################################################################################################################################

# Combine the occurrence data table for S12, COI and S16 data:
FC_occ <- bind_rows(FC_occ_S12, FC_occ_COI)

# Check for duplicate occurrenceIDs:
print(sum(duplicated(FC_occ$occurrenceID))) # should be 0. 

# Save this data table in the obis folder
write_csv(FC_occ, here("obis", "interim_obis", "edna", "FC2022_occ.csv"))

# Create a DNA Derived Data extension: 
# see details: https://rs.gbif.org/extension/gbif/1.0/dna_derived_data_2024-07-11.xml
DNA_extension_12S <- left_join(DNAtable_12S, sequence_meta_12S, by = "materialSampleID")
DNA_extension_16S <- left_join(DNAtable_16S, sequence_meta_16S, by = "materialSampleID")
DNA_extension_COI <- left_join(DNAtable_COI, sequence_meta_COI, by = "materialSampleID")

FC2022_DNA_extension <- bind_rows(DNA_extension_12S,
                                  DNA_extension_16S,
                                  DNA_extension_COI) %>%
  select(eventID, occurrenceID, DNA_sequence, sop,
         samp_collect_device = samp_collec_method, 
         target_gene,
         target_subfragment,
         sample_vol_we_dna_ext,
         pcr_primer_forward,
         pcr_primer_reverse,
         pcr_primer_name_forward,
         pcr_primer_name_reverse,
         env_broad_scale,
         env_local_scale,
         env_medium,
         lib_layout,
         seq_meth,
         otu_class_appr,
         otu_seq_comp_appr,
         otu_db)

# Save this data table in the obis folder
write_csv(DNA_extension, here("obis", "interim_obis", "edna", "FC2022_eMOF.csv"))
