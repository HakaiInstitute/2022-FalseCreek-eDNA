
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
sample_meta <- sample_meta %>% mutate(sample_id = gsub("_", "-", sample_id))

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
         eventDate = paste(year, month, day, sep = "-"),
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
FC2022_event <- FC2022_event %>%
  mutate(filter_start_time = as.character(filter_start_time),
         filter_end_time = as.character(filter_end_time))
FC2022_event[is.na(FC2022_event)] <- ""

# Save this data table in the obis folder
write_csv(FC2022_event, here("obis", "interim_obis", "edna", "FC2022_event.csv"))

#################################################################################################################################

## Occurrence extension

# Create separate data tables for the 12S, COI and 16S data, and join those.
# Start by reading in the 12S data:
FC_ASV_12S <- read_csv(here("data", "metabarcoding", "12S", "FC_ASVs_12s.csv"))
fasta_12S <- seqinr::read.fasta(file = here("data", "metabarcoding", "12S", "FC_12S_ASV_sequences.fasta"), as.string = TRUE, forceDNAtolower = FALSE)
sequence_meta_12S <- read_csv(here("data", "metabarcoding", "12S", "FC12S_Metabarcoding_OBIS_sample_metadata.csv"))
taxonomy_12S <- read_csv(here("data", "metabarcoding", "12S", "FC_taxonomy_12s_worms.csv"))
taxonomy_12S <- taxonomy_12S[, -1]

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
  mutate(eventID = paste("Hakai-FCBB", materialSampleID, sep = "-")) %>%
  group_by(materialSampleID, ASV) %>%
  mutate(occurrenceID = paste(eventID, ASV, "12S", "occ", sep = "-"),
         sampleSizeValue = sum(organismQuantity)) %>%
  ungroup()

# Match ASV to the sequence and taxa information:
ASV <- names(fasta_12S)
DNA_sequence <- unname(unlist(fasta_12S))
fasta_12S <- cbind(ASV, DNA_sequence) %>% as.data.frame()

# Join ASV to samples and taxonomic information:
DNAtable_12S <- dplyr::left_join(taxonomy_12S, FC_ASV_12S, by = "ASV")
DNAtable_12S <- dplyr::left_join(DNAtable_12S, event, by = c("materialSampleID", "eventID"))
DNAtable_12S <- dplyr::left_join(DNAtable_12S, fasta_12S, by = "ASV")

# Occurrence table: each unique sequence by sample combination is considered one occurrence.
sequence_meta_12S <- sequence_meta_12S %>% mutate(materialSampleID = gsub("_", "-", sequence_meta_12S$materialSampleID))

# Join DNA table with sequence metadata:
DNA_occ_12S <- left_join(DNAtable_12S, sequence_meta_12S, by = "materialSampleID") %>%
  select(ASV, eventID, occurrenceID, basisOfRecord, scientificName, scientificNameID,
         taxonRank, scientificNameAuthorship, kingdom, phylum, class, order, family, genus,
         organismQuantity, organismQuantityType, sampleSizeValue, sampleSizeUnit,
         verbatimIdentification, samplingProtocol, associatedSequences, identificationRemarks,
         materialSampleID)

# Join the occurrence data with the WoRMS taxa information:
FC_occ_12S <- DNA_occ_12S %>%
  mutate(occurrenceStatus = "present")

# Save this data table in the obis folder, omitting the first column (ASV)
write_csv(FC_occ_12S[, -1], here("obis", "interim_obis", "edna", "FC2022_occ_12S.csv"))

####################################################################################################################################

# Read in the data files for COI:
FC_ASV_COI <- read_csv(here("data", "metabarcoding", "COI", "FC_ASVs_COI.csv"))
FC_ASV_COI <- FC_ASV_COI %>% select(-FC_NegPCR, -FCXN001, -FCXN002, -FCXN003)
fasta_COI <- seqinr::read.fasta(file = here("data", "metabarcoding", "COI", "FC_CO1_ASV_sequences.fasta"),
                                as.string = TRUE,
                                forceDNAtolower = FALSE)

# Taxa data was claned by Libby Natolia following script 01_FCBB_COI_taxonomy_cleaning.R
sequence_meta_COI <- readxl::read_xlsx(here("data", "metabarcoding", "COI", "FC-COI_Metabarcoding_OBIS_sample_metadata.xlsx"), sheet = "Sheet1")
taxonomy_COI_worms <- read_csv(here("data", "metabarcoding", "COI", "FC_taxonomy_COI_worms.csv"))

# Occurrence table: each unique sequence by sample combination is considered one occurrence.
sequence_meta_COI <- sequence_meta_COI %>% mutate(materialSampleID = gsub("_", "-", sequence_meta_COI$materialSampleID))

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
  mutate(eventID = paste("Hakai-FCBB", materialSampleID, sep = "-")) %>%
  group_by(materialSampleID, ASV) %>%
  mutate(occurrenceID = paste(eventID, ASV, "COI", "occ", sep = "-"),
         sampleSizeValue = sum(organismQuantity)) %>%
  ungroup()

# Match ASV to the sequence and taxa information:
ASV_COI <- names(fasta_COI)
DNA_sequence_COI <- unname(unlist(fasta_COI))
fasta_COI <- cbind(ASV_COI, DNA_sequence_COI) %>% as.data.frame() %>%
  dplyr::rename(ASV = ASV_COI)

# Join ASV to samples and taxonomic information:
DNAtable_COI <- dplyr::left_join(taxonomy_COI_worms, FC_ASV_COI, by = "ASV")
DNAtable_COI <- dplyr::left_join(DNAtable_COI, event, by = c("materialSampleID", "eventID"))
DNAtable_COI <- dplyr::left_join(DNAtable_COI, fasta_COI, by = "ASV")

# Join DNA table with sequence metadata:
DNA_occ_COI <- left_join(DNAtable_COI, sequence_meta_COI, by = "materialSampleID") %>%
  select(ASV, eventID, occurrenceID, basisOfRecord, scientificName, scientificNameID,
         scientificNameAuthorship, kingdom, phylum, class, order, family, genus, 
         taxonRank, verbatimIdentification, identificationQualifier, organismQuantity, 
         organismQuantityType, sampleSizeValue, sampleSizeUnit, verbatimIdentification, 
         samplingProtocol, associatedSequences, identificationRemarks, materialSampleID)

# Extract ASV as well so that we can do manual referencing:
FC_occ_COI <- DNA_occ_COI %>%
  mutate(occurrenceStatus = "present")

# Save this data table in the obis folder
write_csv(FC_occ_COI[, -1], here("obis", "interim_obis", "edna", "FC2022_occ_COI.csv"))

#####################################################################################################################################

## 16S
# Read in the 16S data files:
# Read in the data files for COI:
FC_ASV_16S <- read_csv(here("data", "metabarcoding", "16S", "FC_ASV-table_16S.csv"))
fasta_16S <- seqinr::read.fasta(file = here("data", "metabarcoding", "16S", "FC-16S-ASV-sequences.fasta"),
                                as.string = TRUE,
                                forceDNAtolower = FALSE)

# Taxa data was claned by Libby Natolia following script 01_FCBB_COI_taxonomy_cleaning.R
sequence_meta_16S <- read_csv(here("data", "metabarcoding", "16S", "FC16S_Metabarcoding_OBIS_sample_metadata.csv"))
taxonomy_16S_worms <- read_csv(here("data", "metabarcoding", "16S", "FC_taxonomy_16S_worms.csv"))

# Occurrence table: each unique sequence by sample combination is considered one occurrence.
sequence_meta_16S <- sequence_meta_16S %>% mutate(materialSampleID = gsub("_", "-", sequence_meta_16S$materialSampleID))

#  Pivot the ASV table
FC_ASV_16S <- FC_ASV_16S %>%
  pivot_longer(FC_002:FC_074,
               names_to = "materialSampleID",
               values_to = "organismQuantity") %>%
  filter(organismQuantity != 0) %>%
  mutate(organismQuantityType = "DNA sequence reads",
         sampleSizeUnit = "DNA sequence reads",
         basisOfRecord = "MaterialSample") %>%
  mutate(materialSampleID = gsub("_", "-", materialSampleID)) %>%
  mutate(eventID = paste("Hakai-FCBB", materialSampleID, sep = "-")) %>%
  group_by(materialSampleID, ASV) %>%
  mutate(occurrenceID = paste(eventID, ASV, "16S", "occ", sep = "-"),
         sampleSizeValue = sum(organismQuantity)) %>%
  ungroup()

# Match ASV to the sequence and taxa information:
ASV_16S <- names(fasta_16S)
DNA_sequence_16S <- unname(unlist(fasta_16S))
fasta_16S <- cbind(ASV_16S, DNA_sequence_16S) %>% as.data.frame() %>%
  dplyr::rename(ASV = ASV_16S)

# Join ASV to samples and taxonomic information:
DNAtable_16S <- dplyr::left_join(taxonomy_16S_worms, FC_ASV_16S, by = "ASV")
DNAtable_16S <- dplyr::left_join(DNAtable_16S, event, by = c("materialSampleID", "eventID"))
DNAtable_16S <- dplyr::left_join(DNAtable_16S, fasta_16S, by = "ASV")

# Join DNA table with sequence metadata:
DNA_occ_16S <- left_join(DNAtable_16S, sequence_meta_16S, by = "materialSampleID") %>%
  select(ASV, eventID, occurrenceID, basisOfRecord, scientificName, scientificNameID,
         scientificNameAuthorship, kingdom, phylum, class, order, family, genus, 
         taxonRank, verbatimIdentification, organismQuantity, 
         organismQuantityType, sampleSizeValue, sampleSizeUnit, verbatimIdentification, 
         samplingProtocol, associatedSequences, identificationRemarks, materialSampleID)

# Extract ASV as well so that we can do manual referencing:
FC_occ_16S <- DNA_occ_16S %>%
  mutate(occurrenceStatus = "present")

# Save this data table in the obis folder
write_csv(FC_occ_16S[, -1], here("obis", "interim_obis", "edna", "FC2022_occ_16S.csv"))

#####################################################################################################################################

# Combine the occurrence data table for S12, COI and S16 data:
FC_occ <- bind_rows(FC_occ_S12, FC_occ_COI, FC_occ_16S)
FC_occ[is.na(FC_occ)] <- ""

# Check for duplicate occurrenceIDs:
print(sum(duplicated(FC_occ$occurrenceID))) # should be 0. 

# Save this data table in the obis folder
write_csv(FC_occ[, -1], here("obis", "interim_obis", "edna", "FC2022_occ.csv"))

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
