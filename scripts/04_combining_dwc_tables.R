
# Load packages

library(tidyr)
library(dplyr)
library(janitor)
library(here)
library(readxl)
library(readr)
library(lubridate)
library(hms)

#####################################################################################################################################

# Combine the event table created for the metabarcoding and that has been 
# parsed out of the BoLD DNA barcoding DwC table:
metabarcoding_event <- read_csv(here("obis", "interim_obis", "edna", "FC2022_event.csv"))
dna_barcoding_event <- read_csv(here("obis", "interim_obis", "bold", "FC2022_BoLD_event.csv")) %>%
  mutate(eventDate = as.character(eventDate))

FC2022_event <- bind_rows(metabarcoding_event, dna_barcoding_event)
FC2022_event <- obistools::flatten_event(FC2022_event)

# Remove certain columns that are not part of the event core:
FC2022_event <- FC2022_event %>% select(-c(station_depth, filter_start_time, filter_end_time, filter_vol_ml))

# Change NAs to blank cells:
FC2022_event[] <- lapply(FC2022_event, as.character)
FC2022_event[is.na(FC2022_event)] <- ""

# Save in obis folder
write_csv(FC2022_event, here("obis", "FC2022_event.csv"))

#####################################################################################################################################

# combine the occurrence table created for the metabarcoding and the table
# that has been parsed out of the BoLD DNA barcoding DwC table:
metabarcoding_occ <- read_csv(here("obis", "interim_obis", "edna", "FC2022_occ.csv"))
dna_barcoding_occ <- read_csv(here("obis", "interim_obis", "bold", "FC2022_BoLD_occ.csv"))
FC2022_occ <- bind_rows(metabarcoding_occ, dna_barcoding_occ)
FC2022_occ <- FC2022_occ %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

# Save in obis folder:
write_csv(FC2022_occ, here("obis", "FC2022_occ.csv"))

# Create a species list:
species_list <- FC2022_occ %>%
  select(verbatimIdentification, scientificName, scientificNameID) %>%
  distinct()

# Save in the data folder:
write_csv(species_list, here("data", "species_list.csv"))

#TODO: Manually inspect the occurrence table, look for inconsistencies.

#####################################################################################################################################

# Create the combined DNA derived data extension

metabarcoding_dna_ext <- read_csv(here("obis", "interim_obis", "edna", "FC2022_DNA.csv"))
dna_barcoding_dna_ext <- read_csv(here("obis", "interim_obis", "bold", "FC2022_BoLD_dna_ext.csv"))
FC2022_dna_ext_full <- bind_rows(metabarcoding_dna_ext, dna_barcoding_dna_ext)
FC2022_dna_ext_full <- FC2022_dna_ext_full %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

# Save in obis folder:
write_csv(FC2022_dna_ext_full, here("obis", "FC2022_dna_ext.csv"))

#####################################################################################################################################

# Combine the eMOF extensions

metabarcoding_emof <- read_csv(here("obis", "interim_obis", "edna", "FC2022_emof.csv")) %>%
  mutate(measurementValue = as.character(measurementValue))
dna_barcoding_emof <- read_csv(here("obis", "interim_obis", "bold", "FC2022_BoLD_emof.csv"))
FC2022_emof <- bind_rows(metabarcoding_emof, dna_barcoding_emof)
FC2022_emof <- FC2022_emof %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

# Save in obis folder:
write_csv(FC2022_emof, here("obis", "FC2022_emof.csv"))

#####################################################################################################################################

# The records from BoLD already have an associated materialSampleID included which
# is a unique identifier for the sample on that platform. The Hakai BoLD samples are 
# harvested annually into the iBOL: https://doi.org/10.15468/inygc6 
# We want to make sure users understand that the occurrenceID included in this dataset
# is the same as the occurrenceIDs found in the iBOL dataset. For this we use the resourceRelationship extension:

resourceRelationship <- dna_barcoding_occ %>%
  select(resourceID = occurrenceID, 
         relatedResourceID = materialSampleID) %>%
  mutate(relationshipOfResource = "same as",
         resourceRElationshipID = paste0(resourceID, "-rr-", stringr::str_pad(seq_along(resourceID), width = 3, pad = "0")))

# Save in obis folder:
write_csv(resourceRelationship, here("obis", "FC2022_rr_ext.csv"))

