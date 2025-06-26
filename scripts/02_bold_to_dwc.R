
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

# Read in the DNA barcoding DwC data table from BoLD: https://portal.boldsystems.org/recordset/DS-FCBB
DS_FCBB <- read_csv(here("data", "DS-FCBB.csv"))

# Create parent table:
project <- DS_FCBB %>% 
  mutate(month = as.numeric(format(mdy(eventDate), "%m")),
         day = as.numeric(format(mdy(eventDate), "%d")),
         year = as.numeric(format(mdy(eventDate), "%Y")),
         decimalLatitude = NULL,
         decimalLongitude = NULL,
         eventDate = paste(year, month, day, sep = "-"),
         eventType = samplingProtocol) %>%
  rename(bold_original_parentid = parentEventID) %>%
  select(eventID, bold_original_parentid,
         eventType,
         eventDate, month, day, year,
         samplingProtocol,
         habitat, 
         maximumDepthInMeters = verbatimDepth,
         minimumDepthInMeters = verbatimDepth,
         collectionCode, stateProvince, countryCode, fieldNumber,
         county, municipality, verbatimLocality, institutionCode, fieldNotes) %>%
  mutate(parentEventID = "Hakai-FCBB",
         eventDate = ifelse(bold_original_parentid == "FCRK21", "2022-09-03/05", 
                            ifelse(bold_original_parentid == "FCRK20", "2022-09-03/04", eventDate)),
         day = ifelse(bold_original_parentid %in% c("FCRK21", "FCRK20"), NA, day)) %>%
  mutate(eventRemarks = paste("BoLD original eventID:", eventID),
         institutionID = ifelse(institutionCode == "Hakai Institute", "https://edmo.seadatanet.org/report/5148", NA)) %>%
  group_by(eventID) %>%
  mutate(eventID = paste0(eventID, "-", stringr::str_pad(seq_along(eventID), width = 3, pad = "0"))) %>%
  ungroup() %>%
  select(-bold_original_parentid)

write_csv(project, here("obis", "interim_obis", "bold", "FC2022_BoLD_event.csv"))

## Occurrence extension:

occurrence <- DS_FCBB %>%
  select(eventID, materialSampleID, recordNumber, catalogNumber,
         occurrenceID, kingdom, phylum, class, order, family, subfamily,
         tribe, genus, scientificName, 
         subspecies = subpecies,
         scientificNameAuthorship, 
         identificationRemarks, identifiedBy, sex, lifeStage, occurrenceRemarks,
         organismRemarks, typeStatus, associatedOccurrences, associatedTaxa, recordedBy) %>%
  group_by(eventID) %>%
  mutate(eventID = paste0(eventID, "-", stringr::str_pad(seq_along(eventID), width = 3, pad = "0")),
         occurrenceID = paste0(eventID, "-", occurrenceID),
         basisOfRecord = "MaterialSample") %>%
  ungroup()

write_csv(occurrence, here("obis", "interim_obis", "bold", "FC2022_BoLD_occ.csv"))

## Measurement or Fact extension

emof <- DS_FCBB %>%
  select(eventID, occurrenceID, measurementDeterminedBy,
         measurementDeterminedDate, sex, lifeStage,
         measurementMethod) %>%
  group_by(eventID) %>%
  mutate(eventID = paste0(eventID, "-", stringr::str_pad(seq_along(eventID), width = 3, pad = "0")),
         occurrenceID = paste0(eventID, "-", occurrenceID)) %>%
  ungroup() %>%
  pivot_longer(sex:lifeStage,
               names_to = "measurementType",
               values_to = "measurementValue") %>%
  mutate(measurementDeterminedDate = as.Date(measurementDeterminedDate, format = "%m/%d/%Y"))

cleaned_emof <- emof %>%
  filter(!is.na(measurementValue)) %>%
  mutate(measurementID = paste0(occurrenceID, "-", measurementType)) %>%
  mutate(measurementValueID = case_when(
    measurementValue == "juvenile" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1127/",
    measurementValue == "larva" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1128/",
    measurementValue == "adult" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1116/",
    measurementValue == "zoea" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1137/",
    measurementValue == "nauplius" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1130/",
    measurementValue == "metatrocophore" ~ " "),
         measurementTypeID = case_when(
    measurementType == "lifeStage" ~ "https://vocab.nerc.ac.uk/collection/P01/current/LSTAGE01/"))

write_csv(cleaned_emof, here("obis", "interim_obis", "bold", "FC2022_BoLD_emof.csv"))
  
## DNA Derived 

dna_fcbb <- DS_FCBB %>%
  select(eventID, occurrenceID,
         DNA_sequence = dynamicProperties,
         materialSample = MaterialSample,
         materialSampleID)

write_csv(dna_fcbb, here("obis", "interim_obis", "bold", "FC2022_BoLD_dna_ext.csv"))
  


