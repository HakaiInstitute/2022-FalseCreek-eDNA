
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

# Download the DNA barcoding DwC data table from BoLD: https://portal.boldsystems.org/recordset/DS-FCBB
#TODO: Figure out if this is something that can happen automatically using the BoLD API?

DS_FCBB <- read_csv(here("data", "DS-FCBB.csv"))

# Create parent table:
project <- DS_FCBB %>% 
  mutate(month = as.numeric(format(mdy(eventDate), "%m")),
         day = as.numeric(format(mdy(eventDate), "%d")),
         year = as.numeric(format(mdy(eventDate), "%Y")),
         eventDate = paste(year, month, day, sep = "-"),
         eventType = samplingProtocol) %>%
  rename(bold_original_parentid = parentEventID) %>%
  select(eventID, bold_original_parentid,
         eventType,
         eventDate, month, day, year,
         decimalLatitude,
         decimalLongitude,
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
         institutionID = ifelse(institutionCode == "Hakai Institute", "https://edmo.seadatanet.org/report/5148", 
                                ifelse(grepl("Beaty Biodiversity Museum", DS_FCBB$institutionCode), "http://biocol.org/urn:lsid:biocol.org:col:15106", NA))) %>%
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
         basisOfRecord = "MaterialSample",
         occurrenceStatus = "present") %>%
  ungroup()

# There are entries where scientificName is NA, but there is taxonomic information provided. Where
# scientificName is NA, the lowest level taxa hierarchical information should be used:

tax_cols <- c("kingdom", "phylum", "class", "order", "family", "subfamily", "tribe", "genus", "scientificName")
target_col <- "scientificName"

fill_na_with_last_non_na <- function(occurrence, tax_cols, target_col) {
  
  # Ensure target_col is in cols
  if (!(target_col %in% tax_cols)) stop("target_col must be one of tax_cols")
  
  target_index <- match(target_col, tax_cols)
  
  # Loop through rows where scientificName is NA
  for (i in which(is.na(occurrence[[target_col]]))) {
    row_vals <- unlist(occurrence[i, tax_cols], use.names = FALSE)
    
    # Look for the last non-NA value in columns before scientificName
    if (target_index > 1) {
      prior_vals <- row_vals[1:(target_index - 1)]
      last_val <- tail(na.omit(prior_vals), 1)
      
      if (length(last_val) > 0) {
        occurrence[i, target_col] <- last_val
        }
    }
  }
  
  return(occurrence)
}

occ_filled <- fill_na_with_last_non_na(occurrence, tax_cols, target_col)

# Add scientificNameID, grabbing this information from WoRMS:
get_taxon_info <- function(scientific_name) {
  
  tryCatch({
    record <- worrms::wm_records_name(scientific_name)
    if (length(record) > 0) {
      return(list(
        lsid = record$lsid[1],
        taxonRank = record$rank[1]))
    } else {
      return(list(lsid = NA, taxonRank = NA))
    }
  }, error = function(e) {
    return(list(lsid = NA, taxonRank = NA))
  })
}

taxon_info <- lapply(occurrence$scientificName, get_taxon_info)

occ_filled$scientificNameID <- sapply(taxon_info, `[[`, "lsid")
occ_filled$taxonRank <- sapply(taxon_info, `[[`, "taxonRank")

# For Cumella vulgaris the function could not find an appropriate lsid, add this manually:
#TODO: confirm with Matt Lemay that the LSID is correct, consider changing scientificName in
# the original source files:
occ_filled$scientificNameID <- ifelse(occ_filled$scientificName == "Cumella vulgaris", 
                                    "urn:lsid:marinespecies.org:taxname:182534", 
                                    occ_filled$scientificNameID)
occ_filled$taxonRank <- ifelse(occ_filled$scientificName == "Cumella vulgaris",
                               "Species", occ_filled$taxonRank)
                                
# Select columns where at least one value is not NA (i.e., remove columns that are completely NA)
occurrence_cleaned <- occ_filled %>% select(where (~ any(!is.na(.))))

write_csv(occurrence_cleaned, here("obis", "interim_obis", "bold", "FC2022_BoLD_occ.csv"))

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
  


