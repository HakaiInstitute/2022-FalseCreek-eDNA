
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
lab_sheet <- read_excel(here("data", "DNA_barcoding", "DS-FCBB.xlsx"), sheet = 'Lab Sheet', skip = 2) %>% janitor::clean_names()
voucher_info <- read_excel(here("data", "DNA_barcoding", "DS-FCBB.xlsx"), sheet = 'Voucher Info', skip = 2) %>% janitor::clean_names()
taxonomy <- read_excel(here("data", "DNA_barcoding", "DS-FCBB.xlsx"), sheet = 'Taxonomy', skip = 2) %>% janitor::clean_names()
specimen_detail <- read_excel(here("data", "DNA_barcoding", "DS-FCBB.xlsx"), sheet = 'Specimen Details', skip = 2) %>% janitor::clean_names()
collection_data <- read_excel(here("data", "DNA_barcoding", "DS-FCBB.xlsx"), sheet = 'Collection Data', skip = 2) %>% janitor::clean_names()

# Create parent table:
project <- data.frame(
  eventDate = as.Date(collection_data$collection_date, format = "%d-%b-%Y"),
  eventID = collection_data$sample_id,
  fieldNotes = collection_data$collection_notes,
  maximumDepthInMeters = collection_data$depth,
  minimumDepthInMeters = collection_data$depth,
  habitat = collection_data$habitat,
  decimalLatitude = collection_data$lat,
  decimalLongitude = collection_data$lon,
  collectionCode = voucher_info$collection_code,
  institutionCode = voucher_info$institution_storing,
  catalogNumber = lab_sheet$catalog_num,
  fieldNumber = lab_sheet$field_id,
  recordNumber = lab_sheet$sample_id,
  materialSampleID = lab_sheet$sample_id,
  country = collection_data$country_ocean,
  stateProvince = collection_data$state_province,
  county = collection_data$region,
  municipality = collection_data$sector,
  verbatimLocality = collection_data$exact_site,
  samplingProtocol = collection_data$sampling_protocol,
  verbatimDepth = collection_data$depth,
  georeferenceSources = collection_data$gps_source,
  verbatimElevation = collection_data$elev,
  coordinatePrecision = collection_data$coordinate_accuracy)

project <- project %>% 
  mutate(
    month = as.numeric(format(eventDate, "%m")),
    day = as.numeric(format(eventDate, "%d")),
    year = as.numeric(format(eventDate, "%Y")),
    countryCode = "CA",
    parentEventID = "Hakai-FCBB",
    institutionID = ifelse(institutionCode == "Hakai Institute", "https://edmo.seadatanet.org/report/5148",
                                ifelse(grepl("Beaty Biodiversity Museum", institutionCode), "http://biocol.org/urn:lsid:biocol.org:col:15106", NA))) %>%
  group_by(eventID) %>%
  mutate(eventID = paste0(eventID, "-", stringr::str_pad(seq_along(eventID), width = 3, pad = "0"))) %>%
  mutate(eventDate = as.character(eventDate)) %>%
  ungroup()

# Remove all rows that are completely NA:
project <- project %>% select(where(~ !all(is.na(.))))

# Replace NA values with blank cells
project[is.na(project)] <- ""

write_csv(project, here("obis", "interim_obis", "bold", "FC2022_BoLD_event.csv"))

## Occurrence extension:
occurrence <- data.frame(
  eventID = project$eventID,
  materialSampleID = project$materialSampleID,
  recordNumber = lab_sheet$sample_id,
  catalogNumber = lab_sheet$catalog_num,
  phylum = taxonomy$phylum,
  class = taxonomy$class,
  order = taxonomy$order,
  family = taxonomy$family,
  subfamily = taxonomy$subfamily,
  tribe = taxonomy$tribe,
  genus = taxonomy$genus,
  scientificName = lab_sheet$identification,
  identificationRemarks = taxonomy$identification_method,
  identifiedBy = taxonomy$identifier,
  sex = specimen_detail$sex,
  lifeStage = specimen_detail$life_stage,
  occurrenceRemarks = taxonomy$taxonomy_notes,
  typeStatus = specimen_detail$voucher_status,
  materialSample = specimen_detail$tissue_descriptor,
  associatedOccurrences = specimen_detail$associated_specimens,
  associatedTaxa = specimen_detail$associated_taxa,
  recordedBy = collection_data$collectors,
  organismRemarks = specimen_detail$notes
)

# Some important taxonomic information like kingdom and scientificNameAuthorship is missing:
unique_spp <- unique(occurrence$scientificName)

get_tax_info <- function(name) {
  rec <- tryCatch(worrms::wm_records_name(name), error = function(e) NA)
  if (is.data.frame(rec) && nrow(rec) > 0) {
    return(data.frame(
      name = name,
      lsid = rec$lsid[1],
      authority = rec$authority[1],
      kingdom = rec$kingdom[1],
      rank = rec$rank[1],
      stringsAsFactors = F
    ))
  } else {
    return(data.frame(
      name = name,
      lsid = NA,
      authority = NA,
      kingdom = NA,
      rank = NA,
      stringsAsFactors = F
    ))
  }
}

# Fetch the data:
tax_info <- do.call(rbind, lapply(unique_spp, get_tax_info)) %>%
  dplyr::rename(scientificName = name,
                taxonRank = rank,
                scientificNameAuthorship = authority,
                scientificNameID = lsid)

# Join together:
occurrence <- occurrence %>% left_join(tax_info, by = "scientificName") %>%
  mutate(basisOfRecord = "MaterialSample",
         occurrenceStatus = "present",
         occurrenceID = paste0(eventID, "-", stringr::str_pad(seq_along(eventID), width = 3, pad = "0")))

# For Cumella vulgaris the function could not find an appropriate lsid, add this manually:
#TODO: confirm with Matt Lemay that the LSID is correct:
occurrence$scientificNameID = ifelse(occurrence$scientificName == "Cumella vulgaris", "urn:lsid:marinespecies.org:taxname:182534", occurrence$scientificNameID)
occurrence$scientificNameAuthorship = ifelse(occurrence$scientificName == "Cumella vulgaris", "Hart, 1930", occurrence$scientificNameAuthorship)
occurrence$taxonRank <- ifelse(occurrence$scientificName == "Cumella vulgaris", "Species", occurrence$taxonRank)
occurrence$kingdom <- ifelse(occurrence$scientificName == "Cumella vulgaris", "Animalia", occurrence$kingdom)

# Remove columns with only NA and replace NA with blank cells:
occurrence <- occurrence %>% select(where(~ !all(is.na(.))))
occurrence[is.na(occurrence)] <- ""

# Save locally:
write_csv(occurrence, here("obis", "interim_obis", "bold", "FC2022_BoLD_occ.csv"))

## Measurement or Fact extension

emof <- data.frame(
  eventID = project$eventID,
  occurrenceID = occurrence$occurrenceID,
  sex = specimen_detail$sex,
  lifeStage = specimen_detail$life_stage,
  measurementMethod = taxonomy$identification_method) %>%
  pivot_longer(sex:lifeStage,
               names_to = "measurementType",
               values_to = "measurementValue") %>%
  filter(!is.na(measurementValue)) %>%
  mutate(measurementID = paste0(occurrenceID, "-", measurementType)) %>%
  mutate(measurementValueID = case_when(
    measurementValue == "juvenile" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1127/",
    measurementValue == "larva" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1128/",
    measurementValue == "adult" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1116/",
    measurementValue == "zoea" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1137/",
    measurementValue == "medusa" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1129/",
    measurementValue == "nauplius" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1130/",
    measurementValue == "megalopa" ~ "https://vocab.nerc.ac.uk/collection/S11/current/S1167/",
    measurementValue == "F" ~ "https://vocab.nerc.ac.uk/collection/S10/current/S102/",
    measurementValue == "M" ~ "https://vocab.nerc.ac.uk/collection/S10/current/S103/",
    measurementValue == "metatrocophore" ~ " "),
         measurementTypeID = case_when(
    measurementType == "lifeStage" ~ "https://vocab.nerc.ac.uk/collection/P01/current/LSTAGE01/",
    measurementType == "sex" ~ "https://vocab.nerc.ac.uk/collection/P01/current/ENTSEX01/"))

emof[is.na(emof)] <- ""

write_csv(emof, here("obis", "interim_obis", "bold", "FC2022_BoLD_emof.csv"))
  
## DNA Derived 

# Read lines from file:
lines <- readLines(here("data", "DNA_barcoding", "FCBB.fasta.txt"))

# Extract headers and sequence data:
headers <- lines[grepl("^>", lines)]
sequences <- lines[!grepl("^>", lines)]

header_parts <- strsplit(sub("^>", "", headers), "\\|")
header_df <- do.call(rbind, header_parts)
colnames(header_df) <- c("materialSampleID", "scientificName", "fieldNumber")

# Combine into dataframe:
dna <- data.frame(
  process_id = header_df[, 1],
  scientificName = header_df[, 2],
  materialSampleID = header_df[, 3],
  DNA_sequence = sequences,
  stringsAsFactors = F)

# Join with occurrence, based on materialSampleID:
dna_extension <- occurrence %>%
  left_join(dna, by = "materialSampleID") %>%
  select(eventID, occurrenceID,
         DNA_sequence,
         materialSample, materialSampleID)

# Remove taxa where DNA_sequence is NA. These are already captured in the
# occurrence extension:
dna_extension <- dna_extension[!is.na(dna_extension$DNA_sequence), ]

# Save locally:
write_csv(dna_extension, here("obis", "interim_obis", "bold", "FC2022_BoLD_dna_ext.csv"))
  


