# Initialize the renv package
renv::init()

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

# Read in the sample information data
eDNA_meta <- read_csv(here("data", "FalseCreek_eDNA_Sample_Metadata_OBIS.csv")) %>%
  janitor::clean_names()

eDNA_meta <- eDNA_meta %>% 
  mutate(eventID = paste("Hakai", sample_id, sep = ":"),
         month = case_when(grepl("Sept", eDNA_meta$collection_date) ~ 9),
         day = as.numeric(stringr::str_extract(eDNA_meta$collection_date, "[0-9]+")),
         year = as.numeric(regmatches(eDNA_meta$collection_date, gregexpr("\\d{4}", eDNA_meta$collection_date))),
         eventDate = as.Date(paste(year, month, day, sep = "-")),
         eventType = "sampling",
         language = "en",
         license = "http://creativecommons.org/licenses/by/4.0/legalcode",
         bibliographicCitation = " ",
         accessRights = " ",
         rightsHolder = "Hakai Institute",
         institutionCode = "Hakai Institute",
         institutionID = "https://edmo.seadatanet.org/report/5148",
         modified = lubridate::today(),
         country = "Canada",
         countryCode = "CA",
         geodeticDatum = "WGS84",
         minimumDepthInMeters = depth_m,
         maximumDepthInMeters = depth_m) %>%
  dplyr::rename(
    verbatimLocality = location,
    decimalLatitude = lat,
    decimalLongitude = lon,
    verbatimEventDate = collection_date) %>%
  distinct()

# Save this data table in the obis folder
write_csv(eDNA_meta, here("obis", "eDNA_metadata.csv"))

#  Read in the ASV table:

FC_ASVs_12s <- read_csv(here("data", "FC_ASVs_12s.csv")) %>%
  dplyr::rename(ASV = row_name) %>%
  pivot_longer(FC_002:FC_074,
               names_to = "sample_id",
               values_to = "organismQuantity") %>%
  filter(organismQuantity != 0) %>%
  mutate(organismQuantityType = "DNA sequence reads",
         sampleSizeUnit = "DNA sequence reads",
         basisOfRecord = "MaterialSample") %>%
  group_by(sample_id) %>%
  mutate(sampleSizeValue = sum(organismQuantity),
         occurrenceID = paste(ASV, sample_id, "occ", row_number(), sep = "_")) %>%
  ungroup()

# Sample IDs in the metadata are hyphenated, in the ASV table underscored. Change this to ensure we can nest data:
FC_ASVs_12s$sample_id <- gsub("_", "-", FC_ASVs_12s$sample_id)

# Check for duplicate occurrenceIDs: 
print(sum(duplicated(FC_ASVs_12s$occurrenceID))) # should be 0. 

# Read in the sequence .fasta file: 
fasta_12S <- seqinr::read.fasta(file = here("data", "FC_12S_ASV_sequences.fasta"),
                                as.string = TRUE,
                                forceDNAtolower = FALSE)

ASV <- names(fasta_12S)
DNA_sequence <- unname(unlist(fasta_12S))
fasta_12S <- cbind(ASV, DNA_sequence) %>% as.data.frame()

# Read in the taxonomy table:
taxonomy_12S <- read_csv(here("data", "FC_taxonomy_12s.csv")) 

# Join the data tables so we know what we're working with:
DNAtable <- dplyr::left_join(taxonomy_12S, FC_ASVs_12s, by = "ASV")
DNAtable <- dplyr::left_join(DNAtable, eDNA_meta, by = "sample_id")
DNAtable <- dplyr::left_join(DNAtable, fasta_12S, by = "ASV")

# Occurrence table: each unique sequence by sample combination is considered one occurrence.

# The occurrence core table should contain the following:
# Both scientificName and verbatimIdentification are pulled from the species column - should
# we have to make any changes to the species' name to ensure matching with the WoRMS database
# we can make those changes to the scientificName and leave the verbatimIdentification intact.

# Fields below are taken from https://manual.obis.org/dna_data.html#id_16s-rrna-gene-metabarcoding-data-of-pico--to-mesoplankton 
DNA_occ <- DNAtable %>%
  select(occurrenceID,
         eventDate,
         basisOfRecord,
         scientificName = species,
         decimalLatitude,
         decimalLongitude,
         organismQuantity,
         organismQuantityType,
         sampleSizeValue,
         sampleSizeUnit,
         verbatimIdentification = species) %>%
  mutate(associatedSequences = NA, #TODO: should contain a reference to the URL domain where genetic sequence information associated with the occurrence can be found, e.g., for https://www.ncbi.nlm.nih.gov/bioproject/PRJNA887898/. 
         samplingProtocol = NA, #TODO: link to the relevant methods and processing steps as outlined on protocols.io? 
         identificationRemarks = NA, #TODO: Should be used to record information on how the taxonomic information of the occurrence was reached against which reference database, and with which confidence.
         identificationReferences = NA, #TODO: Should include a link to the bioinformatic pipeline or publication where the identification process is explained. 
         taxonConceptID = NA, #TODO: Can we pull this information from NCBI? For example, Clupea pallasii would be 'NCBI:txid30724'. 
         materialSampleID = NA)

# Add taxonomic information from the WoRMS database to the occurrence table: 
unique_spp <- unique(DNA_occ$scientificName)

# To match with the WoRMS taxonomic database we'll need to remove the ' sp.' from the scientificName
# Additionally, we need to change 'unknown', recommended practice: incertae sedis.
DNA_occ$scientificName <- gsub(" sp\\.", "", DNA_occ$scientificName)
DNA_occ$scientificName <- gsub("unknown", "incertae sedis", DNA_occ$scientificName)
DNA_occ$scientificName <- trimws(DNA_occ$scientificName, which = "right")

DNA_worms <- worrms::wm_records_names(unique(DNA_occ$scientificName), marine_only = F) %>% 
  dplyr::bind_rows() %>%
  filter(status == "accepted") %>%
  dplyr::rename(scientificName = scientificname) %>%
  select(scientificName, authority, rank, kingdom, phylum, class, order, family, genus,
         scientificNameID = lsid)

occ <- left_join(DNA_occ, DNA_worms, by = "scientificName")

# Save this data table in the obis folder
write_csv(occ, here("obis", "eDNA_occ.csv"))

# Create a DNA Derived Data extension: meant to capture information related to the sampled DNA, including sampling, processing, and other bioinformatic methods. The following (free-text) terms are required or highly
# recommended for eDNA and metabarcoding datasets.
# see details: https://rs.gbif.org/extension/gbif/1.0/dna_derived_data_2024-07-11.xml
DNA_extension <- DNAtable %>%
  select(occurrenceID, 
         DNA_sequence) %>%
  mutate(sop = NA, # Standard operating procedures used in assembly and/or annotation of genomes, metagenomes or environmental sequences. 
         samp_collec_device = "Niskin bottle",
         samp_taxon_id = NA, # This can be the NCBI Taxon ID of the sample. 
         target_gene = NA, # Targeted gene or locus name for marker gene studies, e.g. 12S or COI? 
         target_subfragment = NA, # Name of subfragement of a gene or locus.
         sample_vol_we_dna_ext = NA, # Volume (ml) or mass (g) of total collected sample processed for DNA extraction.
         pcr_primer_forward = NA, # Forward PCR Primer that was used to amplify the sequence of the targeted gene.
         pcr_primer_reverse = NA, # Reverse PCR primer that was used to amplify the sequence of the targeted gene.
         pcr_primer_name_forward = NA, # Name of the forward PCR primer
         pcr_primer_name_reverse = NA, # Name of the reverse PCR primer
         pcr_cond = NA, # Description of reaction conditions and components of PCR in the form of 'initial denaturation:94degC_1.5min; annealing:50_1..."
         annealingTemp = NA, # The reaction temperature during the annealing phase of PCR
         annealinTempUnit = NA, # Measurement unit of the reaction temperature during the annealing phase of PCR
         ampliconSize = NA, # The length of the amplicon in basepairs. 
         env_broad_scale = NA, # recommended to use ENVO's biomass classes to describe the major environmental system from which the sample was extracted, e.g. marine biome [ENVO:000000447]
         env_local_scale = NA, # recommended to use ENVO's biomass classes to describe the specific environmental system from which the sample was extracted, e.g. coastal waters [ENVO:00001250]
         env_medium = NA, # Probably this would almost always be: ocean water [ENVO:0002151]. 
         lib_layout = NA, # Specify whether to expect single, paired, or other configuration of reads.
         seq_meth = NA, # Sequencing method used, e.g. Sanger, ABI-solid. 
         otu_class_appr = NA, # Cutoffs and approach used when clustering new UViGs in "species-level" OTUs. Example: 95% ANI;85% AF; greedy incremental clustering
         otu_seq_comp_appr = NA, # Tool and thresholds used to compare sequences when computing "species-level" OTUs
         otu_db = NA # Reference database (i.e. sequences not generated as part of the current study) used to cluster new genomes in "species-level" OTUs, if any. E.g. NCBI Viral RefSeq;83
  )

# Save this data table in the obis folder
write_csv(DNA_extension, here("obis", "eDNA_extension.csv"))
