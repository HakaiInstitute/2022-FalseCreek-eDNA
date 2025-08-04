
# Load packages

library(tidyr)
library(dplyr)
library(janitor)
library(here)
library(readxl)
library(readr)
library(lubridate)
library(hms)

# Combine the event table created for the metabarcoding and that has been 
# parsed out of the BoLD DNA barcoding DwC table:
metabarcoding_event <- read_csv(here("obis", "interim_obis", "edna", "FC2022_event.csv"))
dna_barcoding_event <- read_csv(here("obis", "interim_obis", "bold", "FC2022_BoLD_event.csv"))

FC2022_event <- bind_rows(metabarcoding_event, dna_barcoding_event)
FC2022_event <- obistools::flatten_event(FC2022_event)

# Save in obis folder
write_csv(FC2022_event, here("obis", "FC2022_event.csv"))

# combine the occurrence table created for the metabarcoding and the table
# that has been parsed out of the BoLD DNA barcoding DwC table:
metabarcoding_occ <- read_csv(here("obis", "interim_obis", "edna", "FC2022_occ.csv"))
dna_barcoding_occ <- read_csv(here("obis", "interim_obis", "bold", "FC2022_BoLD_occ.csv"))
FC2022_occ <- bind_rows(metabarcoding_occ, dna_barcoding_occ)
FC2022_occ <- FC2022_occ %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

# Save in obis folder:
write_csv(FC2022_occ, here("obis", "FC2022_occ.csv"))

#TODO: Manually inspect the occurrence table, look for inconsistencies.



