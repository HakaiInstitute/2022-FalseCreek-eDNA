
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


# Combine the event table created for the metabarcoding and that has been 
# parsed out of the BoLD DNA barcoding DwC table:
metabarcoding_event <- read_csv(here("obis", "interim_obis", "edna", "FC2022_event.csv"))
dna_barcoding_event <- read_csv(here("obis", "interim_obis", "bold", "FC2022_BoLD_event.csv"))

FC2022_event <- bind_rows(metabarcoding_event, dna_barcoding_event)
FC2022_event <- obistools::flatten_event(FC2022_event)

# Save in obis folder
write_csv(FC2022_event, here("obis", "FC2022_event.csv"))