## 12S

This folder contains the source files for the 12S metabarcoding data collected during the False Creek Bioblitz eDNA survey.

* `FC_taxonomy_12s.csv` — Detailed taxonomic information for each 12S amplicon sequence variant (ASV)
* `FC_taxonomy_12s_worms.csv` — The taxonomy file above, cleaned and matched to a WoRMS AphiaID (`scientificNameID`) by [scripts/01.1_FCBB_12S_taxonomy_cleaning.R](../../../scripts/01.1_FCBB_12S_taxonomy_cleaning.R).
* `FC_12S_ASV_sequences.fasta` — The 12S sequence for each amplicon sequence variant (ASV).
* `FC_ASVs_12s.csv` — The full matrix of 'samples x ASV' resulting from this survey for the 12S data.
* `FC12S_Metabarcoding_OBIS_sample_metadata.csv` — Additional sampling metadata recommended for publication to OBIS.
