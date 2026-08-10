## 16S

This folder contains the source files for the 16S metabarcoding data collected during the False Creek Bioblitz eDNA survey.

* `FC_taxonomy_16S.csv` — Detailed taxonomic information for each 16S amplicon sequence variant (ASV).
* `FC_taxonomy_16S_cleaned.csv` — The taxonomy file above, cleaned by [scripts/01.3_FCBB_16S_taxonomy_cleaning.R](../../../scripts/01.3_FCBB_16S_taxonomy_cleaning.R) so that the `scientificName` matches the lowest known taxonomic rank for each ASV.
* `FC_taxonomy_16S_worms.csv` — The cleaned taxonomy file above, matched to a WoRMS AphiaID (`scientificNameID`) by [scripts/01.3_FCBB_16S_taxonomy_cleaning.R](../../../scripts/01.3_FCBB_16S_taxonomy_cleaning.R).
* `FC-16S-ASV-sequences.fasta` — The 16S sequence for each amplicon sequence variant (ASV).
* `FC_ASV-table_16S.csv` — The full matrix of 'samples x ASV' resulting from this survey for the 16S data.
* `FC16S_Metabarcoding_OBIS_sample_metadata.csv` — Additional sampling metadata recommended for publication to OBIS.
