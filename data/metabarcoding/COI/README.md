## COI

This folder contains the source files for the COI metabarcoding data collected during the False Creek Bioblitz eDNA survey.

* `FC_taxonomy_COI_Revised_June2025.csv` — Detailed taxonomic information for each COI amplicon sequence variant (ASV).
* `FC_taxonomy_COI_worms.csv` — The taxonomy file above, cleaned and matched to a WoRMS AphiaID (`scientificNameID`) by [scripts/01.2_FCBB_COI_taxonomy_cleaning.R](../../../scripts/01.2_FCBB_COI_taxonomy_cleaning.R).
* `FC_CO1_ASV_sequences.fasta` — The COI sequence for each amplicon sequence variant (ASV).
* `FC_ASVs_COI.csv` — The full matrix of 'samples x ASV' resulting from this survey for the COI data.
* `FC-COI_Metabarcoding_OBIS_sample_metadata.xlsx` — Additional sampling metadata recommended for publication to OBIS.
