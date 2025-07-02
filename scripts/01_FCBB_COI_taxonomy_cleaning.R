# The following is a script by Libby Natola to create the 
# WoRMS taxa file for CO1 marker:

# load libraries ----
library(tidyverse)
library(taxize)
library(microViz)
library(data.table)
library(phyloseq)
library(microViz)

## COI ----
# COI is a marker that amplifies a broader taxonomic group, many of which won't be well annotated within WoRMS

### load/clean data ----
taxtab.COI <- fread("FC_taxonomy_COI.csv", sep=",", header=T, data.table=FALSE)
taxtab.COI <- taxtab.COI %>% 
  column_to_rownames("ASV") %>% 
  select(-V9)

# Convert to taxonomyTable class
tax_tab_COI <- tax_table(as.matrix(taxtab.COI))

# make uncultureds unknowns
tax_tab_COI[,"Species"] <- gsub("^uncultured.*", "unknown", tax_tab_COI[,"Species"])

# clean up unknown and propogate lowest rank out to the species column
tax_tab_COI_fixed <- tax_tab_COI %>% 
  tax_fix(
    min_length = 4,
    unknowns = c("unknown", "no identification", "unknown class", "unknown order", "unknown family", "unknown genus", "uncultured *"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified")

# get list of taxa, remove the taxon ranks (eg "Pleuronectidae family") from the species names
ncbi_taxa_COI <- as.data.frame(tax_tab_COI_fixed) %>%
  transmute(species_clean = str_remove_all(Species, " sp\\.$| Genus| Family| Order| Class| Phylum| Kingdom")) %>%
  pull(species_clean)

# keep the ASV names 
names(ncbi_taxa_COI) <- rownames(tax_tab_COI_fixed)

### convert to WoRMS ----
# get the unique taxa so we don't have to look up the same species in taxize multiple times
tax_lookup_tbl_COI <- tibble(
  ASV = names(ncbi_taxa_COI),
  Species = ncbi_taxa_COI
)

unique_species_COI <- unique(tax_lookup_tbl_COI$Species)

# map over unique species names
# query worms function - babysit this and manually input prompt responses - I generally choose the first one with status == "accepted"
worms_taxonomy_tbl_COI <- map_dfr(unique_species_COI, get_worms_lineage)

# join lineage info back to ASVs
worms_df_COI <- tax_lookup_tbl_COI %>%
  left_join(worms_taxonomy_tbl_COI, by = "Species") %>% 
  column_to_rownames("ASV")

# reorder the columns so the taxonomic ranks are in biological order
cols_order_COI <- c(
  "Kingdom", "Subkingdom", "Infrakingdom", "Phylum", "Phylum (Division)",
  "Subphylum", "Subphylum (Subdivision)", "Infraphylum", "Superclass", "Class", 
  "Subclass", "Infraclass", "Subterclass", "Order", "Suborder", "Infraorder", 
  "Superorder", "Family", "Subfamily", "Superfamily", "Tribe", "Genus", "Species")

# replace NAs with unknown so we can use taxfix to propagate out missing ranks
worms_df_COI <- worms_df_COI %>% 
  select(Species, all_of(setdiff(cols_order_COI, "Species"))) %>% # leave species first to prevent taxfix from erroring over having it alone when worms wasn't able to match taxon
  mutate(across(everything(), ~ replace_na(., "unknown")))

# propagate the taxa ranks out
worms_tax_tab_COI <- tax_table(as.matrix(worms_df_COI)) %>% 
  tax_fix(
    unknowns = c("unknown"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified")

# for anything worms couldn't identify the taxfix will fill them all out as the ASV. Get that ASV, replace all the columns with "no WoRMS ID" except the species column
worms_tax_tab_COI <- as.data.frame(worms_tax_tab_COI) %>%
  select(all_of(cols_order_COI)) %>% # reorder columns in correct taxonomic rank order
  mutate(no_worms = str_detect(Genus, "Species$")) %>% # column IDing rows that worms didn't ID
  mutate(across(
    -c(Species, no_worms),
    ~ if_else(no_worms & str_detect(., "Species$"), "No WoRMS ID", .))) %>% # change the ASVs to no WoRMS ID
  rownames_to_column("ASV") %>% # add the ASV back in so we can join
  left_join(tax_lookup_tbl_COI, by = "ASV", suffix = c("", "_lookup")) %>% # add tax lookup data
  select(ASV, cols_order_COI) %>% # remove those two columns we made
  column_to_rownames("ASV") 

write.csv(worms_tax_tab_COI, "FC_taxonomy_COI_worms.csv", quote = FALSE, row.names = TRUE)