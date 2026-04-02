# author: Libby Natola
# date: 30 May 2025
# script to test converting NCBI taxonomy to WoRMS taxonomy for OBIS upload. Test files are from False Creek bioblitz
# Minor updates made by Tim van der Stap, July 4, 2025, primarly to add some Darwin Core terms.

# Install requirements for microViz and phyloseq:
install.packages(
  "microViz",
  repos = c(
    davidbarnett = "https://david-barnett.r-universe.dev",
    getOption("repos")
  )
)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

# Functions used ----
# To grab authority
get_authority <- function(id) {
  tryCatch(
    {
      rec <- worrms::wm_record(id = id)
      rec$authority
    },
    error = function(e) NA
  )
}

# load libraries ----
library(tidyverse)
library(taxize) # c 0.9.100
library(microViz)
library(data.table)
library(phyloseq)
library(here)

# COI is a marker that amplifies a broader taxonomic group, many of which won't be well annotated within WoRMS
taxtab.COI <- fread(
  here("data", "metabarcoding", "COI", "FC_taxonomy_COI_Revised_June2025.csv"),
  sep = ",",
  header = T,
  data.table = FALSE
)
colnames(taxtab.COI)[-1] <- str_to_title(colnames(taxtab.COI)[-1])
taxtab.COI <- taxtab.COI %>% column_to_rownames("ASV")

# convert to taxonomyTable class
tax_tab_COI <- as.data.frame(tax_table(as.matrix(taxtab.COI)))

# make uncultureds unknowns
tax_tab_COI$Species <- gsub("^uncultured.*", "unknown", tax_tab_COI$Species)

# Convert to taxonomyTable class
tax_tab_COI <- tax_tab_COI |>
  mutate(
    across(everything(), ~ na_if(., "unknown")),
    across(everything(), ~ na_if(., "no identification")),
    across(everything(), ~ na_if(., "unknown class")),
    across(everything(), ~ na_if(., "unknown order")),
    across(everything(), ~ na_if(., "unknown family")),
    across(everything(), ~ na_if(., "unknown genus")),
    across(everything(), ~ na_if(., "uncultured *"))
  )

# Grab the most granular taxonomic data (i.e., the last column populated before
# 'NA' information is provided):
tax_tab_COI$scientificName <- apply(tax_tab_COI, 1, function(row) {
  last_value <- tail(na.omit(row), 1)
  if (length(last_value) == 0) NA else last_value
})

# Add in verbatimIdentification to capture original recorded name:
verbatimIdentification <- tax_tab_COI$scientificName
names(verbatimIdentification) <- rownames(tax_tab_COI)
verbatimIdentification <- as.data.frame(verbatimIdentification)

# Capture identificationQualifiers and remove them from the scientificName.
# Include them in a separate column:
identificationQualifier <- tax_tab_COI |>
  mutate(
    identificationQualifier = case_when(
      grepl("cf. promare", Species) ~ "cf. promare",
      grepl("cf. depressum", Species) ~ "cf. depressum",
      grepl("cf. urtica", Species) ~ "cf. urtica",
      grepl("aff. granulata", Species) ~ "cf. aff. granulata"
    )
  ) |>
  select(identificationQualifier)

tax_tab_COI_fixed <- tax_tab_COI |>
  mutate(
    scientificName = gsub(
      "cf. promare|cf. depressum|cf. urtica|aff. granulata",
      "",
      tax_tab_COI$scientificName
    )
  )

# remove these strings
tax_tab_COI_fixed$scientificName <- gsub(
  "sp\\..*",
  "",
  tax_tab_COI_fixed$scientificName
)
tax_tab_COI_fixed$scientificName <- gsub(
  "CCMP1545|CMC01|EAC-2010|KSL-2016",
  "",
  tax_tab_COI_fixed$scientificName
) |>
  trimws()

# get list of taxa, remove the taxon ranks (eg "Pleuronectidae family") from the species names
ncbi_taxa_COI <- as.data.frame(tax_tab_COI_fixed) |>
  pull(scientificName)

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
## 1   149045  Nitzschia      A.H. Hassall, 1845    accepted
## 1   111022  Amathia        Lamouroux, 1812       accepted
## 2   160576  Chlorella      M.Beijerinck, 1890    accepted
worms_taxonomy_tbl_COI <- map_dfr(unique_species_COI, get_worms_lineage)

# join lineage info back to ASVs
worms_df_COI <- tax_lookup_tbl_COI |>
  left_join(worms_taxonomy_tbl_COI, by = "Species") |>
  column_to_rownames("ASV") |>
  rename(scientificNameID = AphiaID, taxonRank = Rank)

# Append the verbatimIdentification and identificationQualifier, merging based on row names (ASV)
worms_tax_tab_COI_full <- worms_df_COI |>
  rownames_to_column("ASV") |>
  full_join(
    verbatimIdentification |> rownames_to_column("ASV"),
    by = "ASV"
  ) |>
  full_join(
    identificationQualifier |> rownames_to_column("ASV"),
    by = "ASV"
  ) |>
  column_to_rownames("ASV")

# Grab authority data from the AphiaID
worms_tax_tab_COI_full$aphiaID <- as.numeric(sub(
  "urn:lsid:marinespecies.org:taxname:",
  "",
  worms_tax_tab_COI_full$scientificNameID
))

worms_tax_tab_COI_full <- worms_tax_tab_COI_full |>
  mutate(
    scientificNameAuthorship = sapply(
      worms_tax_tab_COI_full$aphiaID,
      get_authority
    )
  )

worms_tax_tab_COI_full <- worms_tax_tab_COI_full |>
  select(
    kingdom = Kingdom,
    phylum = Phylum,
    class = Class,
    order = Order,
    family = Family,
    genus = Genus,
    scientificName = Species,
    scientificNameID,
    taxonRank,
    verbatimIdentification,
    aphiaID,
    scientificNameAuthorship
  )

# Upon inspection, there are 8 unique taxa for which there is no associated WoRMS LSID.
# Three of these taxa can be manually found on WoRMS, and their AphiaIDs have been verified with Libby Natola
# Their taxonomic information, lsid, rank and scientificNameAuthorship are created below and then merged
# into the worms_tax_tab_COI_full data.

names <- c(
  "Nitzschia inconspicua",
  "Prorocentrum pervagatum",
  "Synchaeta kitina"
)
records <- worrms::wm_records_names(name = names, marine_only = F)

# Extract AphiaID, authority, etc.
flat_info <- lapply(seq_along(records), function(i) {
  rec <- records[[i]]
  name <- names[i]

  if (!is.null(rec) && nrow(rec) > 0) {
    tibble(
      scientificName = name,
      kingdom = rec$kingdom[1],
      phylum = rec$phylum[1],
      class = rec$class[1],
      order = rec$order[1],
      family = rec$family[1],
      genus = rec$genus[1],
      aphiaID = rec$AphiaID[1],
      scientificNameAuthorship = rec$authority[1],
      taxonRank = rec$rank[1],
    )
  } else {
    tibble(
      scientificName = name,
      scientificNameID = NA,
      kingdom = NA,
      phylum = NA,
      class = NA,
      order = NA,
      family = NA,
      genus = NA,
      aphiaID = NA,
      scientificNameAuthorship = NA,
      taxonRank = NA
    )
  }
})

df <- bind_rows(flat_info) |>
  mutate(
    scientificNameID = paste0("urn:lsid:marinespecies.org:taxname:", aphiaID),
    verbatimIdentification = scientificName
  )
df <- df[c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "scientificName",
  "scientificNameID",
  "taxonRank",
  "verbatimIdentification",
  "aphiaID",
  "scientificNameAuthorship"
)]

tax_Nitzschia <- df |> filter(scientificName == "Nitzschia inconspicua")
tax_Proro <- df |> filter(scientificName == "Prorocentrum pervagatum")
tax_kitina <- df |> filter(scientificName == "Synchaeta kitina")

# Add into data frame:
worms_tax_tab_COI_full[
  worms_tax_tab_COI_full$scientificName == "Nitzschia inconspicua",
] <- tax_Nitzschia
worms_tax_tab_COI_full[
  worms_tax_tab_COI_full$scientificName == "Prorocentrum pervagatum",
] <- tax_Proro
worms_tax_tab_COI_full[
  worms_tax_tab_COI_full$scientificName == "Synchaeta kitina",
] <- tax_kitina

# Add identificationQualifier column from the tax_tab_COI_fixed df:
worms_tax_tab_COI_full <- merge(
  worms_tax_tab_COI_full,
  identificationQualifier,
  by = "row.names",
  all = TRUE
)

# Still a few AphiaIDs missing -
# !!! Verified with Matt Lemay that these taxa can be omitted from submission to OBIS:
worms_tax_tab_COI_full <- worms_tax_tab_COI_full |>
  filter(
    !scientificName %in%
      c(
        "Eukaryota",
        "Picobiliphyte",
        "Lagenidium caudatum",
        "Austrarchaea",
        "Trieres chinensis"
      )
  ) |>
  select(-aphiaID)

# Change rownames to column, remove NAs"
worms_tax_tab_COI_full <- rownames_to_column(
  worms_tax_tab_COI_full,
  var = "ASV"
)
worms_tax_tab_COI_full <- worms_tax_tab_COI_full |>
  mutate(across(everything(), ~ replace(.x, is.na(.x), "")))

# Remove row names and the first column, and name the ASV column appropriately:
worms_tax_tab_COI_full <- worms_tax_tab_COI_full[, -1]
rownames(worms_tax_tab_COI_full) <- NULL
colnames(worms_tax_tab_COI_full)[
  colnames(worms_tax_tab_COI_full) == "Row.names"
] <- "ASV"

# Add scientificNameAuthorship to scientificName. Trimws as well:
worms_tax_tab_COI_full <- worms_tax_tab_COI_full |>
  mutate(scientificName = paste(scientificName, scientificNameAuthorship))
worms_tax_tab_COI_full$scientificName <- trimws(
  worms_tax_tab_COI_full$scientificName
)

write.csv(
  worms_tax_tab_COI_full,
  here("data", "metabarcoding", "COI", "FC_taxonomy_COI_worms.csv")
)
