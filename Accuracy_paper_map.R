library(sf)
library(maps)
library(ggplot2)
library(dplyr)

# Getting the list of lystae-homogenate sample pairs and their location
M <- read.delim( "~/Desktop/git/paper-metabar-accuracy/matched_lysate_homogenate.tsv")
samples_list<-unique(M$sample)
samples_list

# Download to your local computer the figshare IBA repo and get metadata file"
metadata_path <- "~/Desktop/git/figshare-repos/iba/raw_data/v6/"
#Get the metadata file with all 4K + sample 
Samp_meta <- read.delim(paste0(metadata_path,"samples_metadata_malaise_SE.tsv"))
#Get the trap/sites metadata (geographic location)
Trap_meta <- read.delim(paste0(metadata_path,"sites_metadata_SE.tsv"))

# match all 856 samples (that have both lysates and homogenates) with the location (trap_ID)
## 1) Subset big sample_metadata file to leave only those lysrae-homogenate sample pairs"
subset_samples <- Samp_meta[Samp_meta$sampleID_FIELD %in% samples_list, ]

## 2) Extract unique Trap_IDs for your samples
trap_ids <- unique(subset_samples$trapID)

## 3) Get coordinates for these Trap_IDs
subset_coords <- Trap_meta[Trap_meta$trapID %in% trap_ids, ] 

## 4) Merge coordinates back with subset_samples by Trap_ID
lys_hom_samples <- merge(subset_samples, subset_coords, by = "trapID", all.x = TRUE)

### now add a special status to "15 samples" locations.
samp_15_meta<-read.delim("~/Desktop/git/paper-metabar-accuracy/15_samples_metadata.csv", sep = ";")
samp_15_traps <- unique(samp_15_meta$trap_ID)

# Filter to only those trap IDs, then keep only unique Trap_ID rows
smaller_df <- samp_15_meta[samp_15_meta$trap_ID %in% samp_15_traps, ] %>%
  dplyr::distinct(trap_ID, .keep_all = TRUE)
smaller_df$sample<-NULL
smaller_df$trapID <- smaller_df$trap_ID

subset_coords <- subset_coords %>%
  mutate(
    status = if_else(trapID %in% smaller_df$trapID, "Selected", "All samples")
  )


### Plot a map:

# Sweden outline
sweden_map <- maps::map("world", "Sweden", fill = TRUE, plot = FALSE)
sweden_sf <- st_as_sf(sweden_map)

# Convert your data to sf
points_sf <- st_as_sf(
  subset_coords,
  coords = c("longitude_WGS84", "latitude_WGS84"),
  crs = 4326
)

# Plot
map<-ggplot() +
  geom_sf(data = sweden_sf, fill = "snow", color = "gray40") +
  geom_sf(
    data = points_sf,
    aes(shape = status, fill = trap_habitat, color = trap_habitat),
    size = 3
  ) +
  coord_sf(xlim = c(10, 25), ylim = c(55, 70), expand = FALSE) +
  scale_color_manual(values = c("#DBE8D3", "#FAD510", "#02401B", "#D4E853", "#DD5759", "#46ACC8")) +           # good for up to 8 colors
  scale_shape_manual(values = c(21, 9)) + 
  scale_fill_manual(values=c("#DBE8D3", "#FAD510", "#02401B", "#D4E853", "#DD5759", "#46ACC8")) +
  theme_minimal(base_size = 14)

ggsave(file="Fig_1A_Map_Sampling_sites.jpg", height=7, width=6, plot = map)
 

