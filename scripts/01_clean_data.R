#### 01_clean_data.R ####
# Reads the raw greenhouse data, tidies it and saves a clean copy.
#   Input:  data/raw/Foa4_gh_data.csv
#   Output: data/clean/Foa4_gh_clean.csv
# Run with the RStudio project open, so the file paths work on any computer.

# Read the raw file. fileEncoding removes an invisible marker that Excel adds
# at the start of CSV files, which can break the first column name.
raw <- read.csv("data/raw/Foa4_gh_data.csv", fileEncoding = "UTF-8-BOM")

# Keep the columns we need, with clear names and units.
# - The file has two columns called "Rep"; R renames the second "Rep.1".
#   Rep = plant number in its group (R1-R26), Rep.1 = position in the block.
# - TreatCode, PlantVar and PlantVarTreat are dropped: they only repeat
#   Plant, Variety and Treatment.
clean <- data.frame(
  crop                = raw$Plant,
  variety             = raw$Variety,
  treatment           = ifelse(raw$Treatment == "Control", "Control", "Inoculated"),
  block               = raw$Block,
  position            = raw$Rep.1,
  plant_no            = as.integer(sub("R", "", raw$Rep)),   # "R7" becomes 7
  height_41dai_cm     = raw$height_middle,
  height_final_cm     = raw$height_final,
  total_weight_g      = raw$Total_weight,
  root_weight_g       = raw$Roots_weight,
  discoloration_pct   = raw$Perc_discoloration,
  discoloration_score = raw$Rate_discoloration   # 0 = died without symptoms, 1 = no symptoms
)

# Dead plants: final height AND weight both 0. Weight alone is not enough,
# because small living plants were recorded as 0 g. NA if either is missing.
clean$dead <- clean$height_final_cm == 0 & clean$total_weight_g == 0

# Checks: the script stops with an error if any of these is FALSE.
stopifnot(
  nrow(clean) == 208,                                   # all plants kept
  all(table(clean$variety, clean$treatment) == 26),     # 26 plants per group
  !any(duplicated(clean[, c("variety", "treatment", "plant_no")]))   # one row per plant
)

# Save (row.names = FALSE leaves out R's row numbers).
dir.create("data/clean", showWarnings = FALSE)
write.csv(clean, "data/clean/Foa4_gh_clean.csv", row.names = FALSE)

