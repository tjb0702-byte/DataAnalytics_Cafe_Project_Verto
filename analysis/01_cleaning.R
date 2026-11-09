# ===================================================================
# 01_cleaning.R -- reads analysis/raw_data.csv, writes analysis/dataset.csv
# Run from the repository root
# ===================================================================

## CHECK WE ARE IN THE REPOSITORY ROOT
if (!file.exists("analysis/raw_data.csv")) {
  stop("Not in the repository root, or analysis/raw_data.csv is missing.")
}

library(dplyr)

dir_dat <- "analysis"


# ============================================================
# 1. LOAD AND INSPECT
# ============================================================
df_raw <- read.csv(file.path(dir_dat, "raw_data.csv"), sep = ",", header = TRUE)

glimpse(df_raw)
nrow(df_raw)
colSums(is.na(df_raw))
sum(duplicated(df_raw))


# ============================================================
# 2. SHORTER COLUMN NAMES AND A ROW NUMBER
# ============================================================
# column 23 has accented letters in its name, so rename it by position
names(df_raw)[23] <- "distance_raw"

# sheet_row is the row number in the spreadsheet (+1 for the header row)
df <- df_raw %>%
  mutate(sheet_row = row_number() + 1) %>%
  rename(
    district_raw  = District..number.,
    cafe_name     = Cafe.name,
    address       = Cafe.address..Google.Maps.,
    coordinates   = Coordinates..GMaps...Right.click.the.place...click.the.coordinates.,
    seats_inside  = Seating.Places.Inside,
    rating        = Rating..Google.,
    n_ratings_raw = Number.of.ratings..Google.,
    cafe_type_raw = Cafe.Type,
    is_roastery   = Is.Roastery..0.1.
  )


# ============================================================
# 3. REMOVE EMPTY ROWS
# ============================================================
# DECISION: a row with no cafe name is an empty row of the sheet, not a cafe.

df %>% filter(cafe_name == "") %>% count(Payment, distance_raw, Latitude, Longitude)

df <- df %>% filter(cafe_name != "")
nrow(df)


# ============================================================
# 4. REMOVE DUPLICATE CAFES
# ============================================================
# These pairs are the same business at the same address.
# EMA branches and DECADA/DECAGRAM are different cafes, so they stay.
#
# DECISION: keep the row entered first in the sheet. The rule is fixed in
# advance and does not depend on the values of any variable analysed later.

rows_in_a_pair <- c(3, 145, 10, 104, 45, 98, 48, 159, 50, 160, 58, 197,
                    95, 143, 103, 199, 107, 201, 111, 198, 112, 132)
rows_to_drop   <- c(145, 104, 98, 159, 160, 197, 143, 199, 201, 198, 132)

df %>%
  filter(sheet_row %in% rows_in_a_pair) %>%
  arrange(address) %>%
  select(sheet_row, cafe_name, address, cafe_type_raw, is_roastery, seats_inside, n_ratings_raw)

df <- df %>% filter(!sheet_row %in% rows_to_drop)
nrow(df)


# ============================================================
# 5. DISTRICT
# ============================================================
# Every spelling contains the district number, so keep only the digits.
#
# DECISION: rows 37, 39, 43, 44 are labelled Prague 2 but have Prague 1 or
# Mala Strana addresses. Row 144 is labelled Prague 7 but its address is Praha 5.

df %>% count(district_raw)

df <- df %>%
  mutate(
    district_number = gsub("[^0-9]", "", district_raw),
    district_label  = paste("Prague", district_number),
    district = case_when(
      sheet_row %in% c(37, 39, 43, 44) ~ "Prague 1",
      sheet_row == 144                 ~ "Prague 5",
      TRUE                             ~ district_label
    )
  )

df %>% count(district_raw, district)


# ============================================================
# 6. CAFE TYPE
# ============================================================
# DECISION: a blank type (row 136) is missing, not a category.
# DECISION: row 113 (Paul) is recoded because Paul is one international
# business and row 189 (another Paul branch) is labelled as such.
# Savoy and Cafe-Cafe are resolved by the duplicate rule in section 4.

df %>% count(cafe_type_raw)

df <- df %>%
  mutate(
    cafe_type = case_when(
      cafe_type_raw == "" ~ NA_character_,
      sheet_row == 113    ~ "International Chain",
      TRUE                ~ cafe_type_raw
    )
  )

df %>% count(cafe_type_raw, cafe_type)


# ============================================================
# 7. TEXT TO NUMBERS
# ============================================================
# First list everything that is not a digit or a decimal point
df %>%
  filter(grepl("[^0-9.]", n_ratings_raw) | n_ratings_raw == "") %>%
  select(sheet_row, cafe_name, n_ratings_raw)

df %>%
  filter(grepl("[^0-9.]", distance_raw) | distance_raw == "") %>%
  select(sheet_row, cafe_name, distance_raw)

# Then remove the thousands comma and convert
df <- df %>%
  mutate(n_ratings_no_comma = gsub(",", "", n_ratings_raw)) %>%
  mutate(
    n_ratings   = as.numeric(n_ratings_no_comma),
    distance_km = as.numeric(distance_raw)
  )

# Print the same rows again to check the conversion
df %>%
  filter(grepl("[^0-9.]", n_ratings_raw) | n_ratings_raw == "") %>%
  select(sheet_row, cafe_name, n_ratings_raw, n_ratings)

# DECISION: row 85 keeps 0 ratings. Its Google rating is NA, which is what a
# listing with no reviews looks like, so 0 is a real count and not a gap.
df %>% filter(sheet_row == 85) %>% select(sheet_row, cafe_name, rating, n_ratings)


# ============================================================
# 8. COORDINATES
# ============================================================
# DECISION: latitude and longitude come from the coordinates cell itself,
# because the Latitude/Longitude columns disagree with it in row 39.

df <- df %>%
  mutate(
    lat_text = sub(",.*", "", coordinates),  # everything before the comma
    lon_text = sub(".*,", "", coordinates),  # everything after the comma
    lat = as.numeric(lat_text),
    lon = as.numeric(lon_text)
  )

# Check: which rows disagree with the sheet's own Latitude/Longitude columns?
df %>%
  mutate(
    lat_gap = abs(lat - as.numeric(Latitude)),
    lon_gap = abs(lon - as.numeric(Longitude))
  ) %>%
  filter(lat_gap > 0.0001 | lon_gap > 0.0001) %>%
  select(sheet_row, cafe_name, coordinates, Latitude, Longitude, lat, lon)


# ============================================================
# 9. NEW VARIABLES
# ============================================================
# Formulas from the project guidelines. Positive east_west_km means east of
# the Orloj, positive north_south_km means north of it.
df <- df %>%
  mutate(
    east_west_km    = (lon - 14.4207) * 71,
    north_south_km  = (lat - 50.0871) * 111,
    east_west_abs   = abs(east_west_km),
    north_south_abs = abs(north_south_km)
  )

# The median is taken after duplicates are removed
median_ratings <- median(df$n_ratings, na.rm = TRUE)

# 0/1 indicators. Missing values stay missing. A cafe exactly on the median
# counts as not above it.
df <- df %>%
  mutate(
    is_east         = ifelse(east_west_km > 0, 1L, 0L),
    is_prague1      = ifelse(district == "Prague 1", 1L, 0L),
    is_intl_chain   = ifelse(cafe_type == "International Chain", 1L, 0L),
    is_above_median = ifelse(n_ratings > median_ratings, 1L, 0L)
  )


# ============================================================
# 10. CHECK AND EXPORT
# ============================================================
df_export <- df %>%
  select(
    sheet_row, cafe_name, district, cafe_type, is_roastery, seats_inside,
    n_ratings, distance_km, lat, lon, east_west_km, north_south_km,
    east_west_abs, north_south_abs, is_east, is_prague1, is_intl_chain,
    is_above_median
  )

nrow(df_raw)
nrow(df_export)
colSums(is.na(df_export))
summary(df_export)

write.csv(df_export, file.path(dir_dat, "dataset.csv"), row.names = FALSE)
