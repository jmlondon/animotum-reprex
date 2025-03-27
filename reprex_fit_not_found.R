library(dplyr)
library(sf)
library(geoarrow)
library(arrow)
library(aniMotum)

# PL2018_1002_17A0776

# obs |>
#   mutate(x = st_coordinates(obs)[, 1]
#         ) |>
#   filter(x < -155 | x > 100) |>
#   select(-x) |>
#   tibble::as_tibble() |>
#   write_parquet('spotted_seal.parquet')

locs_obs <- open_dataset('spotted_seal.parquet') |>
  sf::st_as_sf() |>
  sf::st_transform(3571)

locs_fit <- aniMotum::fit_ssm(
  x = locs_obs,
  vmax = 8,
  model = "crw",
  time.step = 0.25,
  id = "speno",
  date = "datetime",
  lc = "quality",
  epar = c(
    "error_semi_major_axis",
    "error_semi_minor_axis",
    "error_ellipse_orientation"
  ),
  map = list(psi = factor(NA)),
  tz = "UTC"
)

predict_pts <-
  aniMotum::grab(locs_fit, what = "predicted", as_sf = TRUE, group = TRUE) |>
  dplyr::rename(speno = id, datetime = date)

predict_lines <- predict_pts %>%
  dplyr::group_by(speno) %>%
  dplyr::summarise(do_union = FALSE) %>%
  sf::st_cast("LINESTRING")

library(ggplot2)

ggplot() + geom_sf(data = locs_obs) + geom_sf(data = predict_lines)

# now, let's simulate
sim_tracks <- sim_fit(x = locs_fit, what = "predicted", reps = 2)

# Error in `sim_fit()`:
# ! object 'fit' not found
# Hide Traceback
#     ▆
#  1. └─aniMotum::sim_fit(x = locs_fit, what = "predicted", reps = 2)
#  2.   └─base::sapply(fit$ssm, function(x) x$pm)
#  3.     └─base::lapply(X = X, FUN = FUN, ...)
