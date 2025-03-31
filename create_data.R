library(DBI)
library(dbplyr)
library(dplyr)
library(RPostgres)
library(sf)
library(arrow)
library(geoarrow)

tryCatch(
  {
    con <- dbConnect(
      RPostgres::Postgres(),
      dbname = 'pep',
      host = Sys.getenv('PEP_PG_IP'),
      user = keyringr::get_kc_account("pgpep_londonj"),
      password = keyringr::decrypt_kc_pw("pgpep_londonj")
    )
  },
  error = function(cond) {
    print("Unable to connect to Database.")
  }
)

tbl_speno <- dplyr::tbl(con, in_schema("capture", "for_telem"))

tbl_deploy <- dplyr::tbl(con, in_schema("telem", "tbl_tag_deployments")) |>
  dplyr::filter(species %in% c('Pv', 'Pl')) |>
  dplyr::mutate(
    species = case_when(
      species == 'Pv' ~ 'Harbor seal',
      species == 'Pl' ~ 'Spotted seal',
      .default = NA
    )
  ) |>
  dplyr::left_join(tbl_speno, by = c("speno", "species")) |>
  dplyr::collect()

locs_obs <- sf::st_read(con, Id("telem", "geo_wc_locs_qa")) |>
  dplyr::left_join(tbl_deploy) |>
  dplyr::select(
    speno,
    deployid,
    ptt,
    instr,
    tag_family,
    type,
    quality,
    locs_dt,
    latitude,
    longitude,
    error_radius,
    error_semi_major_axis,
    error_semi_minor_axis,
    error_ellipse_orientation,
    project,
    species,
    age,
    sex,
    qa_status,
    deploy_dt,
    end_dt,
    deploy_lat,
    deploy_long,
    capture_lat,
    capture_long,
    geom
  ) |>
  dplyr::filter(deployid %in% c('PV2015_1008_14A0822', 'PL2018_1002_17A0776'))

locs_obs <- locs_obs |>
  mutate(x = st_coordinates(locs_obs)[, 1]) |>
  filter(x < -155 | x > 100) |>
  select(-x)

unlink(here::here('locs_obs.parquet'), force = TRUE)

write_parquet(locs_obs, here::here('seal_locs_obs.parquet'))

dbDisconnect(con, disconnect = TRUE)
