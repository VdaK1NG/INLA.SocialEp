library("spdep") # Manejar datos espaciales

us_county_23 <- read_sf('./data-raw/cb_2023_us_county_5m.shp')
temp_data_us_23 <- read.csv("./data-raw/data_temp.csv")
pop_data_us_23 <- read.csv("./data-raw/co-est2023-alldata.csv")
load("./data-raw/ind_ef.Rdata")
sim_sp_ef <- ind.ef

save(us_county_23, file="us_county_23.rda")
save(temp_data_us_23, file="temp_data_us_23.rda")
save(pop_data_us_23, file="pop_data_us_23.rda")
save(sim_sp_ef, file="sim_sp_ef.rda")

usethis::use_data_raw("us_county_23")
usethis::use_data_raw("temp_data_us_23")
usethis::use_data_raw("pop_data_us_23")
usethis::use_data_raw("sim_sp_ef")
