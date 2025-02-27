library(rvc)
library(tidyverse)

# Current list of calibrated species in the Caribbean 2/25/2025 ---------------------------------
calibrated_species_list <- c(
  "BAL VETU",
  "OCY CHRY",
  "SPA VIRI",
  "EPI GUTT",
  "LAC MAXI",
  "ACA COER",
  "SCA ISER",
  "STE PART",
  "CHA CAPI",
  "CAN ROST",
  "STE VARI",
  "HAL GARN",
  "HOL RUFU",
  "SCA TAEN",
  "HOL ADSC",
  "ACA CHIR",
  "HAL MACU",
  "CAR RUBE",
  "PSE MACU",
  "STE PLAN",
  "CEP FULV",
  "HAE FLAV",
  "CEP CRUE",
  "HAE PLUM",
  "LUT APOD",
  "AUL MACU",
  "LUT MAHO",
  "SER TIGR",
  "LUT SYNA",
  "SPA CHRY",
  "HAL BIVI",
  "ELA EVEL",
  "STE ADUS",
  "COR GLAU",
  "HYP UNIC",
  "SPA ATOM",
  "COR PERS",
  "CHR CYAN",
  "HYP PUEL",
  "CHA STRI",
  "ABU SAXA",
  "HAL POEY",
  "SPA RUBR",
  "CHR MULT",
  "STE LEUC",
  "CAN PULL",
  "GNA THOM",
  "MIC CHRY",
  "HAE AURO",
  "HAE SCIU",
  "BOD RUFU",
  "GRA LORE",
  "HOL TRIC",
  "HYP CHLO",
  "ANI VIRG",
  "MUL MART",
  "MYR JACO",
  "HAL RADI",
  "MAL TRIA",
  "HOL CILI",
  "SPH BARR",
  "POM ARCU",
  "SCA VETU",
  "OPI AURI",
  "SYN INTE",
  "LAC TRIQ",
  "OPH MACC",
  "POM PARU",
  "SCO REGA",
  "CLE PARR",
  "SER TABA",
  "PTE VOLI",
  "MEL NIGE",
  "HAL PICT",
  "STE DIEN",
  "NEO MARI",
  "CAL CALA",
  "ACA BAHI",
  "SPA AURO",
  "THA BIFA"
)

# ----------------------------------

SEDAR_triage_occurrence <- function(df, spp_list) {
  
  # subset spp_list to only those species we currently have a calibration factor
  spp <- spp_list %>% filter(SPECIES_CD %in% calibrated_species_list) %>% pull(SPECIES_CD)
  
  # subset further to only species that at least 1% occurrence in any of the NCRMP years 2016 to present
  a <- getDomainOccurrence(df, species = spp, years = c(seq(2016, format(Sys.Date(), "%Y"), 1))) %>%
    select(YEAR:occurrence) %>% 
    arrange(REGION, SPECIES_CD, YEAR) %>% 
    pivot_wider(id_cols = REGION:SPECIES_CD, names_from = YEAR, values_from = occurrence ) %>% 
    replace(., is.na(.), 0) %>% 
  # Change the occurrence cuttoff here
    filter(if_any(starts_with('20'), function(x) x >= 0.01))
  
  # return final list of species that we will provide survey estimates          
  return(df$taxonomic_data %>% filter(SPECIES_CD %in% a$SPECIES_CD) %>% select(SPECIES_CD:COMNAME))
}

SEDAR_triage_output <- function(data, spp) {
  
  dens <- getDomainDensity(data, spp) %>%
    mutate(density_SE = sqrt(var), density_CV = (sqrt(var))/density) %>% 
    select(YEAR:density, n, density_SE, density_CV)
  
  occ <- getDomainOccurrence(data, spp) %>%
    mutate(occurrence_SE = sqrt(var), occurrence_CV = (sqrt(var))/occurrence) %>% 
    select(YEAR:occurrence, n, occurrence_SE, occurrence_CV)
  
  observations <- data$sample_data %>%
    filter(SPECIES_CD %in% spp) %>% 
    group_by(YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, SPECIES_CD) %>% 
    summarise(num_observations = sum(NUM), .groups = "drop") %>% 
    mutate(num_observations = if_else(YEAR < 2016,
                                      ceiling(num_observations),
                                      ceiling(num_observations * 2))) %>% 
    group_by(YEAR, SPECIES_CD) %>% 
    summarise(num_observations = sum(num_observations), .groups = "drop")
  

  return(Reduce(function(x, y) merge(x, y, all=FALSE), list(dens, occ, observations)) %>% arrange(REGION, SPECIES_CD))

  
}

binned_LenFreq = function(d, spp, bin_size, yrs = NULL, st = NULL,colName) {
  
  a = getDomainLengthFrequency(d, species = spp, years = yrs, status = st )
  b = merge(a, data.frame(length_class = seq(1,max(a$length_class),0.5)), all.y =T) %>%
    select(length_class, frequency) %>%
    replace(., is.na(.), 0) %>%
    mutate(bin= as.numeric(cut(length_class, seq(0,max(length_class) + 5,bin_size))))
  c <- b %>%
    group_by(bin) %>%
    summarise(freq = sum(frequency))
  
  colnames(c)[2] <- colName
  
  return(c)
  
}

lenFreq_by_year <- function(df, spp, bin_size) {
  yearList = unique(df$sample_data$YEAR)
  
  dataList <- list()
  for(i in yearList){

    try({a <- binned_LenFreq(d = df, spp = spp, yrs = i, bin_size = bin_size, colName = i)
    dataList[[as.character(i)]] <- a}
    ,silent = TRUE)
  }

  l <- Reduce( function(x, y, ...) merge(x, y, all = TRUE, ...), dataList ) %>%
    replace(., is.na(.), 0) %>% 
    pivot_longer(cols = !c(bin), names_to = "YEAR", values_to = "frequency") %>% 
    mutate(SPECIES_CD = toupper(spp)) %>% 
    arrange(YEAR) %>% 
    select(YEAR, SPECIES_CD, bin, frequency)
  
  return(l)
  
}

multi_spp_lenfreq_by_year <- function(df, spplist, bin_size) {
  
  dataList <- list()
  for(spp in spplist){
    print(paste0("begin lf of ", spp))
    a <- lenFreq_by_year(df = df, spp = spp, bin_size = bin_size)
    dataList[[as.character(spp)]] <- a
    print(paste0(spp, " successfully created length frequency"))
  }
  
  return(Reduce( function(x, y, ...) merge(x, y, all = TRUE, ...), dataList ) %>% arrange(SPECIES_CD, YEAR, bin))
}


