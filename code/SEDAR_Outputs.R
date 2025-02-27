
SEDAR_exploited_phase <- function(data, spp, lc, growthParams = NULL) {
  
  dens <- getDomainDensity(data, spp, length_bins = lc) |> 
    filter(length_class == paste(">=", lc)) |> 
    mutate(density_SE = sqrt(var))
  
  abun <- getDomainAbundance(data, spp, length_bins = lc) |> 
    filter(length_class == paste(">=", lc)) |> 
    mutate(abundance_SE = sqrt(var))
  
  bio <- getDomainTotalBiomass(data, spp, length_bins = lc, growth_parameters = growthParams) |> 
    filter(length_class == paste(">=", lc)) |> 
    mutate(biomass_SE = sqrt(var))
  
  lbar <- getDomainLbar(data, spp, length_bins = lc) |> 
    filter(length_class == paste(">=", lc)) |> 
    mutate(lbar_SE = sqrt(vbar_L))

  return(Reduce(function(x, y) merge(x, y, all=FALSE), list(dens, abun, bio, lbar)))

  
}