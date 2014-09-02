#' Compute multiple segregation indices.
#'
#' Compute different segregation indices between regions,
#'  for multiple counties, and for multiple minority.
#'
#' @param data data.frame, containig at least region ids/names and goups
#'   of population.
#' @param region name or number of the column which contains the ids of
#'   small regions, between which the segregation is computed.
#' @param county Ids/names of bigger regions, for which the given
#'   indices should be computed.
#' @param cent_lng longitude of the center of the region
#' @param cent_lat latitude of the center of the region
#' @param groups character or numeric vector with the names or
#'   numbers of the population groups.
#' @param indices character vector with the functions' names to compute 
#' @param centrum_mode Either "GeoCentrum" or "Precentrum", see Details.
#' @param cent_lng character or numeric
#'   the longitude for each county. The same values
#'   must be present for every region in a county.
#' @param cent_lat character or numeric
#'   the latitude for each county. The same values
#'   must be present for every region in a county.
#' @param mode Computation of self-distance,
#'   "classic" or "realistic", see Details.
#'
#'
#' @details
#'   This function is a wrapper, which will compute the given segregation
#'     indices between region for each county and for each group as 
#'     minority, compared to all the others as majority,
#'     If the county is missing, it will be assumed that every region
#'     is part of the same county.
#'
#'   centrum_mode refers to the way of handling the center of a county.
#'     GeoCentrum will compute the mean of region coordinates,
#'     PreCentrum will use the country_center variables. If in a county
#'       there are more than one center given, the first non-NA will be
#'       used and a warning will be thrown.
#'
#'   distance_mode defines the way of handling a region's distance from
#'     itself. The classic method is computed as sqrt(0.6 * area), while
#'     the realistic method takes the distance of the closest region.
#'
#' @return data.frame with the county-aggregated group numbers,
#'   values for each segregation indices.
#'
segr_inds <- function(data, region, county, groups,  indices, ...,
                      longitude, latitude, centrum_mode,
                      county_center, mode){

  # getting the given arguments
  dots <- list(...)
  if(!length(dots)){
    dots <- list()  
  }

  # checking if county is given. If not, give it:
  if(missing(county) | is.null(county)){
    # add a county variable, but check if it is already there:
    if("added_county_var" %in% names(data)){
      stop("could not add a temporary variable because a variable
            name 'added_county_variable'.\n\tPlease use another name.")
    } else {
      data$added_county_variable <- 0
      county <- as.name("added_county_variable")
    }
  }

  # create a total row if not given:
  if(is.null(dots$total)){
    # check if there's already a total
    if("added_total_variable" %in% names(data)){
      stop("Could not add total population variable.\n
           Please name it in '...' or change the name of your
           'added_total_variable' variable")
    } else {
      data$added_total_variable <- rowSums(data[groups], na.rm=TRUE)
      dots$total <- as.name("added_total_variable")
    }
  }

  # every argument must be a name to evaluate properly
  # if(length(dots)){
  #   dots <- sapply(dots, function(vars){lapply(vars, as.name)})
  # }

  # define the function call
  computing_funct <- expression(
    sapply(indices, function(ind){
      # create a quoted string
      quoted_ind <- paste0("\'", ind, "\'")
      # create a function call
      paste0("do.call(", quoted_ind,
             ", dots[names(dots) %in% names(formals(",
              ind,
              "))])")
    }))

  # the function for each county
  county_compute <- function(county_df){
    # check if county has a unique center (if center is given):
    if(!is.null(dots$cent_lng) | !is.null(dots$cent.lat)){
      if(nrow(unique(county_df[c(cent_lng, cent_lat)]))>1){
        warning(paste("County", unique(county_df[county]),
                      "has more than one centre, the first will be used!"))
        cents <- complete.cases(unique(county_df[
                                       c(cent_lng, cent_lat)]))
        county_df[cent_lng] <- cents[1,cent_lng]
        county_df[cent_lat] <- cents[1,cent_lat]
      }
    }
    
  groups <- sapply(groups, as.name)
  # compute all indices for each group
  ldply(groups, function(group){
    # comitting an evil sin:
    dots$x <- as.name(group)
    # evaluate computing function:
    d <- (lapply(eval(computing_funct), function(fcall){
                 parse(text=fcall)}))
    return(sapply(d, function(fcall){with(county_df, eval(fcall))}))
    },
    .parallel=TRUE)
  }

  # compute "county_compute" for each county
  all_indices <- ddply(data, c(county), county_compute,
                       .parallel=TRUE)
    return(all_indices)
}
