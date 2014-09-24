# Functions for computing segregation indices
# originally written by Gergely Tóth


#------------------------#
#                        #
#    HELPER FUNCTIONS    #
#                        #
#------------------------#


d_ij <- function(area, lng, lat, mode=c("classic", "realistic")){
  # Computing distances between regions

  # area: numeric vector of area for each region
  # lng: numeric vector of the center's longitude for each region
  # lat: numeric vector of the center's latitude for each region
  # mode: character constant for computing of self-distance.
  #  classic: self-distance is (0.6 * area)^0.5
  #  realistic: self-distance is equal to the nearest neighbour

  
  if(!(is.numeric(area) && is.numeric(lng) && is.numeric(lat))){
    stop("area, lng and lat must be numeric!")
  }
  if(!(length(area) == length(lng))){
    stop("area and lng must be the same length!")
  }
  if(!(length(lng) == length(lat))){
    stop("lat and lng must be the same length!")
  }
  if(!(length(area) == length(lat))){
    stop("lat and lng must be the same length!")
  }
  if(!(mode %in% c("classic", "realistic"))){
    stop("Unknown mode! Only classic and realistic is implemented!")
  }

  dij <- as.matrix(dist(data.frame(lng=lng, lat=lat)))

  if (mode == "classic") { 
    dii <- (0.6*area)^0.5
    diag(dij) <- dii 
    }
  if(mode == "realistic"){
    diag(dij) <- apply(dij, 1, min, na.rm =TRUE)
  }

  return(dij)
}


c_ij <- function(area, lng, lat, mode=c("classic", "realistic")){
  # same as exp(-d_ij(area, lng, lat, mode))

  # area: numeric vector of area for each region
  # lng: numeric vector of the center's longitude for each region
  # lat: numeric vector of the center's latitude for each region
  # mode: character constant for computing of self-distance.
  #  classic: self-distance is (0.6 * area)^0.5
  #  realistic: self-distance is equal to the nearest neighbour
    
  dij <- d_ij(area, lng, lat, mode)
  cij <- exp(-dij)
  return(cij)
}


P_gg <- function(per.area, TOTAL, cij, area, lng, lat,
                 mode=c("centrum", "realistic")){
  # per.area: numeric vector containing population size for each area
  # TOTAL: numeric, number of total population
  # cij: result of c_ij(). If not given, will be computed from rest.

  if(missing(cij)){
    cij <- c_ij(area, lng, lat, mode)
  }

  M <- matrix(per.area, ncol=length(per.area), nrow=length(per.area))
  res <- (M * t(M) * cij) / (TOTAL / 2)
  return(res)
}

cent_dist <- function(lng, lat, mode=c("GeoCentrum", "PreCentrum"),
                      cent_lng, cent_lat){
  # computes the distance from each region to the centrum.
  # lng, lat: longitude and latitude for each region
  # mode: GeoCentrum computes the mean of lng and lat,
  #       Precentrum uses predefined centrum coordinates
  # cent_lng, cent_lat: longitude & latitude for the predefined centrum

  if(mode == "PreCentrum" &&
    (!exists(cent_lng) | !exists(cent_lat))){
    mode <- "GeoCentrum"
    warning("For centrum distance, PreCentrum was given, although at
            least one of cent_lng or cent_lat was missing.\n
            \tCentrum is computed as in GeoCentrum mode!")
  }

  if(mode == "GeoCentrum"){
    cent_lng <- mean(lng)
    cent_lat <- mean(lat)
  }


  centr.dist <- sqrt((lng-cent_lng)^2 + (lat-cent_lat)^2)

  return(centr.dist)  
}


#-----------------------#
#                       #
#  SEGREGATION INDICES  #
#                       #
#-----------------------#


abscent_index <- function(x, area, lng, lat, cent_dist, incr_order_cent_dist,
                          centrum_mode=c("GeoCentrum", "PreCentrum"), cent_lng, cent_lat){
  # computes absolute centralization index

  # x, y: minority and majority population for each region
  # incr_order_cent_dist: rank of aech region in the increasing order of the distance from the centre.
  #   If incr_order mitted, will be computed from cent_dist
  # lng, lat: longitude and latitude for each region
  # cent_dist: distance from the centre. If not given, computed according to mode.
  # centrum_mode: type of centrum of the areas.
  #  GeoCentrum: center is the mean of lat and lng
  #  PreCentrum: the centrum is given in cent_lng, cent_lat

  if(missing(cent_dist)){
    cent_dist <- cent_dist(lng, lat, centrum_mode, cent_lng, cent_lat)
  }
  if(missing(incr_order_cent_dist)){
    incr_order_cent_dist <- order(cent_dist)
  }

  x <- x[incr_order_cent_dist]
  area <- area[incr_order_cent_dist]
  sum.C1 <- sum(mapply('*', x[1:length(x-1)], area[2:length(area)]))
  sum.C2 <- sum(mapply('*', x[2:length(x)], area[1:length(area-1)]))

  relcent_value <- sum.C1 - sum.C2

  return(abscent_value)
}


absclust_index <- function(x, X, total, cij, area, lng, lat,
                           mode=c("classic", "realistic")){
  # Computes absolute clustering index

  # x,t: numeric vector of minority and total population
  # X: number of minority people. If missing, sum of x will be used.
  # cij: closeness matrix. Ig not given, computed from are, lng, lat, mode
  # area, lng, lat: numeric vector of each region's area, center longitude and latitude
  # mode: "realistic" and "classic". Former defines self-distance az minimum distance from others,
  #       latter uses (0.6*area)^0.5 to estimate it.

  if(missing(cij)){
    cij <- c_ij(area, lng, lat, mode)
  }

  if(missing(X)){
    X <- sum(x)
  }

  masodik.tag_absC <- (X/(nrow(cij))^2)*sum(cij)
  szamlalo.elsotag__absC <- 0
  nevezo.elsotag__absC <- 0
  for (i_absClu in 1:length(x)) {
    szamlalo.elsotag__absC <- szamlalo.elsotag__absC + ((x[i_absClu]/X)*sum(cij[i_absClu,]*x))
    nevezo.elsotag__absC <- nevezo.elsotag__absC + ((x[i_absClu]/X)*sum(cij[i_absClu,]*total))
    }
  absclust_value <-(szamlalo.elsotag__absC - masodik.tag_absC)/(nevezo.elsotag__absC - masodik.tag_absC)

  return(absclust_value)
}

absconc_index <- function(x, X, total, area, incr_order, decr.order){
  # x,total, area: numeric vectors about minority and total population
  #   and area of each region
  # If X or TOTAL missing, they are computed as sums of x or t respectively.
  # incr.sort: numeric vector, if missing, computed by order(area)
  # decr.sort: numeric vector, if missing, computed by order(area, decreasing=TRUE)

  if(missing(X)){
    X <- sum(x, na.rm=TRUE)
  }
  if(missing(incr_order)){
    incr_order <- order(area)
  }
  if(missing(decr.order)){
    decr.order <- order(area, decreasing=TRUE)
  }
  

  total.sum.inc <- cumsum(total[incr_order])
  n1 <- sum(total.sum.inc < X) + 1
  t1 <- sum(total[incr_order][1:n1])

  total.sum.decr <- cumsum(total[decr.order])
  n2 <- sum(total.sum.decr < X) + 1
  t2 <- sum(total[decr.order][1:n1])




  szamlalo1 <- sum(((x*area)/X), na.rm =TRUE)
  szamlalo2 <- sum(((total[incr_order][1:n1]*area[incr_order][1:n1])/t1), na.rm =TRUE)
  nevezo1 <- sum(((total[decr.order][1:n2]*area[decr.order][1:n2])/t2), na.rm =TRUE)
  nevezo2 <- sum(((total[incr_order][1:n1]*area[incr_order][1:n1])/t1), na.rm =TRUE)
  absconc_value <- 1-((szamlalo1 - szamlalo2)/(nevezo1-nevezo2))
  
  return(absconc_value)
}

atkinson_index <- function(x, X, total, TOTAL, p, P, b){
  # Computes Atkinson index

  # x,t,p: numeric vector of minorityand total population, or proportion of
  #  minority to total population for each area. if p is given, x is not considered.
  # X,T,P: 1 legth numeric, the number of minority and total population, or proportion
  #  of minorities to total. X  can be omitted and it will be counted as sum of x.
  # b: parameter of atkinson index
  
  # check if attributes are given
  if(missing(p)){
    p <- x / total
  }
  if(missing(P)){
    if(missing(X)){
      X <- sum(x, na.rm=TRUE)
    }
    if(missing(TOTAL)){
      TOTAL <- sum(total, na.rm=TRUE)
    }
    P <- X / TOTAL
  }

  szorzat1 <- (P / (1-P))
  szorzat2 <- abs(sum((((1-p)^(1-b))*(p^b)*total)/(P*TOTAL), na.rm = TRUE))^(1/(1-b))
  atkinson_value  <- 1-(szorzat1*szorzat2)

  return(atkinson_value)
}

correlation_index <- function(x, X, total, TOTAL, P){
  # computes correlation index

  # x,t: numeric vector of minorityand total population,
  # X,P: 1 legth numeric, the number of minority and total population, or proportion
  #  of minorities to total. X  can be omitted and it will be counted as sum of x.
  
  # check if attributes are given
  
  if(missing(P)){
    if(missing(X)){
      X <- sum(x, na.rm=TRUE)
    }
    if(missing(TOTAL)){
      TOTAL <- sum(total, na.rm=TRUE)
    }
    P <- X / TOTAL
  }

correlation_value <- (isolation_index(x, X, total)-P)/(1-P)

return(correlation_value)
}

delta_index <- function(x, X, area, A){
  # Comutes Duncan's delta index

  # x, area: numeric vectors for minority population and area for each region
  # X and A can be omitted (sum of area and minority pop.), they will be computed.

  if(missing(X)){
    X <- sum(x)
  }

  if(missing(A)){
    A <- sum(area)
  }

  delta_value <- sum(abs((x / X) - (area / A)), na.rm = TRUE)*0.5
  return(delta_value)
}

dissimilarity_index <- function(x, X, total, TOTAL){
  # Computes dissimilarity index

  # x,total: numeric vectors of minority and total population
  # X,TOTAL: number of all minority and total population - they can be omitted.

  if(missing(X)){
    X <- sum(x, na.rm = TRUE)
  }
  if(missing(TOTAL)){
    TOTAL <- sum(total, na.rm = TRUE)
  }

  b_i_per_B <- x / X
  w_i_per_w <- (total-x) / (TOTAL - X)
  summazat <-  sum(abs(b_i_per_B - w_i_per_w) , na.rm = TRUE)

  dis_value <- summazat / 2

  return(dis_value)
}

distdecinteract_index <- function(x, total, dij, y, area, lng, lat,
                               mode=c("classic", "realistic")){
  # Computes distant decay interaction index

  # x, total, y: numeric vectors for each region for minority, total, majority population
  # dij: object returned from d_ij()- If not supplied, will be computed.

  if(missing(dij)){
    dij <- d_ij(area, lng, lat, mode)
  }
  if(missing(y)){
    y  <- total - x
  }
  if(missing(x)){
    x <- total - y
  }
  if(missing(total)){
    total <- x + y
  }


  kij <- total^(-dij)
  kij <- kij/colSums(kij)

  fst <- x/X
  snd <- colSums(t(kij)*y/t)
  distinteract_value <- fst*snd

  return(distinteract_value)
}

distdecisolation_index <- function(x, X, total, dij, area, lng, lat, 
                              mode=c("classic", "realistic")){
  # Computes distant decay isolation index

  # x, total: minority and total population for each region
  # X: number of all minority population, if not given, computed as sum(x)
  # dij: ersult of dij() function. If not given, computed from the others.
  # area, lng, lat: numeric vectors with each region's area, longitude and latitude center
  # mode: method for computing self-distance

  if(missing(X)){
    X <- sum(X)
  }
  if(missing(dij)){
    dij <- d_ij(area, lng, lat, mode)
  }

  kij <- total^(-dij)
  kij <- kij/colSums(kij)

  fst <- x/X
  snd <- colSums(t(kij)*x/t)
  distinteract_value <- fst*snd

  return(distinteract_value)
}

entropy_index <- function(x, X, total, TOTAL, p, P){
  # Computes entropy index

  # x,total,p: numeric vector of minorityand total population, or proportion of
  #  minority to total population for each area. if p is given, x is not considered.
  # X,TOTAL,P: 1 legth numeric, the number of minority and total population, or proportion
  #  of minorities to total. X  can be omitted and it will be counted as sum of x.
  
  # check if attributes are given
  if(missing(p)){
    p <- x / total
  }
  if(missing(P)){
    if(missing(X)){
      X <- sum(x, na.rm=TRUE)
    }
    if(missing(TOTAL)){
      TOTAL <- sum(TOTAL, na.rm=TRUE)
    }
    P <- X / TOTAL
  }
  
  E <- P*log(1/P)+(1-P)*log(1/(1-P))
  e <- p*log(1/p)+(1-p)*log(1/(1-p))
  e[!is.finite(e)] <- 0

  entropy_value <- sum((total * (E - e)) / (E * TOTAL), na.rm = TRUE)

  return(entropy_value)
}

gini_index <- function(x, total, TOTAL, p, P){
  # Computes Gini index

  # x, t: numeric vectors with minority and total population,
  # p: proportion of minority to total population. If omitted, will be computed.
  # TOTAL, P:number of ptotal population and proportion of total minority.
  # If they omitted, they computed from x and total

  if(missing(p)){
    p <- x / p
  }
  if(missing(TOTAL)){
    TOTAL <- sum(total, na.rm=TRUE)
  }
  if(missing(P)){
    P <- sum(x, na.rm=TRUE) / TOTAL
  }

  cel_matrixPOP <- matrix(as.vector(total), length(total), length(total) )
  cel_matrix_p_i <- matrix(as.vector(p), length(total), length(total))

  options(warn = -1)
  def_matrix <- matrix(c(1:(1+length(total))),length(total),length(total))
  def_matrix[upper.tri(def_matrix)] <-  def_matrix[upper.tri(def_matrix)]-1 
  options(warn = 0)

  p_i_M_p_j <- abs(cel_matrix_p_i[def_matrix] - cel_matrix_p_i )
  t_i_x_t_j <- cel_matrixPOP*cel_matrixPOP[def_matrix]

  szamlalo <- sum((p_i_M_p_j * t_i_x_t_j),  na.rm = TRUE )

  gini_value <- szamlalo / (2*TOTAL^2*P*(1-P))
  
  return(gini_value)
}

interact_index <- function(x, X, y, total){
  # Computes interaction index

  # x, y, total: numeric vectors for minority, amjority and total population for each region.
  # at least two of them must supplied, the third will be computed.
  # X: number of total minority population. If not supplied, will be computed from sum(x)

  if(missing(x)){
    x <- total - y
  }
  if(missing(y)){
    y <- total - x
  }
  if(missing(total)){
    total <- x + y
  }
  if(missing(X)){
    X <- sum(x)
  }

  interact_value <- sum((x / X) * (y / total), na.rm = TRUE)

  return(interact_value)
}

isolation_index  <- function(x, X, total){
  # Computes isolation Index

  # x,t: numeric vector of minorityand total population,
  # X,: 1 legth numeric, the number of minority population.
  #  X  can be omitted and it will be counted as sum of x.
  
  # check if attributes are given
  
  if(missing(X)){
    X <- sum(x, na.rm=TRUE)
  }
  
  isolation_value <- sum((x / X) * (x / total), na.rm = TRUE)

  return(isolation_value)
}

relcent_index <- function(x, y, lng, lat, cent_dist, incr_order_cent_dist,
                          centrum_mode=c("GeoCentrum", "PreCentrum"), cent_lng, cent_lat){
  # Computes relative centralization index

  # x, y: minority and majority population for each region
  # incr_order_cent_dist: rank of aech region in the increasing order of the distance from the centre.
  #   If incr_order mitted, will be computed from cent_dist
  # lng, lat: longitude and latitude for each region
  # cent_dist: distance from the centre. If not given, computed according to mode.
  # mode: type of centrum of the areas.
  #  GeoCentrum: center is the mean of lat and lng
  #  PreCentrum: the centrum is given in cent_lng, cent_lat

  if(missing(cent_dist)){
    cent_dist <- cent_dist(lng, lat, mode, cent_lng, cent_lat)
  }
  if(missing(incr_order_cent_dist)){
    incr_order_cent_dist <- order(cent_dist)
  }

  x <- x[incr_order_cent_dist]
  y <- y[incr_order_cent_dist]

  sum.C1 <- sum(mapply('*', x[1:length(x-1)], y[2:length(y)]))
  sum.C2 <- sum(mapply('*', x[2:length(x)], y[1:length(y-1)]))

  relcent_value <- sum.C1 - sum.C2

  return(relcent_value)
}

relclust_index <- function(x, X, y, Y, cij, area, lng, lat,
                           mode=c("classic", "realistic")){
  # computes relative clustering index
 
  # x,y: minority and majoriti population in each area
  # X, Y: sum of minority and majority population (it can be omitted)
  # cij: result of c_ij(). If omitted, computed from the rest

  if(missing(X)){
    X <- sum(x)
  }
  if(missing(Y)){
    Y <- sum(y)
  }
  if(missing(cij)){
    cij <- c_ij(area, lng, lat, mode)
  }

  relclust_value <- P_gg(x, X, cij) / P_gg(y, Y, cij)

  return(relclust_value)
}

relconc_index <- function(x, X, y, Y, area, total, TOTAL, incr_order, decr.order){
  # Computes relative concentration index

  # x, y, total: minority, majority and total population for each region
  #  at least two of them must be supplied.
  # area: area of each region
  # X, Y, TOTAL: number of minority, majority and total population.
  #  If not supplied, then will be computed from x, y and total
  #  at least two of them must be supplied
  # incr_order, decr.order: increasing and decreasing order of regions by area

  if(missing(x)){
    x <- total - y
  }
  if(missing(y)){
    y <- total - x
  }
  if(missing(total)){
    total <- x + y
  }
  if(missing(incr_order)){
    incr_order <- order(area)
  }
  if(missing(decr.order)){
    decr.order <- order(area, dicreasing = TRUE)
  }

  total.sum.inc <- cumsum(total[incr_order])
  n1 <- sum(total.sum.inc < X) + 1
  t1 <- sum(total[incr_order][1:n1])

  total.sum.decr <- cumsum(total[decr.order])
  n2 <- sum(total.sum.decr < X) + 1
  t2 <- sum(total[decr.order][1:n1])

  szamlalo <- (x*a/X) / (y*a/y) -1
  nevezo1 <- sum((total[incr_order][1:n1] * area[incr_order][1:n1]) / t1)
  nevezo2 <- sum((total[decr.order][1:n2] * area[decr.order][1:n2]) / t2)
  nevezo.rel <- (nevezo1/nevezo2) - 1
  relconc_value <- szamlalo / nevezo

  return(relconc_value)
}

spaprox_index <- function(x, X, y, Y, total, TOTAL, cij, area, lng, lat,
                          mode=c("classic", "realistic")){
  # Computes spatial proximity

  # x, y, total: minority, majority and total population for each region
  #  at least two of them must be given
  # X, Y, TOTAL: number of minority, majority and total population
  #  Can be omitted, will be computes as the sum of x, y, total respectively
  # cij: output of c_ij(). If omitted, area, lng, lat, mode must be supplied.
  # area, lng, lat: area, longitude and latitude for each region
  # mode: computation method of self-distance:
  #  classic: self-distance is (0.6*area)^0.5
  #  realistic: self-distance is equal to the distance of the nearest neighbour

  if(missing(x)){
    x <- total - y
  }
  if(missing(y)){
    y <- total - x
  }
  if(missing(total)){
    total <- x + y
  }
  if(missing(X)){
    X <- sum(x, na.rm=TRUE)
  }
  if(missing(Y)){
    Y <- sum(y, na.rm=TRUE)
  }
  if(missing(TOTAL)){
    TOTAL <- X + Y
  }
  if(missing(cij)){
    cij <- c_ij(area, lng, lat, mode)
  }

  spv1 <- X * P_gg(x, X, cij)
  spv2 <- Y * P_gg(y, Y, cij)
  spv3 <- TOTAL * P_gg(total, TOTAL, cij)
  spaprox_value <- (spv1 + spv2) / spv3

  return(spaprox_value)
}


#-----------------------#
#                       #
#      THE PROGRAM      #
#                       #
#-----------------------#


####
# Computing function
####

computeSegrIndices <- function(x, X, y, Y, total, TOTAL, p, P, area, A, 
                               lng, lat, cent_lat, cent_lng, cent_dist, cij, dij){
  # Computes segregation indices for a dataset with one minority and
  #  in one spatial unit, for example a city, based on its districts.

  # data: data.frame containing the necessary variables for each region
  # x, y, total: coloumn of minority majority and total population.
  #  can be a string or a number. At least two of them must be given.
  # X, Y, TOTAL: number of minority, majority and total population.
  #  If omitted, will be computed
  # area, lng, lat: coloumn of area, longitude and latitude for each region.
  # cent_lat, cent_lng: longitude and latitude of the center of the unit.
  #  If omitted, geographical center will be used
  # cij: value of c_ij(). If omitted, will be computed from the rest

  # it manipulates the exterior indices list.


  ##########################################
  ##                                      ##
  ## THIS FUNCTION FOR INTERNAL USE ONLY! ##
  ##       missing inputs will not        ##
  ##           be computed!               ##
  ##                                      ##

  ##TODO: Kénytelen leszek itt is hibakezelni...
  ##########################################

  segr.res <- list()

  if(is.na(indices$dissimilarity)){
    segr.res$dissimilarity <- dissimilarity_index(x, X, total)
  }
  message("dissimilarity DONE.\n")
  if(is.na(indices$gini)){
    segr.res$gini <- gini_index(x, total, TOTAL, p, P)
  }
  message("gini DONE.\n")
  if(is.na(indices$entropy)){
    segr.res$entropy <- entropy_index(x, X, total, TOTAL, p, P)
  }
  message("entropy DONE.\n")
  if(length(segrInds$atkinson$b)){
    atk.res <- rep(list(NA), length(segrInds$atkinson$b))
    names(atk.res) <- segrInds$atkinson$b
    for(b in segrInds$atkinson$b){
      atk.res[[as.character(b)]] <- atkinson_index(x, X, total, TOTAL, p, P, b)
    }
    segr.res$atkinson <- atk.res
  }
  message("atkinsons DONE.\n")
  if(is.na(indices$interaction)){
    segr.res$interaction <- interact_index(x, X, y, total)
  }
  message("interaction DONE.\n")
  if(is.na(indices$isolation)){
    segr.res$isolation <- isolation_index(x, X, total)
  }
  message("isolation DONE.\n")
  if(is.na(indices$correlation)){
    segr.res$correlation <- correlation_index(x, X, total, TOTAL, P)
  }
  message("correlation DONE.\n")
  if(is.na(indices$duncan.delta)){
    segr.res$duncan.delta <- delta_index(x, X, area, A)
  }
  message("duncan delta DONE.\n")
  if(is.na(indices$absolute.concentration)){
    segr.res$absolute.concentration <- absconc_index(x, X, total, area,
                                                     incr_order, decr.order)
  }
  message("absconc DONE.\n")
  if(is.na(indices$relative.concentration)){
    segr.res$relative.concentration <- relconc_index(x, X, y, Y, area,
                                                     total, TOTAL,
                                                     incr_order, decr.order)
  }
  message("relconc DONE.\n")
  if(is.na(indices$relative.centralization)){
    segr.res$relative.centralization <- relcent_index(x, y, lng, lat, cent_dist,
                                                      incr_order_cent_dist,
                                                      mode, cent_lng, cent_lat)
  }
  message("relcent DONE.\n")
  if(is.na(indices$absolute.centralization)){
    segr.res$absolute.centralization <- abscent_index(x, area, lng, lat,
                                                      cent_dist,
                                                      incr_order_cent_dist,
                                                      mode, cent_lng, cent_lat)
  }
  message("abscent DONE.\n")
  if(is.na(indices$absolute.clustering)){
    segr.res$absolute.clustering <- absclust_index(x, X, total, cij, area,
                                                   lng, lat, mode)
  }
  message("absclust DONE.\n")
  if(is.na(indices$relative.clustering)){
    segr.res$relative.clustering <- relclust_index(x, X, y, Y, cij, area,
                                                   lng, lat, mode)
  }
  message("relclust DONE.\n")
  if(is.na(indices$distance.decay.interaction)){
    segr.res$distance.decay.interaction <- distdecinteract_index(x, total, dij,
                                                                 y, area,
                                                                 lng, lat, mode)
  }
  message("ddint DONE.\n")
  if(is.na(indices$distance.decay.isolation)){
    segr.res$distance.decay.isolation <- distdecisolation_index(x, X, total,
                                                                dij, area,
                                                                lng, lat, mode)
  }
  message("ddisol DONE.\n")
  if(is.na(indices$spatial.proximity)){
    segres$spatial.proximity <- spaprox_index(x, X, y, Y, total, TOTAL, cij,
                                                area, lng, lat, mode)
  }
  message("spaprox DONE.\n")
  return(segr.res)
}


program <- function(dataset, file, output,
                    groups,
                    area,
                    longitude,
                    latitude,
                    areal.unit,
                    areal.center.lng,
                    areal.center.lat,
                    all.indices = TRUE,
                    evenness = FALSE,
                    exposure = FALSE,
                    concentration  = FALSE,
                    centralization = FALSE,
                    clustering     = FALSE,
                    b,
                    centrum.mode = c("GeoCentrum", "PreCentrum"),
                    distance.mode = c("classic", "realistic")){
  # computes segregation indices for multiple areal units and multiple groups.
  # five dimension of segregation can be managed from arguments,
  #  and it has an "all" switch, which owerwrites everything else.
  #  dimensions are defined in Massey, Denton: The Dimensions of Residential Segregation (1988)
  #  dimensions set true will compute the corresponding indices,
  #  but indices set TRUE in segrInds variable will be computed.
  # In the dataset, there should be areal.unit for each entry.
  #  If not, areal.center.lng and lat will be examined.
  #  If not given, the dataset will be assumed for one areal unit.
  #  If only areal.unit is given, the "GeoCentrum" option will be used where necessary.
  #  If areal.unit can be omitted if both areal.center.lng and arel.center.lat is given.
  #  For each unique pair of areal.center.lng and areal.center.lat there will be assumed one unique unit.
  # dataset is a data.frame, can be omitted if file is given.
  # File must be a path to a headered comma or semiclon separated txt (csv), or an SPSS sav file.
  # if output is given, a csv file will be written onto the path given.
  # goups: character or numeric vector giving the name or number of variables about groups of peoples in the population
  # area,longitude, latitude, areal.unit, areal.center.lng, areal.center.lat: character or numeric
  #  giving the name or number of region area, longitude, latitude, container areal unit, its longitude and latitude's variable
  # b: constant(s) for the atkinson index

  # dataset reading
  if(missing(dataset)){
    form <- unlist(strsplit(file))
    form <- paste0(form[(length(form)-2):length(form)], collapse="")
    if(form == "csv"){
      dataset <- read.csv(file, header=TRUE)
      if(ncol(dataset == 1)){
        dataset <- read.csv2(file, header=TRUE)
      }
    }
    if(form == "sav"){
      dataset <- foreign::read.spss(file, to.data.frame = TRUE)
    }
    if(!(form %in% c("csv", "sav"))){
      stop("File extension not recognized, please supply data in 'dataset'!")
    }
  }
  if(!is.data.frame(dataset)){
    dataset <- as.data.frame(dataset)
  }

  # setting indices to compute
  indices <- segrInds

  if(all){
    indices[!(names(indices) == "atkinson")] <- lapply(indices[!(names(indices) == "atkinson")],
                                                       function(ind){ind <- TRUE})
    if(missing(b)){
      warning("b missing: atkinson index will not be computed!")
    } else {
      indices$atkinson$b <- b
    }
  } else {
    if(evenness){
      indices[names(indices) %in% c("dissimilarity", "gini")] <- TRUE
      if(missing(b)){
        warning("b missing: atkinson index will not be computed!")
      } else {
        indices$atkinson$b <- b
      }
    }
    if(exposure){
      indices[c("interaction", "isolation", "correlation")] <- TRUE
    }
    if(concentration){
      indices[c("duncan.delta", "absolute.concentration",
                "relative.concentration")] <- TRUE
    }
    if(centralization){
      indices[c("absolute.centralization", "relative.centralization")] <- TRUE
    }
    if(clustering){
      indices[c("absolute.clustering", "relative.clustering",
                "distance.decay.interaction", "distance.decay.isolation",
                "spatial.proximity")] <- TRUE
    }
  }


  # converting variable names to variable numbers
  var.names <- dimnames(dataset)[[2]]
  if(is.character(groups)){
    n.groups <- sapply(groups, function(g){which(var.names == g)})
  } else {
    n.groups <- groups
  }
  if(!missing(area)){
    if(is.character(area)){
      n.area <- which(var.names == area)
    } else {
      n.area <- area
    }
  }
  if(!missing(longitude)){
    if(is.character(longitude)){
      n.longitude <- which(var.names == longitude)
    } else {
      n.longitude <- longitude
    }
  }
  if(!missing(latitude)){
    if(is.character(latitude)){
      n.latitude <- which(var.names == latitude)
    } else {
      n.latitude <- latitude
    }
  }
  if(!missing(areal.unit)){
    if(is.character(areal.unit)){
      n.areal.unit <- which(var.names == areal.unit)
    } else {
      n.areal.unit <- areal.unit
    }
  }
  if(!missing(areal.center.lng)){
    if(is.character(areal.center.lng)){
      n.areal.center.lng <- which(var.names == areal.center.lng)
    } else {
      n.areal.center.lng <- areal.center.lng
    }
  }
  if(!missing(areal.center.lat)){
    if(is.character(areal.center.lat)){
      n.areal.center.lat <- which(var.names == areal.center.lat)
    } else {
      n.areal.center.lat <- areal.center.lat
    }
  }
        
    

  # mapping areal units
  if(missing(areal.units)){
    if(missing(areal.center.lng) || missing(areal.center.lat)){
      n.areal.units <- NCOL(dataset) + 1
      dataset[n.areal.units] <- numeric(NROW(dataset))
      centrum.mode <- "GeoCentrum"
    } else {
      n.areal.units <- NCOL(dataset) + 1
      dataset[n.areal.units] <- cumsum(!duplicated(dataset[c(n.areal.center.lat,
                                                n.areal.center.lng)]))
    }
  }

  # generate area from longitude and latitude:
  # if(missing(area) && !missing(longitude) && !missing(latitude)){
  #  n.area <- NCOL(dataset) + 1
  #  dataset[n.area] <- compute.areas(dataset[n.longitude], dataset[n.latitude])
  # }

  # chking for indices inputs:
  if(!exists(n.area)){
    indices[c("duncan.delta", "absolute.concentration", "relative.concentration",
              "absolute.centralization", "relative.centralization",
              "absolute.clustering", "relative.clustering",
              "distance.decay.interaction", "distance.decay.isolation",
              "spatial.proximity"] <- FALSE
    warning("Some of the requested indices cannot be computed because of of the lack of area variable!")
  }

  if(!exists(n.longitude) || !exists(n.latitude)){
    indices[c("absolute.centralization", "relative.centralization",
              "absolute.clustering", "relative.clustering",
              "distance.decay.interaction", "distance.decay.isolation",
              "spatial.proximity"] <- FALSE
    warning("Some of the requested indices cannot be computed because of of the lack of longitude or latitude variable!")
  }
  

  ####
  # MAIN PRG
  ####

  # create the data.frame for the indices:
  # first get group names from data:
  if(is.character(groups)){
    group.names <- groups
  } else {
    if(length(dimnames(dataset[[2]]) && is.numeric(groups)){
      group.names <- dimnames(dataset[[2]])[groups]
    }
  }
  # required indices:
  ind.names <- names(sapply(indices[names(indices) != "atkinson"], "==", TRUE))
  if(length(indices$atkinson$b)){
    ind.names <- c(paste("atkinson", segrInds$atkinson$b, sep="_"), ind.names)
  }
  # generate "groupN_indiceM" variable names
  dnames <- outer(group.names, ind.names, paste, sep="_")

  # get areal units for rows, and for loop
  areal.ids <- unique(dataset[areal.units])

  # create data.frame:
  RESULT <- matrix(NA, ncol=length(dnames.names), nrow=length(areal.ids))
  RESULT <- as.data.frame(RESULT)
  dimnames[[RESULT]] <- list(areal.ids, dnames)

  # compute indices for each area
  for(area in areal.ids){

    # create data storage for the areal id
    areal.indices <- indices
    areal.indices[sapply(areal.indices, "==", TRUE)] <- NA


    areal.data <- dataset[which(dataset$areal.units == area), ]

    # chk for areal.lng & areal.lat consistency
    if(!missing(areal.center.lng) && !missing(areal.center.lat)){
      if(NROW(unique(areal.data[c(n.areal.center.lng, n.areal.center.lat)])) > 1){
        warning(paste("There are more combinations of areal.center.lng and areal.center.lat for the",
          area, "area. Please check them. Indices will use 'GeoCentrum' mode."))
        centrum.mode <- "GeoCentrum"
      }
    }
    # chk for NAs in longitude and latitude
    if(!missing(longitude) && !missing(latitude){
      if(sum(is.na(areal.data[n.longitude])) || sum(is.na(areal.data[n.latitude]))){
        areal.indices[c("absolute.centralization", "relative.centralization",
                        "absolute.clustering", "relative.clustering",
                        "distance.decay.interaction", "distance.decay.isolation",
                        "spatial.proximity")] <- FALSE
        warning(paste("Some of the longitude or latitude variable is missing in areal unit",
                      area ". Please check."))
      }
    }
    # chk NAs for area
    if(exists(n.area)){
      if(sum(is.na(areal.data[n.area]))){
        areal.indices[c("duncan.delta", "absolute.concentration", "relative.concentration",
                        "absolute.centralization", "relative.centralization",
                        "absolute.clustering", "relative.clustering",
                        "distance.decay.interaction", "distance.decay.isolation",
                        "spatial.proximity")]
        warning(paste("Some of the area variable is missinf in areal unit",
                      area, ". Please check."))
      }
    }

    # If necesary, compute distance matrices
    if(areal.indices$distance.decay.isolation || indices$distance.decay.interaction){
      dij <- d_ij(area = areal.data[n.area],
                  lng  = areal.data[n.longitude],
                  lat  = areal.data[n.latitude],
                  mode = distance.mode)
    }
    if(areal.indices$spatial.proximity || indices$absolute.clustering ||
       indices$relative.clustering){
      cij <- c_ij(area = areal.data[n.area],
                  lng  = areal.data[n.longitude],
                  lat  = areal.data[n.latitude],
                  mode = distance.mode)
    }
    if(areal.indices$relative.centralization || areal.indices$absolute.centralization){
      if(centrum.mode = "GeoCentrum"){
        areal.center.lng <- NCOL(areal.data) + 1
        areal.center.lat <- NCOL(areal.data) + 1
        areal.data[c(areal.center.lng, areal.center.lat)] <- 0
      }
      cent_dist <- cent_dist(lng = areal.data[n.longitude],
                             lat = areal.data[n.latitude],
                             mode = centrum.mode,
                             cent_lng = areal.data[areal.center.lng]
                             cent_lat = areal.data[areal.center.lat])
    }

    # compute areal total:
    areal.data$total <- rowSums(areal.data[n.groups])
    TOTAL <- sum(areal.data$total)

    # for each group, compute indices:
    for(group in n.groups){
      # compute basic variables
      X <- sum(areal.data[group])
      y <- areal.data$total - areal.data[group]
      Y <- TOTAL - X
      p <- areal.data[group]/areal.data$total
      P <- X/TOTAL
      if(!exists(n.area)){
        n.area <- NCOL(areal.data) + 1
        areal.data[n.area] <- numeric(NROW(areal.data))
      }
      A <- sum(areal.data[n.area])

      if(!exists(n.longitude)){
        n.longitude <- NCOL(areal.data) + 1
        areal.data[n.longitude] <- numeric(NROW(areal.data))
      }
      if(!exists(n.latitude)){
        n.latitude <- NCOL(areal.data) + 1
        areal.data[n.latitude] <- numeric(NROW(areal.data))
      }


      # compute areal data:
      areal.data <- computeSegrIndices(x=areal.data[group],
                                       X=X,
                                       y=y,
                                       Y=Y,
                                       total=areal.data$total,
                                       TOTAL=TOTAL,
                                       p=p,
                                       P=P,
                                       area=areal.data[n.area],
                                       A = A,
                                       lng= areal.data[n.longitude],
                                       lat= areal.data[n.latitude],
                                       cent_lng = areal.data[areal.center.lng],
                                       cent_lat = areal.data[areal.center.lat],
                                       cent_dist = cent_dist,
                                       cij=cij,
                                       dij=dij)
    }

    # TODO: egy-egy lista minden areal.unit minden csoportjára
    # ezt egy táblázatba kéne szuszakolni, és szépen elnevezni, még nem tudom hogy oldom meg?
  }
}




#-----------------------#
#                       #
#   DEFAULT SETTINGS    #
#                       #
#-----------------------#

segrInds <- list(
  atkinson = list(b=seq(0,1,by=0.1)),
  dissimilarity              = TRUE,
  gini                       = TRUE,
  entropy                    = TRUE,
  interaction                = TRUE,
  isolation                  = TRUE,
  correlation                = TRUE,
  duncan.delta               = TRUE,
  absolute.concentration     = TRUE,
  relative.concentration     = TRUE,
  absolute.centralization    = TRUE,
  relative.centralization    = TRUE,
  absolute.clustering        = TRUE,
  relative.clustering        = TRUE,
  distance.decay.interaction = TRUE,
  distance.decay.isolation   = TRUE,
  spatial.proximity          = TRUE)