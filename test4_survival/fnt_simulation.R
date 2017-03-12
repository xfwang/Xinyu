## function: simulate survival data with Weibull parameters
sim_data <- function(alpha, mu, Nobs, Ncen) {
  observed_data <- data.frame(os_status = rep_len('DECEASED', Nobs),
                              os_months = rweibull(n = Nobs, alpha, exp(-(mu)/alpha)),
                              stringsAsFactors = F
  )
  
  censored_data <- data.frame(os_status = rep_len('LIVING', Ncen),
                              os_months = runif(Ncen) * rweibull(Ncen, alpha, exp(-(mu)/alpha)),
                              stringsAsFactors = F
  )
  
  return(observed_data %>% bind_rows(censored_data))
}

## function : simulate proportional harzard survival data
##            alpha     - shape parameter
##            mu        - a sequence of scale parameters
##            idx.event - logical variable, event or censor, same length as mu
sim_ph <- function(alpha, mu, idx.event) {
  if(length(mu) != length(idx.event)) error("mu and idx.event should have matched dimension! ")
  nobs <- sum(idx.event)
  ncen <- sum(!idx.event)
  observed_data <- data.frame(os_status = rep_len('DECEASED', nobs),
                              os_months = rweibull(nobs, alpha, exp(-(mu[idx.event])/alpha)),
                              stringsAsFactors = F
  )
  
  censored_data <- data.frame(os_status = rep_len('LIVING', ncen),
                              os_months = runif(ncen) * rweibull(ncen, alpha, exp(-(mu[!idx.event])/alpha)),
                              stringsAsFactors = F
  )
  
  return(observed_data %>% bind_rows(censored_data))
}


gen_inits <- function() {
  list(
    alpha_raw = 0.01*rnorm(1),
    mu = rnorm(1)
  )
}
