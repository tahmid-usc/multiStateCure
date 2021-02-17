library(gganimate)
library(ggplot2)
library(dplyr)
library(simsurv)

#----------------------------------------------------


# Simulate survival data: Weibull baseline hazard with a γ parameter of 1.5
# maximum follow up time by censoring time = 5 years

#-----------------------------------------------------
  
  # Simulate survival data: Weibull baseline hazard with a γ parameter of 1.5
  # maximum follow up time by censoring time = 5 years
  
  
  n <- 10
  cov <- data.frame(id = 1:n,
                    trt = rbinom(n, 1, 0.5))
  # Simulate the event times
  dat <- simsurv(lambdas = 0.1, 
                 gammas = 1.5, 
                 betas = c(trt = -0.5), 
                 x = cov, 
                 maxt = 5)
  # Merge the simulated event times onto covariate data frame
  dat <- merge(cov, dat)
  
  
  # Function for expanding the survival times to for animation
  # datrow = row of a data frame | id | treatment | eventtime | status |
  # numTime = Number of time frame to breakdown single survival time
  # maxt = maximum survival time 
  
  expandEventTime <- function(datrow, numTime = 10, maxt) {
    datrow <- as.numeric(datrow)
    id <- datrow[1]; trt = datrow[2]; eventtime <- datrow[3]; status <- datrow[4]
    
    # replicate id and treatment variable
    id <- rep(id, numTime)
    trt <- rep(ifelse(trt == 1, 'Treatment', 'Control'), numTime)
    exT <- matrix(eventtime, nrow = numTime, ncol = 1)
    exS <- matrix(ifelse(status == 1, 'Dead', 'Censored'), nrow = numTime, ncol = 1)
    eT <- seq(from = 0, to = eventtime, by = (maxt+1)/(numTime)) #by = this can be improved
    exT[1:length(eT)] <- eT
    exS[1:length(eT)] <- 'Alive'
    exS[length(eT)+1] <- ifelse(status == 1, 'Dead', 'Censored')
    
    time = seq(0, maxt, length.out = numTime)
    
    return(data.frame(id = id, Treatment = trt, eventtime = exT, Status = exS, time = time))
  }
  
  # apply the expandEeventTime function to all rows of the dataset
  # user supplies id, treatment variable, event time, status and number of time frames for animation
  
  createData <- function(id, trt, eventtime, status, numTime = 20){
    maxt <- max(eventtime)
    dat <- data.frame(id = id, trt = trt, eventtime = eventtime, status = status)
    dt <- do.call("rbind", apply(dat, 1, expandEventTime, numTime = numTime, maxt = maxt)) 
    return(dt)
  }
  
  # create expanded dataset
  
  dt <- createData(dat$id, dat$trt, dat$eventtime, dat$status)
  
  # crerate animation object
  
  surv_anim <- dt %>%
    ggplot(aes(x = id, y = eventtime, group = id)) +
    geom_line(aes(color = Treatment), size = 2) +
    geom_point(aes(shape = Status), size = 5, color = 3) + 
    scale_x_continuous(breaks = unique(dt$id)) + 
    labs(title = 'Year: {round(frame_along)}', x = 'ID', y = 'Time') +
    theme_gray(base_size =  20) +
    coord_flip() +
    transition_reveal(time) 
  
  #view animation (requires gifski package)
  
  animate(surv_anim, height = 600, width = 800)
  
  # save animation as gif
  anim_save("./survival_animation.gif")
  