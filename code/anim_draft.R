library(gganimate)
library(ggplot2)
library(dplyr)
library(simsurv)


# Merge the simulated event times onto covariate data frame
dat <- merge(cov, dat) %>% 
  mutate(event = ifelse(status == 1, 'Dead', 'Censored'), 
         time = seq(0, max(eventtime), length.out = nrow(dat)))

p <- dat %>% 
  ggplot(aes(x = id, y = eventtime)) +
  geom_linerange(aes(ymin = 0, ymax = eventtime)) +
  geom_point(aes(shape = event, color = event), stroke = 1, cex = 2) +
  scale_shape_manual(values = c(1, 3, 4)) +
  labs(y = "Time (years)", x = "Subject ID") + coord_flip() + theme_classic()


p +   transition_states(Species,
                        transition_length = 2,
                        state_length = 1)


dt <- data.frame(time = 1:100, x = 1:100, y = (1:100)^2) %>% mutate(y = ifelse(x > 50, NA, y))

anim <- dt %>% ggplot(aes(x = x, y = y)) +
  geom_line() +
  transition_reveal(time) 

anim


dt <- data.frame(time = 1:100, x = 1, y = 100) 

anim <- dt %>% ggplot(aes(x = x, y = y)) +
  geom_linerange(aes(ymin = 0, ymax = y)) +
  transition_reveal(time) 

anim

#--------------------


n <- 10
cov <- data.frame(id = 1:n,
                  trt = rbinom(n, 1, 0.5))
# Simulate the event times
dat <- simsurv(lambdas = 0.1, 
               gammas = 1.5, 
               betas = c(trt = -0.5), 
               x = cov, 
               maxt = 5)


# total time frame

id <- 1
time <- seq(0, max(dat$eventtime), length.out = 100)
y <- seq(0, dat$eventtime[1], length.out = 100) 

dt <- data.frame(id, time, y)

anim <- dt %>% ggplot(aes(x = id, y = y)) +
  geom_linerange(aes(ymin = 0, ymax = y)) +
  transition_reveal(time) +
  coord_flip()

anim




#--------------------



n <- 10
cov <- data.frame(id = 1:n,
                  trt = rbinom(n, 1, 0.5))
# Simulate the event times
dat <- simsurv(lambdas = 0.1, 
               gammas = 1.5, 
               betas = c(trt = -0.5), 
               x = cov, 
               maxt = 5)


# function for expanding event time

maxt <- max(dat$eventtime)

expandEventTime <- function(datrow, numTime = 10) {
  datrow <- as.numeric(datrow)
  id <- datrow[1]; eventtime <- datrow[2]; status <- datrow[3]
  id <- rep(id, numTime)
  exT <- matrix(eventtime, nrow = numTime, ncol = 1)
  exS <- matrix(ifelse(status == 1, 'Dead', 'Censored'), nrow = numTime, ncol = 1)
  eT <- seq(from = 0, to = eventtime, by = maxt/(numTime-1))
  exT[1:length(eT)] <- eT
  exS[1:length(eT)] <- 'Alive'
  exS[length(eT)+1] <- ifelse(status == 1, 'Dead', 'Censored')
  
  time = seq(0, maxt, length.out = numTime)
  
  return(data.frame(id = id, eventtime = exT, status = exS, time = time))
}

dt <- do.call("rbind", apply(dat, 1,expandEventTime)) 


dt %>% 
  ggplot(aes(x = id, y = eventtime, group = id)) +
  geom_line() +
  geom_point(aes(shape = status, color = status), stroke = 1, cex = 2) + 
  coord_flip() +
  transition_reveal(time)


#---------------------------

expandEventTime <- function(datrow, numTime = 10) {
  datrow <- as.numeric(datrow)
  id <- datrow[1]; eventtime <- datrow[2]; status <- datrow[3]
  id <- rep(id, numTime)
  exT <- matrix(eventtime, nrow = numTime, ncol = 1)
  exS <- matrix(ifelse(status == 1, 'Dead', 'Censored'), nrow = numTime, ncol = 1)
  eT <- seq(from = 0, to = eventtime, by = (maxt+1)/(numTime))
  exT[1:length(eT)] <- eT
  exS[1:length(eT)] <- 'Alive'
  exS[length(eT)+1] <- ifelse(status == 1, 'Dead', 'Censored')
  
  time = seq(0, maxt, length.out = numTime)
  
  return(data.frame(id = id, eventtime = exT, status = exS, time = time))
}

createData <- function(id, eventtime, status, numTime = 20){
  maxt <- max(eventtime)
  dat <- data.frame(id = id, eventtime = eventtime, status = status)
  dt <- do.call("rbind", apply(dat, 1, expandEventTime, numTime = numTime)) 
  return(dt)
}

surv_anim <- createData(dat$id, dat$eventtime, dat$status) %>%
  ggplot(aes(x = id, y = eventtime, group = id)) +
  geom_line(size = 2) +
  geom_point(aes(shape = status, color = status), size = 5, cex = 2) + 
  labs(title = 'Year: {round(frame_along)}', x = 'Time', y = 'ID') +
  theme_gray(base_size =  20) +
  coord_flip() +
  transition_reveal(time) 

animate(surv_anim, height = 600, width = 800)
anim_save("./survival_animation.gif")
