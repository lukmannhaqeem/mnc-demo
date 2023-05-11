# Orangutan population estimate using MNC

set.seed(1234)

# Load packages
library(jagsUI)
library(wiqid)
library(raster)

#read the data
ou_data <- read.csv("ou_sdm_analysis.csv", header = TRUE)

#'see' the first six rows of the data frame
head(ou_data)

#plot the data frame objects
plot(ou_data$x, ou_data$y, col = "red")
#plot(ou_data$x, ou_data$y, col = "red", cex = log(ou_data$team1+ou_data$team2))
text(ou_data$x, ou_data$y-2000, labels = ou_data$id, cex = 0.6)

par(mfrow = c(1, 2))	#set the graphical parameters (1 row, 2 columns)

#computes a histogram of the given data frame objects
hist(ou_data$elevation)	#do this for the ecological covariates
hist(ou_data$forest)

#produce the results of various arithmetics of a function
summary(ou_data$elevation)	
summary(ou_data$forest)

hist(ou_data$team1)		#do this for the observational covariates
hist(ou_data$team2)

#returns the sum of all the values present
sum(ou_data$team1)
sum(ou_data$team2)

par(mfrow = c(1, 1))	#reset the graphical parameters to original

# Simple point estimate
( tmp <- colSums(ou_data) )
( nn <- tmp['newnests'] / tmp['period'] / 1.68*100 ) # nests per sq km per day
nn / 1 / 0.85 # ou per sq km


#construct model in BUGS language
#the model vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
model.txt <- "
model {

  # Estimation of probability of detection
	team1 ~ dbin(q, x) 		      # seen by team a
	team2 ~ dbin((1-q) * q, x)	  # missed by team a, seen by team b

	for(i in 1:N) {
    # realised total number of new nests per unit effort
	n[i] ~ dpois(lambda[i] * area[i] * period[i])
	
	# expected total number of new nests
    log(lambda[i]) <- a0 +
					  a1*elevation[i] +
					  a2*forest[i]
	
	# observed total number of new nests
	y[i] ~ dbin(q, n[i])
	
	# ou density
    d[i] <- n[i] / area[i] / period[i]
	}
	
  # Priors:
	x0 ~ dunif(0, 1000)
	x <- trunc(x0)
	q ~ dunif(0, 1)
	a0 ~ dunif(-10, 10)
	a1 ~ dunif(-10, 10)
	a2 ~ dunif(-10, 10)
} "
writeLines(model.txt, "model_2.txt")
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#'scale' or 'centre' the data to:
# 1. ensure that all data are on the same scale
# 2. help models to converge faster
ou_data$elevationS <- (ou_data$elevation - mean(ou_data$elevation)) / sd(ou_data$elevation)
ou_data$forestS <- (ou_data$forest - mean(ou_data$forest)) / sd(ou_data$forest)

#input data for JAGS:
JAGSdata <- list(
				area = ou_data$area,
				y = ou_data$newnests,
				elevation = ou_data$elevationS,
				forest = ou_data$forestS,
				period = ou_data$period,
				N = nrow(ou_data),
				team1 = sum(ou_data$team1),
				team2 = sum(ou_data$team2))
str(JAGSdata)

#set parameters to be monitored 
wanted 	<- c("lambda", "q", "x", "a0", "a1", "a2")

#set the starting values of some variables
init <- function() list(x0 = runif(1, ou_data$team1+ou_data$team2, 1000),
						n = rep(1000, nrow(ou_data)))

#run JAGS from R
out2 <- jags(JAGSdata,
			init,
			wanted,
			"model_2.txt",
			n.chains = 3,
			n.iter = 101000,
			n.burnin = 1000,
			n.thin = 1)
out2

#check r-hat
mean(unlist(out2$Rhat) > 1.1)	# > 0 means there are lack of convergence

#check the distribution of the coefficients
Wout2 <- as.Bwiqid(out2)
plot(Wout2)  			# Nest construction rate, nests / day / ha
plot(Wout2, "lambda")  	# Nest construction rate, nests / day / ha
plot(Wout2, "q")  		# detection probability
plot(Wout2, "a0")  		# intercept coefficient
plot(Wout2, "a1")  		# elevation effect coefficient
plot(Wout2, "a2")  		# forest effect coefficient

( sampSize <- nrow(Wout2) )

#calculate proportion of nest-builders (p) based on output
p_depend <- rbeta(sampSize, 14+1, (78-14)+1)
p <- 1 / (1 + p_depend)	# Proportion of adults in the population = 1 / (1 + p_depend)
plotPost(p)

#calculate nest construction rate (r) based on output
r <- rgamma(sampSize, 602+1, 602+1)
plotPost(r)


#Interpolate orangutan population size along:
#elevation gradient
min.elevationS <- min(ou_data$elevationS)
max.elevationS <- max(ou_data$elevationS)
data_for_prediction <- data.frame(
							elevationS = seq(min.elevationS, max.elevationS, length = 100), 
							forestS = rep(mean(ou_data$forestS), length = 100))
#calculate population size of the given elevation (assuming forest cover effect is constant) 
population_prediction <- exp(
					 out2$mean$a0 +
					 out2$mean$a1*data_for_prediction$elevationS +
					 out2$mean$a2*data_for_prediction$forestS) / mean(p) / mean(r) * 100
mu <- mean(ou_data$elevation)
stdev <- sd(ou_data$elevation)
unscaled.elevation <- (data_for_prediction$elevation * stdev) + mu 
par(mfrow = c(1, 2))
plot(unscaled.elevation, population_prediction, type = "l",
	xlab = "elevation (m a.s.l.)", ylab = "population size")

#forest gradient
min.forestS <- min(ou_data$forestS)
max.forestS <- max(ou_data$forestS)
data_for_prediction <- data.frame(
							elevation = rep(mean(ou_data$elevationS), length = 100), 
							forest = seq(min.forestS, max.forestS, length = 100))
#calculate population size of the given forest cover (assuming elevation effect is constant) 
population_prediction <- exp(
					 out2$mean$a0 +
					 out2$mean$a1*data_for_prediction$elevation +
					 out2$mean$a2*data_for_prediction$forest) / mean(p) / mean(r) * 100
mu <- mean(ou_data$forest)
stdev <- sd(ou_data$forest)
unscaled.forest <- (data_for_prediction$forest * stdev) + mu 
plot(unscaled.forest, population_prediction, type = "l",
	xlab = "forest cover (%)", ylab = "population size")
par(mfrow = c(1, 1))


#interpolate orangutan population across study area
#using species distribution model (sdm)
study_area <- read.csv("study_area.csv", header = TRUE)
head(study_area)

study_area$elevationS <- (study_area$elevation - mean(ou_data$elevation)) / sd(ou_data$elevation)
study_area$forestS <- (study_area$forest - mean(ou_data$forest)) / sd(ou_data$forest)

#calculate population size based on all covariate effects 
study_area$population_prediction <- exp(
		 out2$mean$a0 +
		 out2$mean$a1*study_area$elevationS +
		 out2$mean$a2*study_area$forestS) / mean(p) / mean(r) * 100

#plot elevation map
elevation_map <- rasterFromXYZ(data.frame(x = study_area$x, 
								   y = study_area$y,
								   z = study_area$elevation))
par(mfrow = c(1, 3))
plot(elevation_map, main = "Elevation (m a.s.l.)")

#plot forest cover map
forest_map <- rasterFromXYZ(data.frame(x = study_area$x, 
								   y = study_area$y,
								   z = study_area$forest))
plot(forest_map, main = "Forest cover (%)")

#plot population distribution map
population_map <- rasterFromXYZ(data.frame(x = study_area$x, 
								   y = study_area$y,
								   z = study_area$population_prediction))
plot(population_map, main = "Orang-utan distribution")
par(mfrow = c(1, 1))
