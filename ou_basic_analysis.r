# Orangutan population estimate using MNC

set.seed(1234)

# Load packages
library(jagsUI)
library(wiqid)

#read the data
ou_data <- read.csv("ou_basic_analysis.csv", header = TRUE)

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
    n[i] ~ dpois(lambda * area[i] * period[i])
	d[i] <- n[i] / area[i] / period[i]
	y[i] ~ dbin(q, n[i])

	}
  # Priors:
	x0 ~ dunif(0, 1000)
	x <- trunc(x0)
	lambda ~ dunif(0, 3)
	q ~ dunif(0, 1)
} "
writeLines(model.txt, "model_1.txt")
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#input data for JAGS:
JAGSdata <- list(
				area = ou_data$area,
				y = ou_data$newnests,
				period = ou_data$period,
				N = nrow(ou_data),
				team1 = sum(ou_data$team1),
				team2 = sum(ou_data$team2))
str(JAGSdata)

#set parameters to be monitored 
wanted 	<- c("lambda", "q", "x")

#set the starting values of some variables
init <- function() list(x0 = runif(1, ou_data$team1+ou_data$team2, 1000))

#run JAGS from R
out1 <- jags(JAGSdata,
			 init,
			 wanted,
			 "model_1.txt",
			 n.chains = 3,
			 n.iter = 101000,
			 n.burnin = 1000,
			 n.thin = 1)
out1

#check r-hat
mean(unlist(out1$Rhat) > 1.1)	# > 0 means there are lack of convergence


Wout <- as.Bwiqid(out1)
plot(Wout)  			# Nest construction rate, nests / day / ha
plot(Wout, "lambda")  	# Nest construction rate, nests / day / ha
plot(Wout, "q")  		# detection probability

( sampSize <- nrow(Wout) )

# Other parameters from 
# Ancrenaz, M., Calaque, R., & Lackman-Ancrenaz, I. (2004) Orangutan nesting behavior in disturbed forest of Sabah, Malaysia: implications for nest census. International Journal of Primatology, 25, 983-1000.

# Estimate of proportion of nest builders
# Assuming the individuals observed are a representative sample, and that observations
# of groups (not individuals) are independent:
# 78 adults encountered, of which 14 had dependent young.
# Binomial posterior for prob that an adult has dependent with uniform prior will be
#   Beta(14 + 1, (78-14) + 1)
# Draw 3e5 values from this posterior:
p_depend <- rbeta(sampSize, 14+1, (78-14)+1)
plotPost(p_depend)
# Proportion of adults in the population = 1 / (1 + p_depend)
p <- 1 / (1 + p_depend)
plotPost(p)

# Estimate of nests per day per nest-builder
# During 602 days of follows, construction of 602 new nests was recorded.
# (Many old nests were reused.)
# Assume Poisson distribution, posterior with Gamma(1, 1) prior
r <- rgamma(sampSize, 602+1, 602+1)
plotPost(r)

# Posterior distribution for orangutan density (animals / sq km)
names(Wout)
# lambda1 is high density, lambda2 is medium density
ou_density <- Wout$lambda / p / r * 100 # lambda is per ha (10000m2), convert to per sq km.
plotPost(ou_density, xlim = c(0, 150), xlab = "Density in high-density area")

ou_pop <- ou_density * 85.85   
mean(ou_pop)
hdi(ou_pop)
plotPost(ou_pop, xlab = "Total number of orang-utan")



