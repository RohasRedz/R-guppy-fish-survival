## Cleaning the memory
rm(list=ls(all=TRUE))

## Install packages

## install.packages("package name")
# install.packages(c("devtools","mvtnorm","loo","coda"), repos="https://cloud.r-project.org/",dependencies=TRUE)
# library(devtools)

#install.packages("package name","package name")
# Load packages 
library("lme4")
library(rethinking)
library(readxl)


# get working directory
getwd()

# if necessary you can change the working directory
setwd("data/")


# Read data
Survival <- read_excel("Reznick2002/Survival.xls")


length(complete.cases(Survival))
# inspect the strucure of your data
str(Survival)


#enconding categorical data
Survival$Stream = factor(Survival$Stream)
levels(Survival$Stream) # Three levels of density manipulation

Survival$replicate = factor(Survival$replicate)
levels(Survival$replicate) # Three levels of density manipulation


# inspect a variable

Survival$TREATMENT # print the variable

levels(Survival$TREATMENT) # why null?

Survival$TREATMENT = factor(Survival$TREATMENT)

levels(Survival$TREATMENT) # Three levels of density manipulation
levels(Survival$TREATMENT) = c(1,0.5,2) # Change the levels values for other values

Survival$TREATMENT


# Frequentist approach

# Generalised linear model (binomial family) 
f.mod.1 <-  glm(survival ~ TREATMENT, data= Survival, family = binomial)
summary(f.mod.1)

# install.packages("car")
library(car)
Anova(f.mod.1)
Anova(f.mod.1, type = 'III')



# Generalised linear mixed model (binomial family) 
f.mod.2 <-  glmer(survival ~ TREATMENT + (1|Stream), data= Survival, family = binomial)
summary(f.mod.2)
Anova(f.mod.2)
Anova(f.mod.2, type = 'III')



# Generalised linear mixed model (binomial family) 
f.mod.3 <-  glmer(survival ~ TREATMENT + (1|Stream) + (1|year), data= Survival, family = binomial)
summary(f.mod.3)
Anova(f.mod.3)
Anova(f.mod.3, type = 'III')


# which is the best model?
AIC(f.mod.1, f.mod.2, f.mod.3)

#f.mod.2 is better model with lowest AIC

##########################################################################################################################
## Now,  in a Bayesian model using stan through the rethinking package


### Generalised linear model



glimmer(survival ~ TREATMENT, data= Survival, family = binomial)
str(Survival)

# make sure that the variables have the right format
Survival$survival <- as.numeric(Survival$survival)
Survival$Density <- as.numeric(as.character(Survival$TREATMENT))


# Make the variables for the model
x_data <- as.data.frame(model.matrix(~TREATMENT, data= Survival))


# Make it as a list
d <- list(
  
  survival <- Survival$survival,
  TREATMENT0_5 <- as.numeric(x_data$TREATMENT0.5),
  TREATMENT2 <- as.numeric(x_data$TREATMENT2)
)



# Run the model using a maximun a posteriory (normaly for simpler models)

b.mod.1 <- map(
  
  alist(
    survival ~ dbinom( 1 , p ),
    logit(p) <- Intercept +
      b_TREATMENT0_5*TREATMENT0_5 +
      b_TREATMENT2*TREATMENT2,
    Intercept ~ dnorm(0,10),
    b_TREATMENT0_5 ~ dnorm(0,10),
    b_TREATMENT2 ~ dnorm(0,10)
  ), data = d
)


# summary statistics Bayesian models
precis(b.mod.1, prob = 0.95, digits = 3)

# Interpredations
logistic(0.92)


# summary frequentist model f.mod.2
summary(f.mod.2)
logistic(0.9)



## Get the model in the bayesian model for the rethinking package

glimmer(survival ~ TREATMENT +  (1|Stream) + (1|replicate), data= Survival, family = binomial)
str(Survival)


# Make the data for the model
x_data <- as.data.frame(model.matrix(~TREATMENT, data= Survival))
levels(Survival$Stream) <- as.factor(c(1:5))


# make a data frame
d <- as.data.frame(cbind(
  survival <- Survival$survival,
  TREATMENT0_5 <- as.numeric(x_data$TREATMENT0.5),
  TREATMENT2 <- as.numeric(x_data$TREATMENT2), 
  Stream <- Survival$Stream,
  replicate <- Survival$replicate
))

# set the names of the variables
names(d) <- c("survival", "TREATMENT0_5", "TREATMENT2", "Stream", "replicate")


# Run the model
b.mod.2 <- map2stan(
  
  alist(
    survival ~ dbinom( 1 , p ),
    logit(p) <- Intercept +
      b_TREATMENT0_5*TREATMENT0_5 +
      b_TREATMENT2*TREATMENT2 +
      v_Stream_Intercept[Stream],
    Intercept ~ dnorm(0,10),
    b_TREATMENT0_5 ~ dnorm(0,10),
    b_TREATMENT2 ~ dnorm(0,10),
    v_Stream_Intercept[Stream] ~ dnorm(0,sigma_Stream),
    sigma_Stream ~ dcauchy(0,2)
  ), data = d, chains = 3, cores=3, iter = 3000
  
)




stancode(b.mod.2) # if you want to see the model in Stan code.

# Evaluate if the chains converge
tracerplot(b.mod.2)
graphics.off()

precis(b.mod.2, prob = 0.95)
summary(f.mod.2)

post <- extract.samples(b.mod.2)
str(post)



# plot the posterior probability of survival



# Create a function to extract the predictive values from the posterior
p.link <- function( T_0.5, T_2){
  
  Y.est <- with(post, Intercept + b_TREATMENT0_5 * T_0.5 + b_TREATMENT2 * T_2 )
  return(1/ (1+ exp(-Y.est)))
  
}


# Plot the posterior distribution of the survival probability
dens(p.link(T_0.5 = 1, T_2 = 0), show.HPDI = .95) # Density 0.5
dens(p.link(T_0.5 = 0, T_2 = 0), show.HPDI = .95, add = T, col = "red") # density 1
dens(p.link(T_0.5 = 0, T_2 = 1), show.HPDI = .95, add = T, col = "blue") # density 2


# testing hypothesis
# H.1 = increasing guppy density reduces the probability of survival:

H.1 <- p.link(T_0.5 = 0, T_2 = 1) - p.link(T_0.5 = 0, T_2 = 0)

dens(H.1, show.HPDI = .95, show.zero = T) # Density 0.5
mean(H.1) * 100
abline(v = mean(H.1))
chainmode(H.1) * 100
abline(v = chainmode(H.1), col= "red")


# What is the probability that increasing guppy density by 2x reduces survival?
pp <- length(which(H.1 < 0)) / length(H.1) * 100 
paste(round(pp,3), "%")




# Exercise: Does reducing guppy density to 0.5x increases the probability of survival? what is the probability?

H.1 <- p.link(T_0.5 = 1, T_2 = 0) - p.link(T_0.5 = 0, T_2 = 0)

dens(H.1, show.HPDI = .95, show.zero = T) # Density 0.5
mean(H.1) * 100
abline(v = mean(H.1))
chainmode(H.1) * 100
abline(v = chainmode(H.1), col= "red")


# What is the probability that reducing guppy density to 0.5x survival?
pp <- length(which(H.1 < 0)) / length(H.1) * 100 
paste(round(pp,3), "%")




## Another way of plotting

#make data frame to plot the data
plot.data <- matrix(nrow =3, ncol = 4)
plot.data <- as.data.frame(plot.data)
names(plot.data)  
names(plot.data) <- c("Survival", "LCI", "UCI",  "Density")

# Fill the data frame

plot.data[1,1:3] <- c( chainmode(p.link(T_0.5 = 1, T_2 = 0)), 
                       HPDI(p.link(T_0.5 = 1, T_2 = 0), prob = .95)[1], HPDI(p.link(T_0.5 = 1, T_2 = 0), prob = .95)[2])

plot.data$Density[1] <- "0.5"


# Fill up the other cells

plot.data[2,1:3] <- c( chainmode(p.link(T_0.5 = 0, T_2 = 0)), 
                       HPDI(p.link(T_0.5 = 0, T_2 = 0), prob = .95)[1], HPDI(p.link(T_0.5 = 0, T_2 = 0), prob = .95)[2])

plot.data$Density[2] <- "1"




plot.data[3,1:3] <- c( chainmode(p.link(T_0.5 = 0, T_2 = 1)), 
                       HPDI(p.link(T_0.5 = 0, T_2 = 1), prob = .95)[1], HPDI(p.link(T_0.5 = 0, T_2 = 1), prob = .95)[2])

plot.data$Density[3] <- "2"


plot.data

# Make a plot with error bars

plot <- ggplot(plot.data, aes(x=Density, y=Survival)) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1) + geom_point(cex = 3)

plot


# Make it nicer

plot + xlab("Density treatment") + ylab("% Survival")

graphics.off()


