#### Get reference data ########################################################

## Boys' reference data from ...
##     https://www.who.int/childgrowth/standards/wfa_boys_0_5_zscores.txt
boys <-  read.table("wfaBoys.csv", header = TRUE, sep = ",")
boys$lowerSD <- (boys$SD0 - boys$SD3neg) / 3
boys$upperSD <- (boys$SD3 - boys$SD0) / 3
boys <- subset(boys, Month >= 6)

## Girls' reference data from ...
##    https://www.who.int/childgrowth/standards/wfa_girls_0_5_zscores.txt
girls <- read.table("wfaGirls.csv", header = TRUE, sep = ",")
girls$lowerSD <- (girls$SD0 - girls$SD3neg) / 3
girls$upperSD <- (girls$SD3 - girls$SD0) / 3
girls <- subset(girls, Month >= 6)

#### Find formulae for boys ####################################################

quartz(width = 9, height = 4.5)
par(mfrow = c(1, 3))

## BOYS : Median weight for age
x <- boys$Month
y <- boys$SD0
plot(x, y, pch = 19, cex = 0.5, xlab = "month", ylab = "median weight", main = "Boys")
fit6 <- lm(y ~ poly(x, 6, raw=TRUE))
xx <- seq(6:72)
lines(xx, predict(fit6, data.frame(x=xx)), col = "red")
as.vector(as.vector(fit6$coefficients))

## BOYS : Lower SD for age
y <- boys$lowerSD
plot(x, y, pch = 19, cex = 0.5, xlab = "month", ylab = "Lower SD for weight", main = "Boys")
fit6 <- lm(y ~ poly(x, 6, raw=TRUE))
xx <- seq(6:72)
lines(xx, predict(fit6, data.frame(x=xx)), col = "red")
as.vector(fit6$coefficients)

## BOYS : Upper SD for age
y <- boys$upperSD
plot(x, y, pch = 19, cex = 0.5, xlab = "month", ylab = "Upper SD for weight", main = "Boys")
fit6 <- lm(y ~ poly(x, 6, raw=TRUE))
xx <- seq(6:72)
lines(xx, predict(fit6, data.frame(x=xx)), col = "red")
as.vector(fit6$coefficients)

#### Find formulae for girls ###################################################

quartz(width = 9, height = 4.5)
par(mfrow = c(1, 3))

## GIRLS : Median weight for age
x <-girls$Month
y <- girls$SD0
plot(x, y, pch = 19, cex = 0.5, xlab = "month", ylab = "median weight", main = "Girls")
fit6 <- lm(y ~ poly(x, 6, raw=TRUE))
xx <- seq(6:72)
lines(xx, predict(fit6, data.frame(x=xx)), col = "red")
as.vector(fit6$coefficients)

## GIRLS : Lower SD for age
y <- girls$lowerSD
plot(x, y, pch = 19, cex = 0.5, xlab = "month", ylab = "Lower SD for weight", main = "Girls")
fit6 <- lm(y ~ poly(x, 6, raw=TRUE))
xx <- seq(6:72)
lines(xx, predict(fit6, data.frame(x=xx)), col = "red")
as.vector(fit6$coefficients)

## GIRLS : Upper SD for age
y <- girls$upperSD
plot(x, y, pch = 19, cex = 0.5, xlab = "month", ylab = "Upper SD for weight", main = "Girls")
fit6 <- lm(y ~ poly(x, 6, raw=TRUE))
xx <- seq(6:72)
lines(xx, predict(fit6, data.frame(x=xx)), col = "red")
as.vector(fit6$coefficients)

### z-score function ###################################################

findWAZ <- function(observedWeight, ageMonths, sex, digits = 2)
{
	if(sex == 1) {
	  referenceMedian <- 4.977587e+00 + 6.850076e-01 * ageMonths + -4.140006e-02 * ageMonths^2 +  1.861115e-03 * ageMonths^3 + -4.557401e-05 * ageMonths^4 + 5.633279e-07 * ageMonths^5 + -2.746621e-09 * ageMonths^6
		lowerSD <-         4.715876e-01 + 6.516336e-02 * ageMonths + -4.267320e-03 * ageMonths^2 +  2.070666e-04 * ageMonths^3 + -5.196451e-06 * ageMonths^4 + 6.479001e-08 * ageMonths^5 + -3.173980e-10 * ageMonths^6
		upperSD <-         6.616327e-01 + 7.466573e-02 * ageMonths + -3.956540e-03 * ageMonths^2 +  1.914075e-04 * ageMonths^3 + -4.684480e-06 * ageMonths^4 + 5.667014e-08 * ageMonths^5 + -2.660587e-10 * ageMonths^6
		}
	if(sex == 2) {
	  referenceMedian <- 4.589059e+00 + 6.135219e-01 * ageMonths + -3.509681e-02 * ageMonths^2 + 1.569346e-03 * ageMonths^3 + -3.784446e-05 * ageMonths^4 + 4.607875e-07 * ageMonths^5 + -2.222577e-09 * ageMonths^6
	  lowerSD <-         4.366394e-01 + 6.709415e-02 * ageMonths + -4.060771e-03 * ageMonths^2 + 1.767080e-04 * ageMonths^3 + -4.090310e-06 * ageMonths^4 + 4.921072e-08 * ageMonths^5 + -2.408582e-10 * ageMonths^6
	  upperSD <-         6.281990e-01 + 1.051817e-01 * ageMonths + -5.260342e-03 * ageMonths^2 + 1.927602e-04 * ageMonths^3 + -3.660493e-06 * ageMonths^4 + 3.792593e-08 * ageMonths^5 + -1.696616e-10 * ageMonths^6
	  }
  if(observedWeight < referenceMedian) {
  	z <- (observedWeight - referenceMedian) / lowerSD
  	} else {
			z <- (observedWeight - referenceMedian) / upperSD
			}
	return(round(z, digits = digits))
	}

### Test against z-scorer library #############################################

# Apply new function
svy <- read.table("ugan03.csv", header = TRUE, sep = ",")
z <- vector(mode = "numeric", length = nrow(svy))
pb <- txtProgressBar(min = 0, max = nrow(svy), style = 1)
for(i in 1:nrow(svy))
	{

  z[i] <- findWAZ(observedWeight = svy$weight[i], ageMonths = svy$age[i], sex = svy$sex[i], digits = 2)
  setTxtProgressBar(pb, i)
	}
cat("\n", sep = "")

## Load test data and generate WFAZ using zscorer library
require(zscorer)
svy <- read.table("ugan03.csv", header = TRUE, sep = ",")
svy$age <- svy$age * (365.25 / 12)
svy <- addWGSR(data = svy, sex = "sex", firstPart = "weight",  secondPart = "age", index = "wfa")

## Differences
diff <- svy$wfaz - z
differences <- data.frame(cbind(svy$wfaz, z, diff)); names(differences) <- c("zscorer", "polyWAZ", "difference")
quartz(width = 9, height = 4.5)
par(mfrow = c(1,2))
hist(differences$difference, xlab = "zscorer - polyWAZ", main = "Test of polyWAZ vs. zscorer")
abline(v = median(differences$difference), col = "red")
plot(differences$zscorer, differences$polyWAZ, cex = 0.5, xlab = "WAZ from zscorer", ylab = "WAZ from polyWAZ", main = "Test of polyWAZ vs. zscorer")
abline(a = 0, b = 1, col = "red")
summary(differences)
