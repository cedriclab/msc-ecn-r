#
#	Université de Sherbrooke
#	ECN802, Hiver 2015
#	Prof: Patrick Richard
#	Devoir 1
#	Cédric Levasseur-Laberge
#

#	IMPORTANT: Indiquer le répertoire de travail avec setwd("path/to/directory") avant de commencer

#	setwd("../../ECN802/Devoir 2");
library("foreign");
library("sandwich");
library("locpol");
library("KernSmooth");
library("gmm");


#	UTILITAIRES

#	Définit une fonction locale pour les retards
lagBy <- function(x,k=1){
	lagNAs <- rep(NA, k);
	laggedVector <- x[1:(length(x)-k)]
	
	c(lagNAs, laggedVector)
}


#	EXERCICES

#################
#	QUESTION 1	#
#################

question1 <- read.dta("question1.dta");
q1 <- data.frame(
	cbind(
		"obs"=c(question1["obs"][[1]]),
		"dy"=c(question1["dy"][[1]]),
		"dc"=c(question1["dc"][[1]])
	)
);
q1n <- length(q1$obs);

q1["dyl1"] <- lagBy(q1$dy, 1);
q1MCO <- lm(dc ~ dy + dyl1, data=q1);
summary(q1MCO);

q1["uchap"] <- c(NA, q1MCO$residuals);
q1ar <- ar(q1$uchap[2:238]);
print(q1ar);

q1["uchapl1"] <- lagBy(q1$uchap);
q1["uchapl2"] <- lagBy(q1$uchap,2);
q1["uchapl3"] <- lagBy(q1$uchap,3);
q1["uchapl4"] <- lagBy(q1$uchap,4);
q1["uchapl5"] <- lagBy(q1$uchap,5);

q1ar5 <- lm(uchap ~ uchapl1 + uchapl2 + uchapl3 + uchapl4 + uchapl5, data=q1);
summary(q1ar5);

q1nw1 <- NeweyWest(q1MCO, lag=1);
print(q1nw1);

q1nw2 <- NeweyWest(q1MCO, lag=2);
print(q1nw2);

q1nw5 <- NeweyWest(q1MCO, lag=5);
print(q1nw5);

q1nw10 <- NeweyWest(q1MCO, lag=10);
print(q1nw10);

#	Brute-force le p optimal
q1diffs <- rep(0, 100);
for (i in 1:100) {
	localNW <- NeweyWest(q1MCO, lag=i);
	q1diffs[i] <- abs(localNW[5] - localNW[9]);
}
q1p <- which.min(q1diffs);
print(q1p);
q1opt <- NeweyWest(q1MCO, lag=q1p);
print(q1opt);

q1sigma <- mean(c(sqrt(q1opt[5]), sqrt(q1opt[9])));
print(q1sigma);


#################
#	QUESTION 2	#
#################

question2 <- read.dta("question2v12.dta");
q2 <- data.frame(
	cbind(
		"temps"=c(question2["temps"][[1]]),
		"y"=c(question2["y"][[1]]),
		"c"=c(question2["c"][[1]]),
		"dy"=c(question2["dy"][[1]]),
		"dy1"=c(question2["dy1"][[1]]),
		"dc1"=c(question2["dc1"][[1]]),
		"dy12"=c(question2["dy12"][[1]])
	)
);
q2n <- length(q2$temps);

#	Génère dc
q2["dc"] <- rep(NA, q2n);
for (i in 2:q2n) {
	q2$dc[i] <- q2$c[i] - q2$c[i-1];
}

q2mco <- lm(dc ~ dy + dy1, data=q2);
summary(q2mco);
q2$uchapMCO <- c(NA, NA, q2mco$residuals);
q2validityTest <- lm(uchapMCO ~ dc1 + dy1, data=q2);
summary(q2validityTest);

q2inst <-lm(dy ~ dc1 + dy1, data=q2);
summary(q2inst);

q2["dychap"] <- c(NA, NA, q2inst["fitted.values"][[1]]);

q2iv2 <- lm(dc ~ dychap + dy1, data=q2);
summary(q2iv2);

q2instrumentMatrix <- as.matrix(
	cbind(
		"dc1" = c(q2$dc1),
		"dy1" = c(q2$dy1)
	)
);

q2gmm1 <- gmm(dc ~ dy + dy1, q2instrumentMatrix, type="twoStep", vcov="HAC", prewhite=1, data=q2);
summary(q2gmm1);

q2gmm2 <- gmm(dc ~ dy + dy1, q2instrumentMatrix, type="twoStep", vcov="HAC", prewhite=2, data=q2);
summary(q2gmm2);

q2gmm3 <- gmm(dc ~ dy + dy1, q2instrumentMatrix, type="twoStep", vcov="HAC", prewhite=3, data=q2);
summary(q2gmm3);

q1["dy1"] <- lagBy(q1$dy, 1);
q1["dy2"] <- lagBy(q1$dy, 2);
q1["dc1"] <- lagBy(q1$dc, 1);

q12mco <-lm(dc ~ dy + dy1, data=q1);
q1["instUchap"] <- c(NA, q12mco$residuals);

q12validity <- lm(instUchap ~ dc1 + dy1 + dy2, data=q1);
summary(q12validity);

q12inst <-lm(dy ~ dc1 + dy1 + dy2, data=q1);
summary(q12inst);

q1["dychap"] <- c(NA, NA, q12inst["fitted.values"][[1]]);

q12iv2 <- lm(dc ~ dychap + dy1, data=q1);
summary(q12iv2);

q12instrumentMatrix <- as.matrix(
	cbind(
		"dc1" = c(q1$dc1),
		"dy1" = c(q1$dy1),
		"dy2" = c(q1$dy2)
	)
);

q12gmm1 <- gmm(dc ~ dy + dy1, q12instrumentMatrix, type="twoStep", vcov="HAC", prewhite=1, data=q1);
summary(q12gmm1);

q12gmm2 <- gmm(dc ~ dy + dy1, q12instrumentMatrix, type="twoStep", vcov="HAC", prewhite=2, data=q1);
summary(q12gmm2);

q12gmm3 <- gmm(dc ~ dy + dy1, q12instrumentMatrix, type="twoStep", vcov="HAC", prewhite=3, data=q1);
summary(q12gmm3);

q12gmmIter <- gmm(dc ~ dy + dy1, q12instrumentMatrix, type="iterative", vcov="HAC", prewhite=6, data=q1);
summary(q12gmmIter);

q1["uchapGMM"] <- c(NA, NA, q12gmmIter$residuals);
q12sargan <- lm(uchapGMM ~ dc1 + dy1 + dy2, data=q1);
summary(q12sargan);

#	LIML, fait sur Stata: ivregress liml dc L1.dy (dy = L1.dy L2.dy L1.dc)


#################
#	QUESTION 3	#
#################

epanechnikovKernel <- function(z){
	if (abs(z) >= 2.236068) {
		return(0);
	} else {
		return((3*(1-((z*z)/5)))/8.944272);
	}
}

localEstimator <- function(xVect, yVect, x0, h){
	le <- 0;
	n <- length(xVect);
	for (i in 1:n){
		le <- le + (yVect[i]*epanechnikovKernel((x0-xVect[i])/h));
	}
	return(le/n);
}

question3 <- read.dta("question3.dta");
q3 <- data.frame(
	cbind(
		"obs"=c(question3["obs"][[1]]),
		"x"=c(question3["x"][[1]]),
		"y"=c(question3["y"][[1]])
	)
);
q3n <- length(q3$obs);

q3kdensity1 <- density(q3$x, bw=0.001, kernel="epanechnikov");
print(q3kdensity1);
plot(q3kdensity1$x, q3kdensity1$y);

q3kdensity2 <- density(q3$x, bw=0.1, kernel="epanechnikov");
print(q3kdensity2);
plot(q3kdensity2$x, q3kdensity2$y);

q3kdensity3 <- density(q3$x, bw=0.2, kernel="epanechnikov");
print(q3kdensity3);
plot(q3kdensity3$x, q3kdensity3$y);

q3kdensity4 <- density(q3$x, bw=1, kernel="epanechnikov");
print(q3kdensity4);
plot(q3kdensity4$x, q3kdensity4$y);

q3kdensity5 <- density(q3$x, bw=10, kernel="epanechnikov");
print(q3kdensity5);
plot(q3kdensity5$x, q3kdensity5$y);

hopt <- 1.059*sqrt(var(q3$x))*(q3n^(-0.2));
q3kdensityOpt <- density(q3$x, bw=hopt, kernel="epanechnikov");
print(q3kdensityOpt);
plot(q3kdensityOpt$x, q3kdensityOpt$y);

print(localEstimator(q3$x, q3$y, hopt));

lcePoly <- locpoly(q3$x, q3$y, kernel="epanechnikov", bandwidth=hopt, drv=0);
print(localEstimator(lcePoly$x, lcePoly$y, hopt));

llePoly <- locpoly(q3$x, q3$y, kernel="epanechnikov", bandwidth=hopt, drv=1);
print(localEstimator(llePoly$x, llePoly$y, hopt));
