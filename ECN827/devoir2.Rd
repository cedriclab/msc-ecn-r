#
#	Université de Sherbrooke
#	ECN827, Hiver 2015
#	Prof: Patrick Richard
#	Devoir 2
#	Cédric Levasseur-Laberge
#

#	IMPORTANT: Indiquer le répertoire de travail avec setwd("path/to/directory") avant de commencer

#	setwd("../../ECN827/Devoir 2")
library("foreign");
library("fGarch");

#	UTILITAIRES

#	Définit le 4eme moment
kurtosis <- function(varVector){
	n <- length(varVector);
	if (n != 0) {
		return(sum(varVector^4)/n);
	} else {
		return(0);
	}
};

#	Définit une fonction locale de lag, celle native de R faisant parfois défaut
lagBy <- function(x,k=1){
	lagNAs <- rep(NA, k);
	laggedVector <- x[1:(length(x)-k)]
	
	c(lagNAs, laggedVector)
}

#	Définit le AIC
aic <- function(sigma2, n, m){
	n*log(sigma2, base=exp(1)) + 2*m
}

#	Définit le BIC
bic <- function(sigma2, n, m){
	n*log(sigma2, base=exp(1)) + m*log(n, base=exp(1))
}

#	Définit le Test de Andrews
andrewsTest <- function(indVar, depVar, piVal=0.15){
	andrewsN <- length(indVar);
	andrewsLowerBound <- as.integer(piVal*andrewsN);
	andrewsUpperBound <- as.integer((1-piVal)*andrewsN);	
	
	andrewskF <- rep(0, andrewsN);
	
	for (i in andrewsLowerBound:andrewsUpperBound){
		localDataFrame <- data.frame(cbind(
			"ind" = c(indVar),
			"dep1" = c(depVar),
			"dep2" = c(rep(0, i), depVar[(i+1):andrewsN])
		));
		
		localConstrained <- lm(ind ~ dep1, data=localDataFrame);
		localUnconstrained <- lm(ind ~ dep1 + dep2, data=localDataFrame);
		
		localSCRC <- sum(localConstrained$residuals^2);
		localSCRU <- sum(localUnconstrained$residuals^2);
		
		localFStat <- (localSCRC - localSCRU)/(localSCRU/(andrewsN-2));
		
		andrewskF[i] <- 2*localFStat;
	}

	return(list(
		"KF" = andrewskF,
		"maxKF" = max(andrewskF),
		"maxKFindex" = which.max(andrewskF)
	));
}




#	EXERCICES

################
#	QUESTION 1	#
################

question1 <- read.dta("index.dta");
q1 <- data.frame(
	cbind(
		"no" = c(question1["no"][[1]]),
		"annee" = c(question1["ann\xe9e"][[1]]),
		"mois" = c(question1["mois"][[1]]),
		"index" = c(question1["index"][[1]])
	)
);
n1 <- length(q1$no);

q1arma <- arima(q1$index, order=c(1,0,0));
q1["uchapArma"] <- q1arma$residuals;

uBar <- mean(q1$uchapArma);
q1sigmaTild <- sqrt( (1/(n1 - 1)) * sum((q1$uchapArma - uBar)^2) );
q1["uTild"] <- (q1$uchapArma - uBar)/q1sigmaTild;

q1kurtosis <- kurtosis(q1$uTild);
print(q1kurtosis);

q1$uchapArmasq <- q1$uchapArma^2;
q1$uchapArmasql1 <- lagBy(q1$uchapArmasq);
engleTest <- lm(uchapArmasq ~ uchapArmasql1, data=q1);
summary(engleTest);

q1garch11 <- garchFit(index ~ garch(1,1), data=q1);


################
#	QUESTION 2	#
################

question2 <- read.dta("macro.dta");
q2 <- data.frame(
	cbind(
		"no" = c(question2["no"][[1]])
	)
);
n2 <- length(q3$no);

################
#	QUESTION 3	#
################

question3 <- read.dta("tauxint.dta");
q3 <- data.frame(
	cbind(
		"no" = c(question3["no"][[1]])
	)
);
n3 <- length(q3$no);

################
#	QUESTION 4	#
################

question4 <- read.dta("tar.dta");
q4 <- data.frame(
	cbind(
		"obs" = c(question4["obs"][[1]]),
		"y" = c(question4["y"][[1]])
	)
);
n4 <- length(q4$obs);
q4["yl1"] <- lagBy(q4$y, 1);

q4AR1 <- lm(y ~ yl1, data=q4);
print(summary(q4AR1));

q4Andrews <- andrewsTest(q4$y, q4$yl1);
plot(1:30, q4Andrews$KF);
print(q4Andrews$maxKFindex);

q4["yl1split"] <- c(rep(0,20), q4$yl1[21:30]);
q4ARandrews <- lm(y ~ yl1 + yl1split, data=q4);
summary(q4ARandrews);
