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
library("sandwich");
library("fGarch");
library("vars");

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

#	Définit la méthode de Chan
chanMethod <- function(indVar, depVar, piVal=0.15){
	localOrder <- order(indVar);
	localInd <- indVar[localOrder];
	localDep <- depVar[localOrder];
	
	chanN <- length(localInd);
	chanLowerBound <- as.integer(piVal*chanN);
	chanUpperBound <- as.integer((1-piVal)*chanN);
	
	chanSCR <- rep(NA, chanN);
	chanModels <- rep(NA, chanN);
	
	for (i in chanLowerBound:chanUpperBound){
		
		It <- as.numeric(localDep > localDep[i]);
		Itinv <- as.numeric(It == 0);

		localDataFrame <- data.frame(cbind(
			"ind" = c(localInd),
			"dep1" = c(localDep * It),
			"dep2" = c(localDep * Itinv)
		));
		
		localModel <- lm(ind ~ dep1 + dep2, data=localDataFrame);
		chanModels[i] <- localModel;
		chanSCR[i] <- sum(localModel$residuals^2);
	}
	
	return(list(
		"SCR" = chanSCR,
		"minSCR" = min(chanSCR, na.rm=TRUE),
		"threshold" = localDep[which.min(chanSCR)],
		"minSCRindex" = which.min(chanSCR),
		"model" = chanModels[which.min(chanSCR)]
	));
}



#	EXERCICES

#################
#	QUESTION 1	#
#################

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

#	4e moment non-conditionnel
q1arma <- arima(q1$index, order=c(1,0,0));
q1["uchapArma"] <- q1arma$residuals;

uBar <- mean(q1$uchapArma);
q1sigmaTild <- sqrt( (1/(n1 - 1)) * sum((q1$uchapArma - uBar)^2) );
q1["uTild"] <- (q1$uchapArma - uBar)/q1sigmaTild;

q1kurtosis <- kurtosis(q1$uTild);
print(q1kurtosis);

#	Test de Engle
q1$uchapArmasq <- q1$uchapArma^2;
q1$uchapArmasql1 <- lagBy(q1$uchapArmasq);
engleTest <- lm(uchapArmasq ~ uchapArmasql1, data=q1);
summary(engleTest);

#	ARMA-GARCH
q1garch11 <- garchFit(index ~ garch(1,1), data=q1);
q1["epsilonGarch11"] <- c(q1garch11@residuals);
q1["sigmaGarch11"] <- c(q1garch11@sigma.t);
q1["zGarch11"] <- c(q1$epsilonGarch11/q1$sigmaGarch11);
print(mean(q1$zGarch11));
print(var(q1$zGarch11));
hist(q1$zGarch11);

 
#################
#	QUESTION 2	#
#################

question2 <- read.dta("macrov12.dta");
q2 <- data.frame(
	cbind(
		"y" = c(question2["y"][[1]]),
		"r" = c(question2["r"][[1]]),
		"p" = c(question2["p"][[1]])
	)
);
n2 <- length(q2$y);

#	Trouver la value optimale de p
maxP <- 12;
aics <- rep(NA, maxP);
bics <- rep(NA, maxP);
SCRs <- rep(NA, maxP);

for (i in 4:maxP) {
	localq2var <- VAR(q2, p=i);
	SCRs[i] <- sum(localq2var$varresult$y$residuals^2 + localq2var$varresult$r$residuals^2 + localq2var$varresult$p$residuals^2)
	aics[i] <- aic((SRCs[i]/n2), n2, ((3*i)+1));
	bics[i] <- bic((SRCs[i]/n2), n2, ((3*i)+1));
}
print(which.min(aics));
print(which.min(bics));
print(which.min(SCRs));

#	Tests de causalité de Granger
goodVar <- VAR(q2, p=5);
causY <- causality(goodVar, cause="y");
print(causY$Granger);
causR <- causality(goodVar, cause="r");
print(causR$Granger);
causP <- causality(goodVar, cause="p");
print(causP$Granger);



#################
#	QUESTION 4	#
#################

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

#	Applique la méthode de Andrews
q4Andrews <- andrewsTest(q4$y, q4$yl1);
plot(1:30, q4Andrews$KF);
print(q4Andrews$maxKFindex);

q4["yl1split"] <- c(rep(0,20), q4$yl1[21:30]);
q4ARandrews <- lm(y ~ yl1 + yl1split, data=q4);
summary(q4ARandrews);

#	Applique la méthode de Chan
q4yorder <- order(q4$y);
q4["ysorted"] <- q4$y[q4yorder];
q4["yl1sorted"] <- q4$yl1[q4yorder];

q4Chan <- chanMethod(q4$y[2:30], q4$yl1[2:30]);
plot(q4$y[2:30], q4Chan$SCR);
print(q4Chan$threshold);

q4It <- c(0, as.numeric(q4$yl1[2:30] > q4Chan$threshold));
q4Itinv <- as.numeric(q4It==0);
q4["yl1Chan1"] <- q4$yl1 * q4It;
q4["yl1Chan2"] <- q4$yl1 * q4Itinv;

q4modelChan <- lm(y ~ yl1Chan1 + yl1Chan2, data=q4);
summary(q4modelChan);