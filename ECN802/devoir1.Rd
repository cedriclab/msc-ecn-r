#
#	Université de Sherbrooke
#	ECN802, Hiver 2015
#	Prof: Patrick Richard
#	Devoir 1
#	Cédric Levasseur-Laberge
#

#	IMPORTANT: Indiquer le répertoire de travail avec setwd("path/to/directory") avant de commencer

library(foreign);

#	Définit une fonction locale pour les retards
lagBy <- function(x,k=1){
	lagNAs <- rep(NA, k);
	laggedVector <- x[1:(length(x)-k)]
	
	c(lagNAs, laggedVector)
}


#################
#	QUESTION 1	#
#################

q1 <- read.dta("question1.dta");
q1N <- length(q1$x1);

q1Data <- data.frame(
	cbind(
		x1=c(q1$x1), 
		x2=c(q1$x2), 
		x2l1=lagBy(q1$x2,1), 
		x2l4=lagBy(q1$x2,4), 
		x2l5=lagBy(q1$x2,5), 
		y=c(q1$y), 
		yl1=lagBy(q1$y,1), 
		yl4=lagBy(q1$y,4), 
		yl5=lagBy(q1$y,5)
	)
);

#	Modèle non-linéaire

nonLinearFunction <- function(dataFrame, a1, a4, b1, b2){
	ychap <- c(((1-a1-a4-a1*a4)*b1) + dataFrame["x2"]*b2 + a1*(dataFrame["yl1"] - b2*dataFrame["x2l1"]) + a4*(dataFrame["yl4"] - b2*dataFrame["x2l4"]) + a1*a4*(dataFrame["yl5"] + b2*dataFrame["x2l5"]));
};
nonLinearModel <- y ~ ((1-a1-a4-a1*a4)*b1) + x2*b2 + a1*(yl1 - b2*x2l1) + a4*(yl4 - b2*x2l4) + a1*a4*(yl5 + b2*x2l5);
startValues <- c(a1=0, a4=0, b1=0, b2=0);

q1nlm <- nls(formula=nonLinearModel, data=q1Data, start=startValues, na.action=na.omit, trace=TRUE);

#	Retire les 5 premiers éléments
q1DataNoNA <- data.frame(
	cbind(
		x1=q1Data$x1[6:250], 
		x2=q1Data$x2[6:250], 
		x2l1=q1Data$x2l1[6:250], 
		x2l4=q1Data$x2l4[6:250], 
		x2l5=q1Data$x2l5[6:250], 
		y=q1Data$y[6:250], 
		yl1=q1Data$yl1[6:250], 
		yl4=q1Data$yl4[6:250], 
		yl5=q1Data$yl5[6:250]
	)
);

#	Estimer le modèle par MCO

linearModel <- y ~ x2 + yl1 + yl4 + yl5 + x2l1 + x2l4 + x2l5;
q1Mco <- lm(formula=linearModel, data=q1Data);
q1McoCoefs <- coef(q1Mco);
summary(q1Mco);

#	Relance l'estimation du modèle non-linéaire avec les coefficients des MCOs comme valeurs de départ
alpha1 <- 0.5092815;
alpha4 <- -0.1536069;
beta1 <- 0.5966792;
beta2 <- 0.5416752;
startValues <- c(a1=alpha1, a4=alpha4, b1=beta1, b2=beta2);
q1nlm2 <- nls(formula=nonLinearModel, data=q1Data, start=startValues, na.action=na.omit, trace=TRUE);
summary(q1nlm2);

#	Génère les données pour la régression artificielle
q1Data$ychapArt <- rep(NA, q1N);
q1Data$uchapArt <- rep(NA, q1N);
q1Data$bb2Art <- rep(NA, q1N);			# la valeur artificielle qui a pour coefficient b_beta2
q1Data$ba1Art <- rep(NA, q1N);			# la valeur artificielle qui a pour coefficient b_alpha1
q1Data$ba4Art <- rep(NA, q1N);			# la valeur artificielle qui a pour coefficient b_alpha4

for (t in 6:q1N) {
	
	q1Data$ychapArt[t] <- nonLinearFunction(q1Data[t, c("x2", "x2l1", "x2l4", "x2l5", "yl1", "yl4", "yl5")], alpha1, alpha4, beta1, beta2)[[1]];
	q1Data$uchapArt[t] <- q1Data$ychapArt[t] - q1Data$y[t];
	
	q1Data$bb2Art[t] <- q1Data$x2[t] - alpha1*q1Data$x2l1[t] - alpha4*q1Data$x2l4[t] + alpha1 * alpha4 * q1Data$x2l5[t];
	q1Data$ba1Art[t] <- (alpha4 - 1)*beta1 + q1Data$yl1[t] - beta2*q1Data$x2l1[t] - alpha4*q1Data$yl5[t] + alpha4*beta2*q1Data$x2l5[t];
	q1Data$ba4Art[t] <- (alpha1 - 1)*beta1 + q1Data$yl4[t] - beta2*q1Data$x2l4[t] - alpha1*q1Data$yl5[t] + alpha1*beta2*q1Data$x2l5[t];
}

artReg <- lm(uchapArt ~ bb2Art + ba1Art + ba4Art, data=q1Data);
summary(artReg);


#################
#	QUESTION 3	#
#################

q3 <- read.dta("count.dta");
q3N <- length(q3$med);

q3Data <- data.frame(q3);
poisson <- glm(med ~ child + access + health, data=q3Data, family="poisson");
summary(poisson);

q3Data$ychapPoisson <- c(predict(poisson));
q3Data$uchapPoisson <- rep(0, q3N);

for (i in 1:q3N){
	q3Data$uchapPoisson[i] <- q3Data$ychapPoisson[i] - q3Data$med[i];
}


#################
#	QUESTION 4	#
#################

q4x <- c(85, 76, 58, 38, 47, 68, 75, 21, 60, 61);
q4mu <- mean(q4x);
q4sigma2 <- var(q4x);