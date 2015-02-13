#
#	Université de Sherbrooke
#	ECN827, Hiver 2015
#	Prof: Patrick Richard
#	Devoir 1
#	Cédric Levasseur-Laberge
#

#	IMPORTANT: Indiquer le répertoire de travail avec setwd("path/to/directory") avant de commencer
#setwd("../../ECN827/Devoir 1");

library(foreign);

#	UTILITAIRES

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



#	EXERCICES

#################
#	QUESTION 3	#
#################

question3 <- read.dta("question3.dta");
y <- question3$y;

#	Estimation du modèle ARMA(1,1)
q3_AR_n <- arima(y, order=c(1,0,1));
summary(q3_AR_n);
coeffs <- coef(q3_AR_n);
alpha <- coeffs["ar1"];
theta <- coeffs["ma1"];
eps0 <- coeffs["intercept"];
epsilons <- q3_AR_n$residuals;

#	Utilise les coefficients obtenus dans la dernière régression pour 
#	obtenir les valeurs prédites de Y via le processus AR
ychap <- rep(0, length(y));
ychap[1] <- y[1];
for (t in 2:length(y)) {
	for (j in 1:(t-1)) {
		ychap[t] <- ychap[t] + ((-1)^(j+1) * theta^(j-1) * (alpha + theta) * y[t-j]);
	}
	ychap[t] <- ychap[t] + ((-1)^(t+1) * ((theta^(t+1) * alpha * y[1]) + (theta^t * epsilons[1])));
}

q3ar5 <- ar(x = y[1:50], aic = FALSE, order.max = 5);
q3ar9 <- ar(x = y[1:100], aic = FALSE, order.max = 9);
q3ar16 <- ar(x = y[1:200], aic = FALSE, order.max = 16);
q3ar22 <- ar(x = y[1:300], aic = FALSE, order.max = 22);
q3ar26 <- ar(x = y[1:400], aic = FALSE, order.max = 26);



#################
#	QUESTION 4	#
#################

question4 <- read.dta("question4.dta");
y <- question4$y;
x <- question4$x;
n <- length(y);

#	Définit les lags à utiliser lors de la régression par les MCOs
question4["y_l1"] <- lagBy(y, 1);

question4["x_l1"] <- lagBy(x, 1);
question4["x_l2"] <- lagBy(x, 2);
question4["x_l3"] <- lagBy(x, 3);
question4["x_l4"] <- lagBy(x, 4);

question4["ychap_dl4"] <- rep(0, n);
question4["ychap_adl11"] <- rep(0, n);
question4["ychap_adl12"] <- rep(0, n);

question4["uchap_dl4"] <- rep(0, n);
question4["uchap_adl11"] <- rep(0, n);
question4["uchap_adl12"] <- rep(0, n);

#	Modèle DL(4)
dl4 <- lm(y ~ x_l1 + x_l2 + x_l3 + x_l4, data=question4)
dl4_coefs <- coef(dl4);
sigma2_dl4 <- 0;

for (t in 5:n) {
	question4$ychap_dl4[t] <- dl4_coefs["(Intercept)"] + (dl4_coefs["x_l1"] * question4$x_l1[t]) + (dl4_coefs["x_l2"] * question4$x_l2[t]) + (dl4_coefs["x_l3"] * question4$x_l3[t]) + (dl4_coefs["x_l4"] * question4$x_l4[t]);
	question4$uchap_dl4[t] <- question4$ychap_dl4[t] - y[t];
	sigma2_dl4 <- sigma2_dl4 + (question4$uchap_dl4[t] * question4$uchap_dl4[t]);
}
sigma2_dl4 <- sigma2_dl4/(n-4-2);
aic_dl4 <- aic(sigma2_dl4, (n-4), 6);
bic_dl4 <- bic(sigma2_dl4, (n-4), 6);
dl4_res <- c(s2=sigma2_dl4, aic=aic_dl4, bic=bic_dl4);

#	Modèle ADL(1,1)
adl11 <- lm(y ~ x_l1 + y_l1, data=question4)
adl11_coefs <- coef(adl11);
sigma2_adl11 <- 0;

for (t in 2:n) {
	question4$ychap_adl11[t] <- adl11_coefs["(Intercept)"] + (adl11_coefs["x_l1"] * question4$x_l1[t]) + (adl11_coefs["y_l1"] * question4$y_l1[t]);
	question4$uchap_adl11[t] <- question4$ychap_adl11[t] - y[t];
	sigma2_adl11 <- sigma2_adl11 + (question4$uchap_adl11[t] * question4$uchap_adl11[t]);
}
sigma2_adl11 <- sigma2_adl11/(n-1-2);
aic_adl11 <- aic(sigma2_adl11, (n-1), 4);
bic_adl11 <- bic(sigma2_adl11, (n-1), 4);
adl11_res <- c(s2=sigma2_adl11, aic=aic_adl11, bic=bic_adl11);

#	Modèle ADL(1,2)
adl12 <- lm(y ~ x_l1 + x_l2 + y_l1, data=question4);
adl12_coefs <- coef(adl12);
sigma2_adl12 <- 0;

for (t in 3:n) {
	question4$ychap_adl12[t] <- adl12_coefs["(Intercept)"] + (adl12_coefs["x_l1"] * question4$x_l1[t]) + (adl12_coefs["x_l2"] * question4$x_l2[t]) + (adl12_coefs["y_l1"] * question4$y_l1[t]);
	question4$uchap_adl12[t] <- question4$ychap_adl12[t] - y[t];
	sigma2_adl12 <- sigma2_adl12 + (question4$uchap_adl12[t] * question4$uchap_adl12[t]);
}
sigma2_adl12 <- sigma2_adl12/(n-2-2);
aic_adl12 <- aic(sigma2_adl12, (n-2), 5);
bic_adl12 <- bic(sigma2_adl12, (n-2), 5);
adl12_res <- c(s2=sigma2_adl12, aic=aic_adl12, bic=bic_adl12);

q4Models <- c(dl4=dl4_res, adl11=adl11_res, adl12=adl12_res);
print(q4Models);

#	Tests de Diebold-Mariano

forecastErrors <- rep(NA, 20);
DM_d <- rep(0, 20);

for (i in 1:20) {
	forecastIndex <- n-20+i;
	
	forecastErrorDL4 <- question4$y[forecastIndex] - question4$ychap_dl4[forecastIndex];
	forecastErrorADL12 <- question4$y[forecastIndex] - question4$ychap_adl12[forecastIndex];
	
	forecastErrors[i] <- c(dl4=forecastErrorDL4, adl12=forecastErrorADL12);
	
	DM_d[i] <- forecastErrorDL4^2 - forecastErrorADL12^2;
}

DM_dBar <- mean(DM_d);
DMStat <- DM_dBar/sqrt(var(DM_d));
print(DMStat);



#################
#	QUESTION 5	#
#################
question5 <- read.dta("question5.dta");

q5 <- data.frame(
	cbind(
		x=question5$x,
		Xl1=lagBy(question5$x, 1),
		y=question5$y,
		Yl1=lagBy(question5$y, 1),
		deltaY=rep(0, length(question5$y)),
		deltaX=rep(0, length(question5$y)),
		t=rep(0, length(question5$y))
	)
);
#	Génère les delta X et delta Y
for (i in 2:length(q5$y)) {
	q5$t[i] <- i;
	q5$deltaY[i] <- q5$y[i] - q5$y[i-1];
	q5$deltaX[i] <- q5$x[i] - q5$x[i-1];
}
q5$deltaYl1 <- lagBy(q5$deltaY, 1);

#	Modèle sans tendance déterministe
adfMCO <- lm(deltaY ~ Yl1 + deltaYl1, data=q5);
summary(adfMCO);

#	Modèle avec tendance déterministe
adfDetMCO <- lm(deltaY ~ t + Yl1 + deltaYl1, data=q5);
summary(adfDetMCO);

#	Estimation simple du modèle par les MCO
basicMCO <- lm(y ~ t + x, data=q5);
summary(basicMCO);

for (i in 1:200){
	q5$ychap[i] <- -1.061e+01 - 1.179e-01*i + 1.665e+00*q5$x[i];
	q5$uchap[i] <- q5$ychap[i] - q5$y[i];
	q5$deltaUchap[i] <- q5$uchap[i] - q5$uchap[i-1];
} 