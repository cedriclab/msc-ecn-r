#
#	Université de Sherbrooke
#	ECN827, Hiver 2015
#	Prof: Patrick Richard
#	Devoir 1
#	Cédric Levasseur-Laberge
#

#	IMPORTANT: Indiquer le répertoire de travail avec setwd("path/to/directory") avant de commencer


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


#	Utilisation de la fonction ar() pour déterminer l'ordre p du processus autorégressif
q3AR <- ar(y);
print(q3AR);

#	Définit une fonction qui résout le système d'équations permettant de trouver la valeur de alpha-chapeau et theta-chapeau à partir des beta-chapeau de chaque modèle AR
getAlphaThetaHats <- function(beta1, beta2){
	thetaChap <- (beta2*(-1))/beta1;
	alphaChap <- beta1 - thetaChap;
	
	print(paste("alphaChap: ", alphaChap, " | thetaChap: ", thetaChap));
}

#	Processus AR et ARMA

#	Définit une fonction générale
compareARandARMA <- function(p,size){
	q3ar <- ar(x = y[1:size], aic = FALSE, order.max = p);
	q3ARMA <- arima(y[1:size], order=c(p,0,1));
	
	print(paste("ARMA(", p, ",1 ): "));
	print(paste("alphaChap: ", coef(q3ARMA)["ar1"], " | thetaChap: ", coef(q3ARMA)["ma1"]));
	print(paste("AR(", p, "): "));
	getAlphaThetaHats(q3ar$ar[1], q3ar$ar[2]);
	print(" ");
}

compareARandARMA(5,50);
compareARandARMA(9,100);
compareARandARMA(16,200);
compareARandARMA(22,300);
compareARandARMA(26,400);




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
q5N <- length(q5$y);

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

#	Vérifie la forme du processus autorégressif pour Y
arTestY <- ar(x = q5$deltaY, order.max=10);
print(arTestY);

#	Génère les lags pour Y
q5$deltaYl1 <- lagBy(q5$deltaY, 1);
q5$deltaYl2 <- lagBy(q5$deltaY, 2);
q5$deltaYl3 <- lagBy(q5$deltaY, 3);
q5$deltaYl4 <- lagBy(q5$deltaY, 4);
q5$deltaYl5 <- lagBy(q5$deltaY, 5);
q5$deltaYl6 <- lagBy(q5$deltaY, 6);
q5$deltaYl7 <- lagBy(q5$deltaY, 7);
q5$deltaYl8 <- lagBy(q5$deltaY, 8);

#	Modèle sans tendance déterministe
adfMCOY <- lm(deltaY ~ Yl1 + deltaYl1 + deltaYl2 + deltaYl3 + deltaYl4 + deltaYl5 + deltaYl6 + deltaYl7 + deltaYl8, data=q5);
summary(adfMCOY);

#	Modèle avec tendance déterministe
adfDetMCOY <- lm(deltaY ~ t + Yl1 + deltaYl1 + deltaYl2 + deltaYl3 + deltaYl4 + deltaYl5 + deltaYl6 + deltaYl7 + deltaYl8, data=q5);
summary(adfDetMCOY);

#	Vérifie la forme du processus autorégressif pour Y
arTestX <- ar(x = q5$deltaX, order.max=10);
print(arTestX);

#	Génère les lags pour X
q5$deltaXl1 <- lagBy(q5$deltaX, 1);
q5$deltaXl2 <- lagBy(q5$deltaX, 2);
q5$deltaXl3 <- lagBy(q5$deltaX, 3);
q5$deltaXl4 <- lagBy(q5$deltaX, 4);
q5$deltaXl5 <- lagBy(q5$deltaX, 5);
q5$deltaXl6 <- lagBy(q5$deltaX, 6);
q5$deltaXl7 <- lagBy(q5$deltaX, 7);
q5$deltaXl8 <- lagBy(q5$deltaX, 8);

#	Modèle sans tendance déterministe
adfMCOX <- lm(deltaX ~ Xl1 + deltaXl1 + deltaXl2 + deltaXl3 + deltaXl4 + deltaXl5 + deltaXl6 + deltaXl7 + deltaXl8, data=q5);
summary(adfMCOX);

#	Modèle avec tendance déterministe
adfDetMCOX <- lm(deltaY ~ t + Xl1 + deltaXl1 + deltaXl2 + deltaXl3 + deltaXl4 + deltaXl5 + deltaXl6 + deltaXl7 + deltaXl8, data=q5);
summary(adfDetMCOX);



#	Estimation simple du modèle par les MCO
basicMCO <- lm(y ~ t + x, data=q5);
summary(basicMCO);

for (i in 1:200){
	q5$ychap[i] <- -1.061e+01 - 1.179e-01*i + 1.665e+00*q5$x[i];
	q5$uchap[i] <- q5$ychap[i] - q5$y[i];
	q5$deltaUchap[i] <- q5$uchap[i] - q5$uchap[i-1];
} 

#	Vérifie la forme du processus autorégressif pour û
arTestUchap <- ar(x = q5$deltaUchap, order.max=10);
print(arTestUchap);

q5$uchapl1 <- lagBy(q5$uchap, 1);

#	Génère les lags pour delta û
q5$deltaUchapl1 <- lagBy(q5$deltaUchap, 1);
q5$deltaUchapl2 <- lagBy(q5$deltaUchap, 2);
q5$deltaUchapl3 <- lagBy(q5$deltaUchap, 3);
q5$deltaUchapl4 <- lagBy(q5$deltaUchap, 4);
q5$deltaUchapl5 <- lagBy(q5$deltaUchap, 5);
q5$deltaUchapl6 <- lagBy(q5$deltaUchap, 6);
q5$deltaUchapl7 <- lagBy(q5$deltaUchap, 7);
q5$deltaUchapl8 <- lagBy(q5$deltaUchap, 8);

#	Modèle avec tendance déterministe
adfDetMCOUchap <- lm(deltaUchap ~ t + uchapl1 + deltaUchapl1 + deltaUchapl2 + deltaUchapl3 + deltaUchapl4 + deltaUchapl5 + deltaUchapl6 + deltaUchapl7 + deltaUchapl8, data=q5);
summary(adfDetMCOUchap);


#	Modèle à correction d'erreur

q5EcmLm <- y ~ t + Yl1 + x;
q5ECM <- lm(q5EcmLm, data=q5);
summary(q5ECM);

#	Test de Breusch-Godfrey

q5$ychapECM1 <- c(0, predict(q5ECM));
q5$uchapECM1 <- rep(0, q5N);
for (t in 1:q5N) {
	q5$uchapECM1[t] <- q5$ychapECM1[t] - q5$y[t];
}
q5$uchapECM1l1 <- lagBy(q5$uchapECM1, 1);

q5BG <- lm(uchapECM1 ~ t + Yl1 + x + uchapECM1l1, data=q5);
summary(q5BG);

#	Régression de Y sur X
q5Easy <- lm(y ~ x, data=q5);
summary(q5Easy);

q5$ychapEasy <- predict(q5Easy);
q5$uchapEasy <- rep(0, q5N);
for (t in 1:q5N) {
	q5$uchapEasy[t] <- q5$ychapEasy[t] - q5$y[t];
}
q5$uchapEasyl1 <- lagBy(q5$uchapEasy, 1);

q5BG2 <- lm(uchapEasy ~ x + uchapEasyl1, data=q5);
summary(q5BG2);

q5MCOdet <- lm(y ~ t + x, data=q5);
summary(q5MCOdet);

q5$ychapMCOdet <- predict(q5MCOdet);
q5$uchapMCOdet <- rep(0, q5N);
for (t in 1:q5N) {
	q5$uchapMCOdet[t] <- q5$ychapMCOdet[t] - q5$y[t];
}
q5$uchapMCOdetl1 <- lagBy(q5$uchapMCOdet, 1);

q5BG3 <- lm(uchapMCOdet ~ x + t + uchapMCOdetl1, data=q5);
summary(q5BG3);