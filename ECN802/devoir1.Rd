#
#	Université de Sherbrooke
#	ECN802, Hiver 2015
#	Prof: Patrick Richard
#	Devoir 1
#	Cédric Levasseur-Laberge
#

#	IMPORTANT: Indiquer le répertoire de travail avant de commencer

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

x1 <- c(q1$x1);
x2 <- c(q1$x2);
x2l1 <- lagBy(x2,1);
x2l4 <- lagBy(x2,4);
x2l5 <- lagBy(x2,5);
y <- c(q1$y);
yl1 <- lagBy(y,1);
yl4 <- lagBy(y,4);
yl5 <- lagBy(y,5);

#	Retire les 5 premiers éléments
q1Data <- data.frame(
	cbind(
		x1=x1[6:250], 
		x2=x2[6:250], 
		x2l1=x2l1[6:250], 
		x2l4=x2l4[6:250], 
		x2l5=x2l5[6:250], 
		y=y[6:250], 
		yl1=yl1[6:250], 
		yl4=yl4[6:250], 
		yl5=yl5[6:250]
	)
);

#	Estimer le modèle par MCO
q1_mco <- lm(formula = y ~ x2 + yl1 + yl4 + yl5 + x2l1 + x2l4 + x2l5, data = q1Data);
summary(q1_mco);

#	Génère une série temporelle multivariée
q1DataTS <- ts(data.frame(
	cbind(
		x1=x1, 
		x2=x2, 
		x2l1=x2l1, 
		x2l4=x2l4, 
		x2l5=x2l5, 
		y=y, 
		yl1=yl1, 
		yl4=yl4, 
		yl5=yl5
	)
));

#	Estimation du modèle ARMA(1,1)
q1_arima <- arima(ts(y), order=c(1,0,1));
summary(q1_arima);
coef(q1_arima);