#
#	Université de Sherbrooke
#	ECN827, Hiver 2015
#	Prof: Patrick Richard
#	Devoir 1
#	Cédric Levasseur-Laberge
#

#	IMPORTANT: Indiquer le répertoire de travail avant de commencer

library(foreign);


#	Question 3
question3 <- read.dta("question3.dta");
y <- question3$y;

q3_AR_n <- arima(y, order=c(1,0,1));
summary(q3_AR_n);
coeffs <- coef(q3_AR_n);
alpha <- coeffs["ar1"];
theta <- coeffs["ma1"];

#	Utilise les coefficients obtenus dans la dernière régression pour 
#	obtenir les valeurs prédites de Y via le processus AR
ychap <- rep(0, length(y));
for (i in 1:length(y)) {
	if (i==1) {
		ychap[i] <- y[i];
	} else {
		epsilon <- ychap[i-1] - y[i-1];
		ychap[i] <- alpha^(i-1) * theta * epsilon;
	}
}


arimaFromSampleSize <- function(sample, size, p){
	localRegResult <- arima(sample[1:size], order=c(p,0,1));
	paste("Sample size is: ", size, " and AR process is of p=", p, " order.")
	print(coef(localRegResult))
}

sizes <- c(50, 100, 200, 300, 400);
pOrders <- c(5, 9, 16, 22, 26);
for (i in 1:5){arimaFromSampleSize(y, sizes[i], pOrders[i]);}


#	Question 4
question4 <- read.dta("question4.dta");

#	Question 5
#question5 <- read.dta("question5.dta");
