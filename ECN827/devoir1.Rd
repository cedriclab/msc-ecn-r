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
coef(q3_AR_n);

arimaFromSampleSize <- function(sample, size, p){
	localRegResult <- arima(sample[1:size], order=c(p,0,1));
	paste("Sample size is: ", size, " and AR process is of p=", p, " order.")
	print(coef(localRegResult))
}

sizes <- c(50, 100, 200, 300, 400);
pOrders <- c(5, 9, 16, 22, 26);
for (i in 1:5){arimaFromSampleSize(y, sizes[i], pOrders[i]);}


#	Question 4
#question4 <- read.dta("question4.dta");

#	Question 5
#question5 <- read.dta("question5.dta");
