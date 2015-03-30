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

question2 <- read.dta("question2.dta");
q2 <- data.frame(
	cbind(
		"obs"=c(question2["obs"][[1]])
	)
);
q2n <- length(q2$obs);

#################
#	QUESTION 3	#
#################

question3 <- read.dta("question3.dta");
q3 <- data.frame(
	cbind(
		"obs"=c(question3["obs"][[1]]),
		"x"=c(question3["x"][[1]]),
		"y"=c(question3["y"][[1]])
	)
);
q3n <- length(q3$obs);