#include <Rcpp.h>
#include <ctime>
#include <iostream>
#include <vector>
#include <math.h>
#include <float.h>
using namespace Rcpp;


const double epsilon = std::numeric_limits<double>::min();
const double two_pi = 2.0*3.14159265358979323846;

double * mvrnormC(double mu0, double mu1, double sigma[4]) {
	static double mrn[2];
	double r[2];
	double l[3];

	// Cholesky decomposition sigma=l%%t(l)
	l[0] = sqrt(sigma[0]);
	l[1] = sigma[1] / l[0];
	l[2] = sqrt(sigma[3] - l[1] * l[1]);

	// u1,u2 ~ runif [0,1]
	double u1, u2;
	do
	{
		u1 = rand()*(1.0 / RAND_MAX);
		u2 = rand()*(1.0 / RAND_MAX);
	} while (u1 <= epsilon);

	// Box–Muller transform r[1],r[2]~rnorm[0,1]
	r[0] = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	r[1] = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);

	// mrn=mu+l%*%r
	mrn[0] = mu0 + l[0] * r[0];
	mrn[1] = mu1 + l[1] * r[0] + l[2] * r[1];
	return(mrn);
}

//should work
double lnRGdp(double *y, int yn,double ymax, double Xtemp1, double Xtemp2, double X1, double X2) {
	if (Xtemp1 < (-exp(Xtemp2) / ymax)) {
		return(-DBL_MAX);
	}
	double logLikeXtemp = 0;
	double logLikeX = 0;
	double invExpXtemp2 = exp(-Xtemp2);
	double invExpX2 = exp(-X2);

	//logNormalPdf=log(f(Xtemp1)*f(Xtemp2)/(f(X1)*f(X2))), where f is normal
	//var(x1)=100, var(x2)=1000 =>0.01 and 0.0001, mu(x1)=mu(x2)=0
	double logNormalPdf = -0.5*((Xtemp1 - X1)*(Xtemp1 + X1)*0.01 + (Xtemp2 - X2)*(Xtemp2 + X2)*0.0001);

	// xi=0:  l(sigma)=-k*log(sigma)-sigma^(-1)*sum(yi) p80, Stuart Coles, An Introduction...
	// xi!=0: l(sigma)=-k*log(sigma)-(1+1/xi)*sum(log(1+xi*yi/sigma)) p80, Stuart Coles, An Introduction...
	if (Xtemp1 == 0) {
		for (int i = 0; i < yn; i++) {
			logLikeXtemp += y[i];
		}
		logLikeXtemp = -logLikeXtemp * invExpXtemp2;
		logLikeXtemp -= yn*Xtemp2;
	}
	else {
		for (int i = 0; i < yn; i++) {
			logLikeXtemp += log(1 + Xtemp1*y[i] * invExpXtemp2);
		}
		logLikeXtemp = -logLikeXtemp * (1 + 1 / Xtemp1);
		logLikeXtemp -= yn*Xtemp2;
	}
	if (X1 == 0) {		
		for (int i = 0; i < yn; i++) {
			logLikeX += y[i];
		}
		logLikeX = -logLikeX * invExpX2;
		logLikeX -= yn*X2;
	}
	else {
		for (int i = 0; i < yn; i++) {
			logLikeX += log(1 + X1*y[i] * invExpX2);
		}
		logLikeX = -logLikeX * (1 + 1 / X1);
		logLikeX -= yn*X2;
	}
	return(logLikeXtemp - logLikeX + logNormalPdf);
}




// [[Rcpp::export]]
Rcpp::List mcmcGpdC(NumericVector y, int nstart, int n, NumericVector start, NumericVector muR, NumericMatrix varR, double tau, double a, double gam, double lamda)
{
	srand((unsigned)time(NULL));
	//double total1, total2;
	double mu[2], var[4];
	double lamdaxVar[4];
	double R;
	double gamC = gam;
	double uni;
	double *Xtemp;
	double ymax= -DBL_MAX;
	int yn = y.size();
	double yC[yn];
	//vector (y) to array (yC)
	std::copy(y.begin(), y.end(), yC);
	std::copy(muR.begin(), muR.end(), mu);
	std::copy(varR.begin(), varR.end(), var);
	//double *mrn;
	std::vector<double> thetaXi(n), thetaPhi(n);

	thetaXi[0]= start[0], thetaPhi[0] = start[1];

	//find largest y (ymax)
	for (int j = 0; j < yn; j++) {
		if (ymax < yC[j]) {
			ymax = yC[j];
		}
	}

	if (lamda == 0) {
		lamda = pow(2.38, 2.0) / 2;
	}

	for (int i = 1; i <n; i++) {
		if (tau != 0) {
			gamC = 0.5*exp(-(i + nstart) / tau);// improve 1/sqrt(i)??
		}

		for (int k = 0; k < 4; k++) lamdaxVar[k] = var[k] * lamda;
		Xtemp = mvrnormC(thetaXi[i - 1], thetaPhi[i - 1], lamdaxVar);

		R = lnRGdp(yC, yn, ymax, *(Xtemp), *(Xtemp + 1), thetaXi[i - 1], thetaPhi[i - 1]);
		//R = -0.3566;
		uni = log(rand()*(1.0 / RAND_MAX));
		if (uni < R) {
			thetaXi[i] = *(Xtemp);
			thetaPhi[i] = *(Xtemp+1);
		}
		else {
			thetaXi[i] = thetaXi[i-1];
			thetaPhi[i] = thetaPhi[i-1];
		}
		if (gamC > DBL_EPSILON) {
			if (a != 0) {
				if (R >= 0) {
					R = 1;
				}
				else {
					R = exp(R);
				}
				lamda = lamda*exp(gamC*(R - a));
			}
			// Matrix calculation: var=var+gam*((theta[,i]-mu)%*%t(theta[,i]-mu)-var)
			var[0] = var[0] + gamC*((thetaXi[i] - mu[0])*(thetaXi[i] - mu[0]) -var[0]);
			var[1] = var[1] + gamC*((thetaXi[i] - mu[0])*(thetaPhi[i] - mu[1]) - var[1]);
			var[3] = var[3] + gamC*((thetaPhi[i] - mu[1])*(thetaPhi[i] - mu[1]) - var[3]);
			var[2] = var[1];
			// Matrix calculation: mu<-mu+gam*(theta[,i]-mu)
			mu[0] = mu[0] + gamC*(thetaXi[i] - mu[0]);
			mu[1] = mu[1] + gamC*(thetaPhi[i] - mu[1]);
		}
	}

	std::vector<double> muC(mu, mu + 2);
	std::vector<double> varC(var, var + 4);

	return Rcpp::List::create(Rcpp::Named("data")= NAN, Rcpp::Named("data") = NAN, Rcpp::Named("theta") = NAN, 
		Rcpp::Named("thetaXi") = thetaXi, Rcpp::Named("thetaPhi") = thetaPhi, Rcpp::Named("mu") = muC, 
		Rcpp::Named("var") = varC, Rcpp::Named("n") = n, Rcpp::Named("burnin") = NAN, Rcpp::Named("aRate") = NAN, 
		Rcpp::Named("MLE") = NAN, Rcpp::Named("MLEest") = NAN, Rcpp::Named("tau") = tau, Rcpp::Named("a") = a, 
		Rcpp::Named("gam") = gam, Rcpp::Named("lamda") = lamda);
}

/*
	mrn = mvrnormC(mu[0],mu[1], var);
	total1 = *(mrn);
	total2 = *(mrn+1);
	std::vector<double> mvrnorm(2);
	mvrnorm[0] = total1;
	mvrnorm[1] = total2;

	return Rcpp::List::create(Rcpp::Named("real") = total1, Rcpp::Named("total") = total2, Rcpp::Named("vec1") = mvrnorm,Rcpp::Named("lamda")=lamda);
*/