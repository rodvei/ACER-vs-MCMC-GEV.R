#include <Rcpp.h>
#include <ctime>
#include <iostream>
using namespace Rcpp;
using namespace std;

const double epsilon = numeric_limits<double>::min();
const double two_pi = 2.0*3.14159265358979323846;

double * mvrnormC(double mu[2], double sigma[4]) {
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
	mrn[0] = mu[0] + l[0] * r[0];
	mrn[1] = mu[1] + l[1] * r[0] + l[2] * r[1];
	return(mrn);
}




/*Rcpp::List sumC(NumericVector x)
{
	int n = x.size();
	double total=0.0;
	for (int i = 0; i < n; i++) {
		total += x[i];
	}

	return Rcpp::List::create(Rcpp::Named("vec")=x, Rcpp::Named("total")= total);
}*/

// var=c(a11,a12,a21,a22)
// [[Rcpp::export]]
Rcpp::List mcmcgpdC(int n, NumericVector muR, NumericMatrix varR) //,NumericVector y, NumericVector start, double tau, double a, double gamma, double lamda)
{
	srand((unsigned)time(NULL));
	double total1, total2;
	double mu[2];
	copy(muR.begin(), muR.end(), mu);
	double var[4] = { varR(0,0),varR(0,1), varR(1,0),varR(1,1) };
	double *mrn;
	mrn = mvrnormC(mu, var);
	total1 = *(mrn);
	total2 = *(mrn+1);
	return Rcpp::List::create(Rcpp::Named("real") = total1, Rcpp::Named("total") = total2);
}