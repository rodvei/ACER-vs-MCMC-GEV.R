#include <iostream>
#include <ctime>
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
l[2] = sqrt(sigma[3] - l[1]*l[1]);

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
mrn[0] = mu[0]+l[0]*r[0];
mrn[1] = mu[1]+l[1]*r[0]+l[2]*r[1];
return(mrn);
}

int main()
{
	srand((unsigned)time(NULL));
	double mu[2] = {1,2};
	double sigma[4] = { 1,2,2,5 };
	double *mrn;
	cout << time(NULL) << endl;
	mrn = mvrnormC(mu, sigma);
	cout << *(mrn)<<endl;
	cout << *(mrn+1)<<endl<<endl;

	mrn = mvrnormC(mu, sigma);
	cout << *(mrn) << endl;
	cout << *(mrn + 1) << endl;
	cin.get();
	return 0;
}


