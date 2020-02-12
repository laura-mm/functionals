// lotka-volterra dynamics with eigen array class, RK4
// modified for condor
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <fstream>
#include <omp.h>
#include <eigen3/Eigen/Dense>
#include <iomanip>
using namespace Eigen;
using namespace std;

random_device generator;
mt19937 twist(generator());

double a;
int N;
int T;
double dt;
double mu;
IOFormat p(20, DontAlignCols, ",", ",", "", "", "", "");
IOFormat c(6, DontAlignCols, ",", ",", "", "", "", "");


double norm()
{
	normal_distribution<double> distributionn(0.0, 1.0);
	double num = distributionn(twist);
	return num;
}
double rando()
{
	uniform_real_distribution<double> distributionr(0.0, 1.0);
	double num = distributionr(twist);
	return num;
}


//ArrayXd functiom(MatrixXd m) {return m.array();}

ArrayXd functiom(MatrixXd m) //piecewise
{
	ArrayXd v = m.array();
	for (int i = 0; i < v.size(); i++)
	{
		if (v(i) >= a) v(i) = a;
		else if (v(i) <= -a) v(i) = -a;
	}
	return v;
}


/*ArrayXd functiom(MatrixXd m) //holling
{
	ArrayXd v = m.array();
	return (2.0*a*v)/(a + (2.0*abs(v)));
}*/

class simulation
{
	public:
	MatrixXd A;
	ArrayXd x; ArrayXd y;
	bool fixed; bool unique; bool diverge;
	vector<ArrayXd> trajx; vector<ArrayXd> trajy; //
	simulation(double sigma, double gamma)
	{
		diverge = false;
		x = ArrayXd(N); y = ArrayXd(N);
		for (int i = 0; i < N; i++) {x(i) = rando(); y(i) = rando();}
		//trajx.push_back(x); trajy.push_back(y); //
		MatrixXd AA(N,N);
		AA(0,0) = 0.0;
		for (int i = 1; i < N; i++)
		{
			AA(i,i) = 0.0;
			for (int j = 0; j < i; j++)
			{
				double z1 = norm();
				double z2 = norm();
				double y1 = (gamma*z1) + (sqrt(1.0-(gamma*gamma))*z2);
				AA(i,j) = (mu/N) + ((sigma*z1)/sqrt(N));
				AA(j,i) = (mu/N) + ((sigma*y1)/sqrt(N));
			}
		}
		A = AA;
	}
	~simulation() {}
	void check() // this is for traj only last 1% and this does each species induvidually!!
	{
		fixed = true;
		for (int t = 2; t < (0.01*T)+2; t++)
		{
			if (((abs(trajx[t] - trajx[t-1]) < 0.0001) || (abs((trajx[t] - trajx[t-1])/(trajx[t-1] - trajx[t-2])) < 1.0)).minCoeff() == 0) {fixed = false; break;}
			if (((abs(trajy[t] - trajy[t-1]) < 0.0001) || (abs((trajy[t] - trajy[t-1])/(trajy[t-1] - trajy[t-2])) < 1.0)).minCoeff() == 0) {fixed = false; break;}
		}
		unique = true;
		for (int t = 2; t < (0.01*T)+2; t++)
		{
			if (((abs(trajx[t] - trajy[t]) < 0.0001) || ((abs(trajx[t] - trajy[t])/abs(trajx[t-1] - trajy[t-1])) < 1.0)).minCoeff() == 0) {unique = false; break;}
		}
	}
	void run()
	{
		for (int t = 0; t < T; t++)
		{
			ArrayXd k1 = dt*x*(1.0 - x + functiom(A*x.matrix())); //
			ArrayXd k2 = dt*(x + (k1/2.0))*(1.0 - (x + (k1/2.0)) + functiom(A*(x + (k1/2.0)).matrix())); //
			ArrayXd k3 = dt*(x + (k2/2.0))*(1.0 - (x + (k2/2.0)) + functiom(A*(x + (k2/2.0)).matrix())); //
			ArrayXd k4 = dt*(x + k3)*(1.0 - (x + k3) + functiom(A*(x + k3).matrix())); //
			x += (k1 + (2.0*k2) + (2.0*k3) + k4)/6.0;

			k1 = dt*y*(1.0 - y + functiom(A*y.matrix())); //
			k2 = dt*(y + (k1/2.0))*(1.0 - (y + (k1/2.0)) + functiom(A*(y + (k1/2.0)).matrix())); //
			k3 = dt*(y + (k2/2.0))*(1.0 - (y + (k2/2.0)) + functiom(A*(y + (k2/2.0)).matrix())); //
			k4 = dt*(y + k3)*(1.0 - (y + k3) + functiom(A*(y + k3).matrix())); //
			y += (k1 + (2.0*k2) + (2.0*k3) + k4)/6.0;

			if (x.maxCoeff() > pow(10.0, 20.0) || y.maxCoeff() > pow(10.0, 20.0)) {diverge = true; break;} // linear only

			if (t >= (0.99*T)-2) {trajx.push_back(x); trajy.push_back(y);} //

		}
		check();
	}
};


void plot(double sigma, double gamma, int plots)
{

	ofstream datax;
	datax.open ("trajectoriesx.txt");
	ofstream datay;
	datay.open ("trajectoriesy.txt");
	
	simulation sim(sigma, gamma);
	sim.run();
	cout << sim.fixed << endl;
	cout << sim.unique << endl;
		
	for (int i = 0; i < plots; i++)
	{
		for (int t = 0; t < 0.01*T; t++)
		{
			datax << sim.trajx[t](i);
			datay << sim.trajy[t](i);
			if (i != plots-1 || t != T-1) {datax << ","; datay << ",";}
		}
	}
}

int main(int argc, char** argv)	
{
	if (argc != 3) {return 1;}

	mu = 0.0;
	N = 200;
	T = 200000;
	dt = 0.001;
	a = 2.0;


	int v1 = atoi(argv[1]);
	int v2 = atoi(argv[2]);


	string filename = "p2_" + to_string(v1) + "_" + to_string(v2) + ".txt"; //
	ofstream file; file.open(filename);

	//file << v1 << ", " << v2 << endl;

	double gamma = (0.02*(double)v1) - 1.0;
	double po = (0.04*(double)v2) - 1.0;
	double sigma = pow(10.0, po);


	int count0 = 0;
	int count1 = 0;
	int count2 = 0;

	for (int r = 0; r < 200; r ++)
	{
		simulation sim(sigma, gamma);
		sim.run();
		if (sim.diverge == false) // linear only
		{ //
		if (sim.fixed) {if (sim.unique) count0 ++; else count1 ++;}
		else count2 ++;
		} //
		//cout << sim.fixed << sim.unique << endl;
	}

	//file << count0 << ", " << count1 << ", " << count2 << endl;

	file << (double)count0/200.0 << "," << (double)count1/200.0 << "," << (double)count2/200.0;

	file.close();


	return 0;
}

/*int main()	
{
	mu = 0.0;
	a = 0.5;
	N = 200;
	T = 200000;
	dt = 0.001;
	int grid = 100;
	int runs = 100;
	int plots = 20;


	simulation sim(100, 0.75); // sigma, gamma
	sim.run();
	cout << sim.fixed << endl;
	cout << sim.unique << endl;

	int v1 = 50; // gamma
	int v2 = 30; // sigma

	double gamma = (0.02*(double)v1) - 1.0;
	double po = (0.04*(double)v2) - 1.0;
	double sigma = pow(10.0, po);

	plot(sigma, gamma, 20);




	return 0;
}*/












