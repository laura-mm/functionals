// lotka-volterra dynamics with eigen array class, RK4



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

ArrayXd functiom(MatrixXd m) // piecewise
{
	ArrayXd v = m.array();
	for (int i = 0; i < v.size(); i++)
	{
		if (v(i) >= a) v(i) = a;
		else if (v(i) <= -a) v(i) = -a;
	}
	return v;
}


/*ArrayXd functiom(MatrixXd m) // holling
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

			// if (x.maxCoeff() > pow(10.0, 20.0) || y.maxCoeff() > pow(10.0, 20.0)) {diverge = true; break;} // linear only

			//trajx.push_back(x); trajy.push_back(y); //

			if (t >= (0.99*T)-2) {trajx.push_back(x); trajy.push_back(y);} //
		}
		check();
	}
	ArrayXd measures()
	{
		ArrayXd m = ArrayXd::Zero(8); // top, middle, bottom, M, q, diversity, d, h the last 3 are relative things
		
		ArrayXd sum1 = ArrayXd::Zero(N);
		ArrayXd sum2 = ArrayXd::Zero(N);

		double Msq = 0.0;

		for (int t = 0; t < 0.01*T; t++)// this is not exactly last 1% is off by two, but this is negligable
		{
			sum1 += trajx[t];
			sum2 += trajx[t]*trajx[t];

			double M = trajx[t].mean();
			double q = (trajx[t]*trajx[t]).mean();
			for (int i = 0; i < N; i++)
			{
				if (abs(trajx[t](i) - 1.0 - a) < 0.0001) m(0) ++;
				else if (abs(trajx[t](i) - 1.0 + a) < 0.0001 || abs(trajx[t](i)) < 0.0001) m(2) ++;
				else m(1) ++;
			}

			m(3) += M;
			m(4) += q;
			m(5) += M*M/q;
			//m(6) += (((trajx[t] - trajy[t])*(trajx[t] - trajy[t])).mean())/(M*M); // comes from previous d
			m(6) += (((trajx[t] - trajy[t])*(trajx[t] - trajy[t])).mean());
			Msq += M*M;
		}
		

		sum1 /= (0.01*(double)T);
		sum2 /= (0.01*(double)T);

		m(9) = m(8)/m(4);

		m(0) /= 0.01*(double)N*(double)T;
		m(1) /= 0.01*(double)N*(double)T;
		m(2) /= 0.01*(double)N*(double)T;
		m(3) /= 0.01*(double)T;
		m(4) /= 0.01*(double)T;
		m(5) /= 0.01*(double)T;
		//m(6) /= 0.01*(double)T; // comes from previous d
		m(6) /= Msq;

		sum2 -= (sum1*sum1);
		//for (int i = 0; i < N; i ++) if (abs(sum1(i)) < pow(10.0, -10.0)) sum2(i) = 0.0; // what is this for??
		//m(7) = (sum2/(sum1*sum1)).mean();
		m(7) = sum2.mean()/(sum1*sum1).mean();
		//m(8) = (sum2/(sum1*sum1)).mean();

		return m;
	}

};

void fiveplot(int grid, int runs)
{
	string filename = "measpV2_" + to_string((int)(10*a)) + "_" + to_string((int)(10*mu)) + ".txt"; //
	ofstream file; file.open(filename);
	for (int g = 0; g <= 4; g++)
	{
		//int g = 2; //
		double gamma = (0.5*(double)g) - 1.0; 
		for (int s = 0; s <= grid; s++)
		{
			//int s = 15; //
			double po = (2.0*s/(double)grid) - 1.0;
			double sigma = pow(10.0, po);

			ArrayXd meas = ArrayXd::Zero(8);
		
			for (int r = 0; r < runs; r++)
			{
				cout << "g " << g << ", s " << s << ", r " << r << endl;
				simulation sim(sigma, gamma);
				sim.run();
				meas += sim.measures();
			}
			meas /= (double)runs;
			file << meas.format(c);
			if (s != grid) file << ",";
		}
		if (g != 4) file << ",";
	}
}
		
				


/*
void histogram(int runs)
{
	string filename = "hist2_" + to_string((int)(10*a)) + "_" + to_string((int)(10*mu)) + ".txt"; //
	ofstream file; file.open(filename);
	for (int gi = 0; gi <= 2; gi++)
	{
		double gamma = (double)gi - 1.0;
		for (int si = 0; si <= 2; si++)
		{
			double sigma = pow(10.0, ((double)si*0.5) - 1.0);
		
			for (int r = 0; r < runs; r++)
			{
				cout << gi << ", " << si << ", " << r << endl;
				simulation sim(sigma, gamma);
				sim.run();
				file << sim.x.format(c);
				if (r != runs-1) file << ",";
			}
			if (si != 2) file << ",";
		}
		if (gi != 2) file << ",";
	}
	file.close();
}


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
		for (int t = 0; t < T; t++)
		{
			datax << sim.trajx[t](i);
			datay << sim.trajy[t](i);
			if (i != plots-1 || t != T-1) {datax << ","; datay << ",";}
		}
	}
}

*/

			



int main()	
{
	mu = 0.0;
	//double sigma = 0.5; //pow(10.0, 3.0);
	a = 0.5;
	N = 200;
	T = 200000;
	dt = 0.001;
	int grid = 20; // for fiveplots
	int runs = 20;
	//double gamma = 0.0;
	int plots = 20;

	// function h in fiveplots()
	// function2 in histogram()
	// traj matrix 1%
	// check fixed and unique at end
	// functions containing traj

	/*simulation sim(pow(10.0, 1.0), 0.0);
	sim.run();
	cout << sim.fixed << endl;
	cout << sim.unique << endl;*/









	
	//histogram(runs);

	fiveplot(grid, runs);

	//plot(pow(10.0, 1.0), 0.0, 20);



	return 0;
}




