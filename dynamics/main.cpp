// nr3 headers
#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"
#include "nr_headers/fourier.h"
#include "nr_headers/correl.h"
#include "nr_headers/stepper.h"
#include "nr_headers/stepperbs.h"
#include "nr_headers/odeint.h"
#include "nr_headers/stepperdopr853.h"
#include "nr_headers/stepperdopr5.h"
#include "nr_headers/stepperstoerm.h"

// other headers
#include "headers/derivatives.h"
#include "headers/write_matrix.h"
#include "headers/box_muller.h"
#include "headers/read_input.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <math.h>

using namespace::std;

int main()
{

	int N;
	double meanS,meanA;
	double varS,varA;
	int seed;

	int tf;
	int tsave;
	
	read_input(N,meanS,varS,meanA,varA,tf,tsave,seed, "input.txt");
	Ran r(seed);

	vector<double> temp(N,.0);
	vector<vector<double> > w(N,temp);
	VecDoub x(N,0.0);
	double mean =0;
	for(int i=0;i<N;++i){
		x[i] = bm_transform(r);
		for( int j=0;j<i;++j) {
			double a = meanA+varA*bm_transform(r);
			double s = meanS+varS*bm_transform(r);
			w[i][j] = a+s;
			w[j][i] = -1.*a+s;
			mean += 2*s;
		}
		double s= meanS+varS*bm_transform(r);
		w[i][i] = s;
		mean += s;
	}



	double atol = 1.e-9;
	double rtol = atol;
	double h1 = 0.01;
	double hmin = 0.;

	NW nw(w,N);
	Output out(tsave);
	Odeint<StepperDopr5<NW> > ode(x,0,tf,atol,rtol,h1,hmin,out,nw);
	ode.integrate();

	write_matrix(out.xsave,out.count,"t.csv");
	write_matrix(out.ysave,N,out.count,"x.csv");



	return 0;
}


