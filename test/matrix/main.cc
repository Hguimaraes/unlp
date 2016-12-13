#include <iostream>
#include <vector>
#include "matrix.h"

using namespace std;

int main(){
	vector<vector<double> > x;
	vector<double> y;
	y.push_back(0);
	y.push_back(0);
	y.push_back(0);

	x.push_back(y);
	y.clear();
	y.push_back(1);
	y.push_back(1);
	y.push_back(1);

	x.push_back(y);
	y.clear();
	y.push_back(0);
	y.push_back(0);
	y.push_back(0);
	x.push_back(y);

	unlp::Matrix tmp(x);

	

	tmp.printMatrix();
	unlp::Matrix tmp2 = tmp.t();
	tmp2.printMatrix();
	return 0;
}