#include <iostream>
#include <vector>
#include "matrix.h"

using namespace std;

int main(){
	unlp::Matrix test_a(2,2,0);
	test_a.set(1,1,1);
	test_a.set(2,2,1);
	test_a.printMatrix();

	vector<double> x;
	x.push_back(0.3);
	x.push_back(0.5);
	unlp::Matrix test_b(x);
	test_b.t().printMatrix();

	unlp::Matrix test_c = test_b * test_b.t();
	test_c.printMatrix();

	return 0;
}