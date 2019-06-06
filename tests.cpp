#include<iostream>
#include<vector>
#include<locale>
#include<cmath>
#include<cassert>
#include<cstring>
#include"CalcHeader.h"
using namespace mcal;
using namespace std;



//tests
//####################################################################################
const double eps = 0.5;
void eps_assert(const double &ans, double &expected)
{
	assert(ans > expected - eps && ans < expected + eps);
}
void vec_assert(const vector<double> &ans, vector<double> &expected)
{
	assert(ans.size() == expected.size());
	for (int i = 0; i < ans.size(); ++i)
		eps_assert(ans[i], expected[i]);
}
void matrix_assert(const vector<vector<double>> &ans, vector<vector<double>> &expected)
{
	assert(ans.size() == expected.size());
	for (int i = 0; i < ans.size(); ++i)
		vec_assert(ans[i], expected[i]);
}
void test_cond()
{
	vector<vector<double>> M;
	double expected;
	M = { {1,2,3},{4,5,6},{7,8,8} };
	expected = 121.531;
	eps_assert(cond(M), expected);
}
void test_gauss()
{
	vector<vector<double>> M;
	double expected;
	vector<double> ans, expected_m;
	//1
	M = { {1,2,3},{4,5,6} };
	expected = 1;
	expected_m = { -1, 2 };
	eps_assert(Gauss(M, ans), expected);
	vec_assert(ans, expected_m);
	//2
	vec_assert(MainElem(M), expected_m);
	//3
	M = { {1,2,0},{2,4,0} };
	expected = 1e9;
	eps_assert(Gauss(M, ans), expected);
}
void test_eigenvalues()
{
	vector<vector<double>> M;
	double expected;
	//1
	M = { {1,2,3},{4,5,6},{7,8,9} };
	expected = 16.1;
	eps_assert(Scal(M), expected);
	eps_assert(Sim(M), expected);
	//2
	M = { {1} };
	expected = 1.0;
	eps_assert(Scal(M), expected);
	eps_assert(Sim(M), expected);
}
void test_det()
{
	vector<vector<double>> M;
	double expected;
	//1
	M = { {1,2,3},{4,5,6},{7,8,9} };
	expected = 0;
	eps_assert(Determinant(M), expected);
	//2
	M = { {1} };
	expected = 1;
	eps_assert(Determinant(M), expected);
}
void test_minor()
{
	vector<vector<double>> M;
	vector<vector<double>> expected;
	//1
	M = { {1,2},{3,4} };
	expected = { {4,-3},{-2,1} };
	matrix_assert(Minor(M), expected);
	//2
	M = { {1} };
	expected = { {0} };
	matrix_assert(Minor(M), expected);
}
void test_trans()
{
	vector<vector<double>> M;
	vector<vector<double>> expected;
	//1
	M = { {1,2},{3,4} };
	expected = { {1,3},{2,4} };
	matrix_assert(Transpose(M), expected);
	//2
	M = { {1} };
	expected = { {1} };
	matrix_assert(Transpose(M), expected);
}
void test_inver()
{
	vector<vector<double>> M;
	vector<vector<double>> expected;
	M = { {1, 2},{3,4} };
	expected = { {-2, 1},{1.5, -0.5} };
	matrix_assert(Invertible(M), expected);
}
void test_norm()
{
	vector<vector<double>> M;
	double expected;
	M = { {1, 2},{3,4} };
	expected = 5.47723;
	eps_assert(Norm(M), expected);
}
void test_multi()
{
	vector<vector<double>> M;
	vector<vector<double>> M1;
	M = { {1, 2}, {3, 4} };
	M1 = { {6, 3}, {2, 8} };
	vector<vector<double>> expected;
	expected = { {10, 19}, {26, 41} };
	matrix_assert(multi(M, M1), expected);
}
void test_multiV()
{
	vector<vector<double>> M;
	vector<double> v;
	M = { {1, 2}, {3, 4} };
	v = { 9, 2 };
	vector<double> expected;
	expected = { 13, 35 };
	vec_assert(multiV(M, v), expected);
}


void tests()
{
	test_cond();
	test_gauss();
	test_eigenvalues();
	test_det();
	test_minor();
	test_trans();
	test_inver();
	test_norm();
	test_multi();
	test_multiV();
}
//####################################################################################

int main()
{
	//cout << "Matrix calc v 2.6.4" << endl;
	tests();
	cout << "tests completed"<<endl;
}
