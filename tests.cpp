#include<iostream>
#include<vector>
#include<locale>
#include<cmath>
#include<cassert>
#include<cstring>
using namespace std;

vector<double> Gauss(vector<vector<double>> A)
{
	int n = (int)A.size();
	vector<vector<double>> a;
	vector<double> tmptmp;
	tmptmp.assign(n + 1, 0);
	a.assign(n, tmptmp);
	vector<double> b;
	b.assign(n, 0);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n + 1; j++)
			a[i][j] = A[i][j];

	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
			for (int y = i; y < n + 1; y++)
				a[j][y] = a[j][y] - A[j][i] * a[i][y] / A[i][i];

	for (int i = n - 1; i >= 0; i--)
		for (int j = i - 1; j >= 0; j--)
			for (int y = i + 1; y >= 0; y--)
				a[j][y] = a[j][y] - a[j][i] * a[i][y] / a[i][i];

	for (int i = 0; i < n; i++)
	{
		a[i][n] = a[i][n] / a[i][i];
		a[i][i] = a[i][i] / a[i][i];
		b[i] = a[i][n];
	}
	/*for (auto i : a)
	{
		for (auto j : i)
			cout << j << ' ';
		cout << endl;
	}*/

	return b;
}

vector<double> MainElem(vector<vector<double>> A)
{
	int n = (int)A.size();
	vector<vector<double>> a;
	vector<double> tmptmp;
	tmptmp.assign(n + 1, 0);
	a.assign(n, tmptmp);
	vector<double> b;
	b.assign(n, 0);
	vector<int> c;
	c.assign(n, 0);
	double tmp = 0;
	for (int i = 0; i < n; i++)
	{
		c[i] = i;
		for (int j = 0; j < n + 1; j++)
			a[i][j] = A[i][j];
	}

	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
			if (abs(a[i][j]) >= abs(a[i][i]))
				for (int y = 0; y < n + 1; y++)
				{
					tmp = c[i];
					c[i] = c[j];
					c[j] = (int)tmp;
					tmp = a[i][y];
					a[i][y] = a[j][y];
					a[j][y] = tmp;
				}
	b = Gauss(a);
	for (int i = 0; i < n; i++)
	{
		tmp = b[i];
		b[i] = b[c[i]];
		b[c[i]] = (int)tmp;
	}
	return b;
}

double Determinant(vector<vector<double>> A)
{
	int n = (int)A.size();
	if (n > 1)
	{
		vector<vector<double>> a;
		vector<double> tmptmp;
		tmptmp.assign(n, 0);
		a.assign(n, tmptmp);
		vector<vector<double>> c;
		tmptmp.assign(n - 1, 0);
		c.assign(n - 1, tmptmp);
		double b = 0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				a[i][j] = A[i][j];


		for (int i = 0; i < n; i++)
		{
			for (int y = 0; y < n - 1; y++)
				for (int j = 0; j < n - 1; j++)
					if (j < i)
						c[y][j] = A[y + 1][j];
					else
						c[y][j] = A[y + 1][j + 1];


			if (i % 2 == 0) b += a[0][i] * Determinant(c);
			else b -= a[0][i] * Determinant(c);
		}
		return b;
	}
	else return A[0][0];
}

vector<vector<double>> Minor(vector<vector<double>> A)
{
	int n = (int)A.size();
	if (n > 1)
	{
		vector<vector<double>> a;
		vector<double> tmptmp;
		tmptmp.assign(n, 0);
		a.assign(n, tmptmp);
		vector<vector<double>> c;
		tmptmp.assign(n - 1, 0);
		c.assign(n - 1, tmptmp);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int y = 0; y < n - 1; y++)
				{
					for (int u = 0; u < n - 1; u++)
					{
						if (y < i)
							if (u < j) c[y][u] = A[y][u];
							else c[y][u] = A[y][u + 1];
						else
							if (u < j)  c[y][u] = A[y + 1][u];
							else c[y][u] = A[y + 1][u + 1];
					}
				}
				if ((i + j) % 2 == 0) a[i][j] = Determinant(c);
				else a[i][j] = -Determinant(c);
			}
		}
		return a;
	}
	else
	{
		vector<vector<double>> a;
		a = { {0} };
		return a;
	}
}

vector<vector<double>> Transpose(vector<vector<double>> A)
{
	int n = (int)A.size();
	vector<vector<double>> a;
	vector<double> tmptmp;
	tmptmp.assign(n, 0);
	a.assign(n, tmptmp);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a[i][j] = A[j][i];
	return a;
}

vector<vector<double>> Invertible(vector<vector<double>> A)
{
	int n = (int)A.size();
	vector<vector<double>> a = Transpose(Minor(A));
	double c = Determinant(A);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a[i][j] = a[i][j] / c;
	return a;
}

double Norm(vector<vector<double>> A)
{
	double a = 0;
	for (int i = 0; i < (int)A[0].size(); i++)
		for (int j = 0; j < (int)A[1].size(); j++)
			a += A[i][j] * A[i][j];
	return sqrt(a);
}

double Norm(vector<double> A)
{
	double a = 0;
	for (int i = 0; i < (int)A.size(); i++)
	{
		a += A[i] * A[i];
		//cout << a << ' ' << A[i] << endl;
	}
	return sqrt(a);
}

double cond(vector<vector<double>> A)
{
	return Norm(A) * Norm(Invertible(A));
}

vector<vector<double>> multi(vector<vector<double>> a, vector<vector<double>> b)
{
	int n = a.size();
	vector<vector<double>> c;
	vector<double> tmptmp;
	tmptmp.assign(n, 0);
	c.assign(n, tmptmp);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				c[i][j] += a[i][k] * b[k][j];
	return c;
}

vector<double> multiV(vector<vector<double>> a, vector<double> b)
{
	int n = a.size();
	vector<double> c;
	c.assign(n, 0);
	for (int i = 0; i < n; i++)
		for (int k = 0; k < n; k++)
			c[i] += a[i][k] * b[k];
	return c;
}

vector<double> sumV(vector<double> a, vector<double> b)
{
	int n = a.size();
	vector<double> c;
	c.assign(n, 0);
	for (int i = 0; i < n; i++)
		c[i] += a[i] + b[i];
	return c;
}

vector<double> miV(vector<double> a, vector<double> b)
{
	int n = a.size();
	vector<double> c;
	c.assign(n, 0);
	for (int i = 0; i < n; i++)
		c[i] = a[i] - b[i];
	return c;
}

vector<double> multV(vector<double> a, double b)
{
	int n = a.size();
	vector<double> c;
	c.assign(n, 0);
	for (int i = 0; i < n; i++)
		c[i] = a[i] * b;
	return c;
}

double scal(vector<double> a, vector<double> b)
{
	int n = a.size();
	double v = 0;
	for (int i = 0; i < n; i++)
	{
		v += a[i] * b[i];
	}
	return v;
}

double Scal(vector<vector<double>> a)
{
	int n = a.size();
	vector<double> b;
	for (int i = 0; i < n; ++i)
		b.push_back(((int)(a[0][i] * 10)) / 10.0);
	vector<double> v = multiV(a, b), w = multiV(a, v);
	double c = scal(w, v) / scal(v, v), d = 0;
	for (int i = 0; i < 10000; i++)
	{
		v = w;
		w = multiV(a, v);
		d = c;
		c = scal(w, v) / scal(v, v);
		if (abs(c - d) < 0.001)
			return c;
	}
	return c;
}

double Sim(vector<vector<double>> a)
{
	int n = a.size();
	vector<double> b;
	for (int i = 0; i < n; ++i)
		b.push_back(((int)(a[0][i] * 10)) / 10.0);
	vector<double> v = multiV(a, b), w = multiV(a, v);
	double c = w[0] / v[0], d = 0;
	for (int i = 0; i < 10000; i++)
	{
		v = w;
		w = multiV(a, v);
		d = c;
		c = w[0] / v[0];
		//cout << i << ") " << c << endl;
		if (abs(c - d) < 0.001)
			return c;
	}
	return c;
}

void writeM(vector<vector<double>> a)
{
	int n = a.size();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << endl;
	}
}


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
	vector<double> expected;
	M = { {1,2,3},{4,5,6} };
	expected = { -1, 2 };
	vec_assert(Gauss(M), expected);
	vec_assert(MainElem(M), expected);
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
