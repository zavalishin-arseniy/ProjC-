#include <vector>
#include<iostream>
#include"CalcHeader.h"
namespace mcal
{
	std::vector<double> Gauss(std::vector<std::vector<double>> A)
	{
		int n = (int)A.size();
		std::vector<std::vector<double>> a;
		std::vector<double> tmptmp;
		tmptmp.assign(n + 1, 0);
		a.assign(n, tmptmp);
		std::vector<double> b;
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

	std::vector<double> MainElem(std::vector<std::vector<double>> A)
	{
		int n = (int)A.size();
		std::vector<std::vector<double>> a;
		std::vector<double> tmptmp;
		tmptmp.assign(n + 1, 0);
		a.assign(n, tmptmp);
		std::vector<double> b;
		b.assign(n, 0);
		std::vector<int> c;
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

	double Determinant(std::vector<std::vector<double>> A)
	{
		int n = (int)A.size();
		if (n > 1)
		{
			std::vector<std::vector<double>> a;
			std::vector<double> tmptmp;
			tmptmp.assign(n, 0);
			a.assign(n, tmptmp);
			std::vector<std::vector<double>> c;
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

	std::vector<std::vector<double>> Minor(std::vector<std::vector<double>> A)
	{
		int n = (int)A.size();
		if (n > 1)
		{
			std::vector<std::vector<double>> a;
			std::vector<double> tmptmp;
			tmptmp.assign(n, 0);
			a.assign(n, tmptmp);
			std::vector<std::vector<double>> c;
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
			std::vector<std::vector<double>> a;
			a = { {0} };
			return a;
		}
	}

	std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>> A)
	{
		int n = (int)A.size();
		std::vector<std::vector<double>> a;
		std::vector<double> tmptmp;
		tmptmp.assign(n, 0);
		a.assign(n, tmptmp);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				a[i][j] = A[j][i];
		return a;
	}

	std::vector<std::vector<double>> Invertible(std::vector<std::vector<double>> A)
	{
		int n = (int)A.size();
		std::vector<std::vector<double>> a = Transpose(Minor(A));
		double c = Determinant(A);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / c;
		return a;
	}

	double Norm(std::vector<std::vector<double>> A)
	{
		double a = 0;
		for (int i = 0; i < (int)A[0].size(); i++)
			for (int j = 0; j < (int)A[1].size(); j++)
				a += A[i][j] * A[i][j];
		return sqrt(a);
	}

	double Norm(std::vector<double> A)
	{
		double a = 0;
		for (int i = 0; i < (int)A.size(); i++)
		{
			a += A[i] * A[i];
			//cout << a << ' ' << A[i] << endl;
		}
		return sqrt(a);
	}

	double cond(std::vector<std::vector<double>> A)
	{
		return Norm(A) * Norm(Invertible(A));
	}

	std::vector<std::vector<double>> multi(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b)
	{
		int n = a.size();
		std::vector<std::vector<double>> c;
		std::vector<double> tmptmp;
		tmptmp.assign(n, 0);
		c.assign(n, tmptmp);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < n; k++)
					c[i][j] += a[i][k] * b[k][j];
		return c;
	}

	std::vector<double> multiV(std::vector<std::vector<double>> a, std::vector<double> b)
	{
		int n = a.size();
		std::vector<double> c;
		c.assign(n, 0);
		for (int i = 0; i < n; i++)
			for (int k = 0; k < n; k++)
				c[i] += a[i][k] * b[k];
		return c;
	}

	std::vector<double> sumV(std::vector<double> a, std::vector<double> b)
	{
		int n = a.size();
		std::vector<double> c;
		c.assign(n, 0);
		for (int i = 0; i < n; i++)
			c[i] += a[i] + b[i];
		return c;
	}

	std::vector<double> miV(std::vector<double> a, std::vector<double> b)
	{
		int n = a.size();
		std::vector<double> c;
		c.assign(n, 0);
		for (int i = 0; i < n; i++)
			c[i] = a[i] - b[i];
		return c;
	}

	std::vector<double> multV(std::vector<double> a, double b)
	{
		int n = a.size();
		std::vector<double> c;
		c.assign(n, 0);
		for (int i = 0; i < n; i++)
			c[i] = a[i] * b;
		return c;
	}

	double scal(std::vector<double> a, std::vector<double> b)
	{
		int n = a.size();
		double v = 0;
		for (int i = 0; i < n; i++)
		{
			v += a[i] * b[i];
		}
		return v;
	}

	double Scal(std::vector<std::vector<double>> a)
	{
		int n = a.size();
		std::vector<double> b;
		for (int i = 0; i < n; ++i)
			b.push_back(((int)(a[0][i] * 10)) / 10.0);
		std::vector<double> v = multiV(a, b), w = multiV(a, v);
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

	double Sim(std::vector<std::vector<double>> a)
	{
		int n = a.size();
		std::vector<double> b;
		for (int i = 0; i < n; ++i)
			b.push_back(((int)(a[0][i] * 10)) / 10.0);
		std::vector<double> v = multiV(a, b), w = multiV(a, v);
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

	void writeM(std::vector<std::vector<double>> a)
	{
		int n = a.size();
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				std::cout << a[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}
}