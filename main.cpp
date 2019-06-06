#include<iostream>
#include<vector>
#include<locale>
#include<cmath>
#include<cassert>
#include<cstring>
#include"CalcHeader.h"
using namespace std;
using namespace mcal;

bool exitflag = false;
int user;
int n;
vector<double> str;
vector<vector<double>> matr, matr1, matrans;
vector<double> vec, vec1;
double doub;

int main()
{
	cout << "Matrix calc v 3.6.4" << endl;
	while (!exitflag)
	{
		cout << "__________________" << endl;
		cout << "1. Condition number" << endl;
		cout << "2. Gauss method" << endl;
		cout << "3. Main element method" << endl;

		cout << "4. Power method" << endl;
		cout << "5. Scalar product method" << endl;

		cout << "6. Determinant" << endl;
		cout << "7. Minor" << endl;

		cout << "8. Transposition" << endl;
		cout << "9. Invertible matrix" << endl;

		cout << "10. Matrix norm" << endl;

		cout << "11. Matrix multiplication" << endl;
		cout << "12. Matrix and vector multilication" << endl;

		cout << "13. Exit" << endl;
		cout << "__________________" << endl;
		int x; matr.clear();
		cin >> user;
		switch (user)
		{
		case 1:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			doub = cond(matr);
			cout << endl;
			cout << doub;
			cout << endl;
			break;

		case 2:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n + 1, 0);
				for (int j = 0; j < n + 1; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			cout << Gauss(matr, vec);
			cout << endl;
			for (auto i : vec)
				cout << i << ' ';
			cout << endl;
			break;

		case 3:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n + 1, 0);
				for (int j = 0; j < n + 1; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			vec = MainElem(matr);
			cout << endl;
			for (auto i : vec)
				cout << i << ' ';
			cout << endl;
			break;

		case 4:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			doub = Sim(matr);
			cout << endl;
			cout << doub;
			cout << endl;
			break;
		case 5:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			doub = Scal(matr);
			cout << endl;
			cout << doub;
			cout << endl;
			break;

		case 6:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			doub = Determinant(matr);
			cout << endl;
			cout << doub;
			cout << endl;
			break;

		case 7:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			matrans = Minor(matr);
			cout << endl;
			writeM(matrans);
			break;

		case 8:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			matrans = Transpose(matr);
			cout << endl;
			writeM(matrans);
			break;

		case 9:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			matrans = Invertible(matr);
			cout << endl;
			writeM(matrans);
			break;

		case 10:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			doub = Norm(matr);
			cout << endl;
			cout << doub;
			cout << endl;
			break;

		case 11:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matrans.push_back(str);
			}

			matrans = multi(matr, matrans);
			cout << endl;
			writeM(matrans);
			break;

		case 12:
			cin >> n;
			for (int i = 0; i < n; ++i)
			{
				str.assign(n, 0);
				for (int j = 0; j < n; ++j)
					cin >> str[j];
				matr.push_back(str);
			}
			str.assign(n, 0);
			for (int i = 0; i < n; ++i)
				cin >> str[i];
			vec = multiV(matr, str);
			cout << endl;
			for (auto i : vec)
				cout << i << ' ';
			cout << endl;
			break;

		case 13:
			exitflag = true;
			break;

		default:
			cout << "enter valid data" << endl;
			break;
		}
	}
}
