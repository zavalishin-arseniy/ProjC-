#include<iostream>
#include<vector>
namespace mcal
{

	int Gauss(std::vector<std::vector<double>> a, std::vector<double> &ans);

	std::vector<double> MainElem(std::vector<std::vector<double>> A);

	double Determinant(std::vector<std::vector<double>> A);

	std::vector<std::vector<double>> Minor(std::vector<std::vector<double>> A);

	std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>> A);


	std::vector<std::vector<double>> Invertible(std::vector<std::vector<double>> A);

	double Norm(std::vector<std::vector<double>> A);


	double Norm(std::vector<double> A);


	double cond(std::vector<std::vector<double>> A);


	std::vector<std::vector<double>> multi(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b);


	std::vector<double> multiV(std::vector<std::vector<double>> a, std::vector<double> b);


	std::vector<double> sumV(std::vector<double> a, std::vector<double> b);

	std::vector<double> miV(std::vector<double> a, std::vector<double> b);


	std::vector<double> multV(std::vector<double> a, double b);


	double scal(std::vector<double> a, std::vector<double> b);


	double Scal(std::vector<std::vector<double>> a);

	double Sim(std::vector<std::vector<double>> a);


	void writeM(std::vector<std::vector<double>> a);

}

