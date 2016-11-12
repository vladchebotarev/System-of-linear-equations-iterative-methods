// Lab_7.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <Windows.h>
#include <math.h>

using namespace std;

const int n = 4;
const int il_petli = 100;
const double eps = 1e-10;

void printMatrix(double **mac)
{

	for (int i = 0; i<n; i++)
	{
		cout << "| ";
		for (int j = 0; j<n; j++)
		{
			cout.width(4);
			cout << mac[i][j];
		}
		cout << " |" << endl;
	}
}

void printVector(double *wek)
{
	for (int i = 0; i<n; i++)
	{

		cout << "| ";
		cout.width(4);
		cout << wek[i];
		cout << "| ";
		cout << endl;
	}
}

void insertData(double **A, double *b, double *x)
{
	A[0][0] = 100.0; A[0][1] = -1.0; A[0][2] = 2.0; A[0][3] = -3.0;
	A[1][0] = 1.0; A[1][1] = 200.0; A[1][2] = -4.0; A[1][3] = 5.0;
	A[2][0] = -2.0; A[2][1] = 4.0; A[2][2] = 300.0; A[2][3] = -6.0;
	A[3][0] = 3.0; A[3][1] = -5.0; A[3][2] = 6.0; A[3][3] = 400.0;

	b[0] = 116.0; b[1] = -226.0; b[2] = 912.0; b[3] = -1174.0;

	x[0] = 2.0; x[1] = 2.0; x[2] = 2.0; x[3] = 2.0;
}

double est(double *x, double *x_nowe)
{
	x[0] = fabs(x[0] - x_nowe[0]);
	x[1] = fabs(x[1] - x_nowe[1]);
	x[2] = fabs(x[2] - x_nowe[2]);
	x[3] = fabs(x[3] - x_nowe[3]);

	double max = x[0];
	if (x[1] > max) max = x[1];
	if (x[2] > max) max = x[2];
	if (x[3] > max) max = x[3];

	return max;
}

double residuum(double **A, double *b, double *x_nowe)
{
	double Ax[n];

	Ax[0] = fabs((A[0][0] * x_nowe[0] + A[0][1] * x_nowe[1] + A[0][2] * x_nowe[2] + A[0][3] * x_nowe[3]) - b[0]);
	Ax[1] = fabs((A[1][0] * x_nowe[0] + A[1][1] * x_nowe[1] + A[1][2] * x_nowe[2] + A[1][3] * x_nowe[3]) - b[1]);
	Ax[2] = fabs((A[2][0] * x_nowe[0] + A[2][1] * x_nowe[1] + A[2][2] * x_nowe[2] + A[2][3] * x_nowe[3]) - b[2]);
	Ax[3] = fabs((A[3][0] * x_nowe[0] + A[3][1] * x_nowe[1] + A[3][2] * x_nowe[2] + A[3][3] * x_nowe[3]) - b[3]);

	double max = Ax[0];
	if (Ax[1] > max) max = Ax[1];
	if (Ax[2] > max) max = Ax[2];
	if (Ax[3] > max) max = Ax[3];

	return max;
}

void metoda_Jacobiego(double **A, double *b, double *x)
{
	double *x_nowe = new double[n]; //nowe przyblizenia
	double suma = 0.0;
	double EST = 0.0, RESIDUUM = 0.0;

	cout << endl << endl << "Metoda Jacobiego" << endl;
	cout << "  n |            x1 |            x2 |            x3 |            x4 |            EST |     RESIDIUM |" << endl;
	cout << "-----------------------------------------------------------------------------------------------------" << endl;

	for (int iter = 0; iter < il_petli; iter++)
	{
		for (int i = 0; i < n; i++)
		{
			suma = 0.0;
			for (int j = 0; j < n; j++)
			if (j != i)
				suma += A[i][j] * x[j];

			x_nowe[i] = (1.0 / A[i][i]) * (b[i] - suma);

		}

		EST = est(x, x_nowe);
		RESIDUUM = residuum(A, b, x_nowe);

		for (int i = 0; i < n; i++)
			x[i] = x_nowe[i];

		cout.width(4);
		cout << iter << "|";
		cout.width(15);
		cout << x_nowe[0] << "|";
		cout.width(15);
		cout << x_nowe[1] << "|";
		cout.width(15);
		cout << x_nowe[2] << "|";
		cout.width(15);
		cout << x_nowe[3] << "|";
		cout.width(15);
		cout << EST << "|";
		cout.width(15);
		cout << RESIDUUM << "|" << endl;

		if (EST < eps && RESIDUUM < eps)
			break;
	}
	cout << "-----------------------------------------------------------------------------------------------------" << endl;
}


void metoda_Gaussa_Seidela(double **A, double *b, double *x)
{
	double *x_poprz = new double[n]; //stare wart
	double suma = 0.0;
	double EST = 0.0, RESIDUUM = 0.0;

	cout << endl << endl << "\t Metoda Gaussa_Seidela" << endl;
	cout << "  n |            x1 |            x2 |            x3 |            x4 |            EST |     RESIDIUM |" << endl;
	cout << "-----------------------------------------------------------------------------------------------------" << endl;

	for (int iter = 0; iter < il_petli; iter++)
	{
		for (int i = 0; i < n; i++)
		{
			suma = 0.0;
			for (int j = 0; j < n; j++){
				if (j != i)
					suma += A[i][j] * x[j];
			}
			x_poprz[i] = x[i];
			x[i] = (1.0 / A[i][i]) * (b[i] - suma);
		}

		EST = est(x_poprz, x);
		RESIDUUM = residuum(A, b, x);


		cout.width(4);
		cout << iter << "|";
		cout.width(15);
		cout << x[0] << "|";
		cout.width(15);
		cout << x[1] << "|";
		cout.width(15);
		cout << x[2] << "|";
		cout.width(15);
		cout << x[3] << "|";
		cout.width(15);
		cout << EST << "|";
		cout.width(15);
		cout << RESIDUUM << "|" << endl;

		if (EST < eps&&RESIDUUM < eps)
			break;


	}
	cout << "-----------------------------------------------------------------------------------------------------" << endl;
}

void metoda_SOR(double **A, double *b, double *x)
{
	double *x_nowe = new double[n]; //nowe przyblizenia
	double *x_poprz = new double[n]; //stare wartosc
	double suma = 0.0, omega = 1.0;
	double EST = 0.0, RESIDUUM = 0.0;

	cout << endl << endl << "\t Metoda SOR" << endl;

	cout << "  n |            x1 |            x2 |            x3 |            x4 |            EST |     RESIDIUM |" << endl;
	cout << "-----------------------------------------------------------------------------------------------------" << endl;

	for (int iter = 0; iter < il_petli; iter++)
	{
		for (int i = 0; i < n; i++)
		{
			suma = 0.0;
			for (int j = 0; j < n; j++)
			if (j != i)
				suma += A[i][j] * x[j];

			x_poprz[i] = x[i];
			x_nowe[i] = (1.0 - omega) * x[i] + (omega / A[i][i]) * (b[i] - suma);
			x[i] = x_nowe[i];
		}

		EST = est(x_poprz, x_nowe);
		RESIDUUM = residuum(A, b, x_nowe);


		cout.width(4);
		cout << iter << "|";
		cout.width(15);
		cout << x[0] << "|";
		cout.width(15);
		cout << x[1] << "|";
		cout.width(15);
		cout << x[2] << "|";
		cout.width(15);
		cout << x[3] << "|";
		cout.width(15);
		cout << EST << "|";
		cout.width(15);
		cout << RESIDUUM << "|" << endl;

		if (EST < eps && RESIDUUM < eps)
			break;


	}
	cout << "-----------------------------------------------------------------------------------------------------" << endl;
}

void rozwiazanie()
{
	double **A, *b, *x;

	A = new double*[n];

	for (int i = 0; i < n; i++)
		A[i] = new double[n];

	b = new double[n];

	x = new double[n];

	insertData(A, b, x);

	cout << "Macierz A:" << endl;
	printMatrix(A);
	cout << endl;
	cout << " Wektor b:" << endl;
	printVector(b);
	cout << " Wektor x:" << endl;
	printVector(x);

	metoda_Jacobiego(A, b, x);
	system("pause");
	insertData(A, b, x);
	metoda_Gaussa_Seidela(A, b, x);
	system("pause");
	insertData(A, b, x);
	metoda_SOR(A, b, x);

	for (int i = 0; i < n; i++)
	{
		delete[] A[i];
	}

	delete[] A;
	delete[] b;
	delete[] x;

}

int _tmain(int argc, _TCHAR* argv[])
{
	rozwiazanie();
	system("pause");
	return 0;
}


