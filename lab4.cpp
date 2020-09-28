#include <iostream>
#include <cmath>

using namespace std;

void fill(double **a, double *b, int n){ // Заполнение матрицы коэффициентов
    a[0][0] = 1;
    b[0] = 1;
    for(int i = 1; i < n; i++)
        a[0][i] = 0;
    
    for(int i = 1; i < n-1; i++) {
        for(int j = 0; j < n; j++){
            if(i == j) a[i][j] = -2;
            else if(j == i-1 || j == i+1) a[i][j] = 1;
            else a[i][j] = 0;
        }
        b[i] = 2./pow(i+1, 2);
    }
    a[n-1][0] = 1;
    a[n-1][n-1] = 1;
    b[n-1] = -19./3.;
    for(int i = 1; i < n-1; i++)
        a[n-1][i] = 2;
}

void gauss(double **a, double *b, double *x, int m){ // Метод Гаусса
	for(int k = 1; k < m; k++){
		for(int j = k; j < m; j++){
			double h = a[j][k-1]/a[k-1][k-1];
			for(int i = 0; i < m; i++)
				a[j][i] -= h*a[k-1][i];
			b[j] -= h*b[k-1];
		}
		for(int i = m-1; i >= 0; i--){
			x[i] = b[i]/a[i][i];
			for(int c = m-1; c > i; c--)
				x[i] -= a[i][c]*x[c]/a[i][i];
		}
	}
}


void g_z(double **a, double *b, double *x, int m, double eps){ // Метод Гаусса-Зейделя
	double p[m];
	double delta;
	int iter = 0;
	do{
	    for(int i = 0; i < m; i++)
	        p[i] = x[i];

	    for(int i = 0; i < m; i++){
	        double var = 0.;
	        for(int j = 0; j < i; j++)
	            var += a[i][j]*x[j];

	        for(int j = i+1; j < m; j++)
	            var += a[i][j]*p[j];

	        x[i] = (b[i] - var)/a[i][i];
	    }
	    delta = 0;
	    for(int i = 0; i < m; i++)
	        delta += (x[i] - p[i])*(x[i] - p[i]);
		iter++;
	} while(sqrt(delta) >= eps);
	cout << "iterations: " << iter << endl;
}

void invert(double **a, int n){ // Поиск обратной матрицы методом Гаусса-Жордана
	double **b = new double*[n]; // Единичная матрица
	for(int i = 0; i < n; i++){
		b[i] = new double[n];
		b[i][i] = 1.;
	}
	double **ab = new double*[n];
	for(int i = 0; i < n; i++) ab[i] = new double[2*n];

	for(int i = 0; i < n; i++){
	    for(int j = 0; j < n; j++){
	        ab[i][j] = a[i][j];
	        ab[i][j + n] = b[i][j];
	    }
	}
	// Зануление нижнего левого угла
	for(int k = 0; k < n; k++){
	    for(int i = 0; i < 2*n; i++)
	        ab[k][i] = ab[k][i]/a[k][k];

		for(int i = k + 1; i < n; i++){
	        double K = ab[i][k]/ab[k][k];
			for(int j = 0; j < 2*n; j++)
	            ab[i][j] -= ab[k][j]*K;
	    }
	    for(int i = 0; i < n; i++)
	        for(int j = 0; j < n; j++)
	            a[i][j] = ab[i][j];
	}
	// Зануление верхнего правого угла
	for(int k = n - 1; k > -1; k--){
	    for(int i = 2*n - 1; i > -1; i--)
	        ab[k][i] = ab[k][i]/a[k][k];

	    for(int i = k - 1; i > -1; i--){
	        double K = ab[i][k]/ab[k][k];
	        for(int j = 2*n - 1; j > -1; j--)
	            ab[i][j] -= ab[k][j]*K;
	    }
	}

	for(int i = 0; i < n; i++)
	    for(int j = 0; j < n; j++)
	        a[i][j] = ab[i][j + n];
	delete[] b;
	delete[] ab;
}

double norm(double* x, int n) {
    double s = 0;
    for(int i = 0; i < n; i++)
        s += x[i]*x[i];
    return sqrt(s);
}

int main(){
	int n = 20;
	double **a = new double*[n];
	for(int i = 0; i < n; i++) a[i] = new double[n];
	double b[n];
	double x[n];
	for(int i = 0; i < n; i++) x[i] = 0;

	fill(a, b, n);
	cout << "Metod Gaussa-Zeidelya:\n";

	g_z(a, b, x, n, 1e-4);
	// Вычисление вектора невязки
	double r[n] /*= {}*/;
	for(int i = 0; i < n; i++){
		double s = 0.;
		for(int j = 0; j < n; j++)
			s += a[i][j]*x[j];
		r[i] = b[i] - s;
	}
	cout << "Norma vektora nevyazki: " << norm(r, n) << endl;

	gauss(a, b, x, n);

	cout << "\nMetod Gaussa\n";

	fill(a, b, n);
	// Вычисление вектора невязки
	for(int i = 0; i < n; i++){
		double s = 0.;
		for(int j = 0; j < n; j++)
			s += a[i][j]*x[j];
		r[i]= b[i] - s;
	}
	cout << "Norma vektora nevyazki: " << norm(r, n) << endl;

	double mu = 1.; // Число обусловленности
	double ax_norm = 0.;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++)
			ax_norm += pow(a[i][j]*x[j], 2);
	}
	mu *= sqrt(ax_norm)/pow(norm(x, n), 2);

	double inv_ax_norm = 0.;
	invert(a, n);
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			inv_ax_norm += pow(a[i][j]*x[j], 2);
	mu *= sqrt(inv_ax_norm);

	cout << "\nChislo obuslovlennosti: " << mu << endl;

	delete[] a;
	return 0;
}
