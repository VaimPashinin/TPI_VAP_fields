#include <iostream>
#include <vector>
#include <string>
#include <particle.h>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>

using namespace std;

double _bSpline(double t, double step) // Базовый кубический сплайн
{
	double result = 0;
	if (t < step) result = 4/6 - pow(abs(t), 2.0) + 0.5 * pow(abs(t), 3.0);
	else if (t < step * 2) result = pow(2 - abs(t), 3.0) / 6;
	return result;
}

void GetInput(string filename, vector<particle>& partilceList)
{
	ifstream input;
	input.open(filename);
	string line;
	if (input.is_open())
	{
		while (getline(input, line))
		{
			cout << line << endl;
			// TODO сделать адекватное считывание из файла
		}
	}
	input.close();
}

void FlowFit_fg(vector<particle>& partilceList, vector<vector<Eigen::VectorXd>>& velocityWeights, vector<vector<Eigen::VectorXd>>& accelerationWeights)
{
	// Первое поколение FlowFit. Вычисляет поля скоростей v и ускорений a
	// В точке пространства в момент времени t вектор выражается в виде линейной комбинации взвешенных свёрток базовых сплайнов 
	// Веса сплайнов вычисляются аналогично TrackFit: для ограниченного объёма пространства строится равномерная сетка с шагом h
	// Условия на h таковы, что отношение количесва частиц к произведению числа узлов сетки по координатам должно лежать в пределах [0.05, 0.33]
	// В узлах сетки вычисляются скорости и ускорения точек с использованием вычисленных траекторий точек
	// В качестве скорости в узле выбирается скорость ближайшей к узлу точки
	// Для удобства вычислений набор узлов c_ijk представляется в виде вектора c_l, l = i + I*j + I*J*k
	// Поиск вектора c_l осуществляется с помощью метода наименьших квадратов: F(c(t)) = sum_(l = 0)^(I+J+K) |v_l(t) - u_l(t)|^2
	// где u_l(t) - скорости в узлах, v_l(t) - взвешенные сплайны, F(c(t)) - функция потерь
	// Вычисления коэфицентов v и a абсолютно аналогичны, поиск осуществляется для каждого кадра

	double min_x = 100; double min_y = 100; double min_z = 100;
	double max_x = -100; double max_y = -100; double max_z = -100;
	for (int i = 0; i < partilceList.length(); i++)
	{
		partilceList[i].TraceFit();
		min_x = min_x > partilceList[i].GetMin(1) ? partilceList[i].GetMin(1) : min_x; // Необходимые границы области для вычисления значений I, J, K и h
		min_y = min_y > partilceList[i].GetMin(2) ? partilceList[i].GetMin(2) : min_y;
		min_z = min_z > partilceList[i].GetMin(3) ? partilceList[i].GetMin(3) : min_z;
		max_x = max_x > partilceList[i].GetMax(1) ? partilceList[i].GetMax(1) : max_x;
		max_y = max_y > partilceList[i].GetMax(2) ? partilceList[i].GetMax(2) : max_y;
		max_z = max_z > partilceList[i].GetMax(3) ? partilceList[i].GetMax(3) : max_z;
	}

	double h; // TODO вычисление шага и IJK
	int I; int J; int K;

	for (int coord = 1; coord < 4; coord++)
	{
		for (double t = 0.0; t < 11,9; t += 0.1)
		{
			Eigen::MatrixXd A((I+1)*(J+1)*(K+1), (I+1)*(J+1)*(K+1));
            Eigen::VectorXd bv((I+1)*(J+1)*(K+1));
            Eigen::VectorXd ba((I+1)*(J+1)*(K+1));
            for (int i = 0; i < (I+1)*(J+1)*(K+1); i++)
            {
                for (int j = 0; j < (I+1)*(J+1)*(K+1); j++) // TODO исправить кашу в индексах
                {
                    double aij = 0;
                    for (int l = 0; k < (I+1)*(J+1)*(K+1); l++)
                    {
						double x = 
                        aij += _bSpline(GetPosition(1,0.1*l) - GetPosition(1,0.1*j), step) * _bSpline(GetPosition(1,0.1*l) - GetPosition(1,0.1*i));
                    }
                    double biv = 0;
					double bia = 0;
                    for (int l = 0; k < (I+1)*(J+1)*(K+1); l++)
                    {
                        biv += _bSpline(GetPosition(1, 0.1*l) - GetPosition(1, 0.1*i)) * _bSpline(GetPosition(1,0.1*l) - GetPosition(1,0.1*i)) * _bSpline(GetPosition(1,0.1*l) - GetPosition(1,0.1*i)) * GetVelocity(coord, GetPosition(1,0.1*l));
						bia += _bSpline_bSpline(GetPosition(1,0.1*l) - GetPosition(1,0.1*i))*_bSpline(GetPosition(1,0.1*l) - GetPosition(1,0.1*i))*_bSpline(GetPosition(1,0.1*l) - GetPosition(1,0.1*i)) * GetAcceleration(coord, GetPosition(1,0.1*l)); 
                    }
                    A << aij;
                    bv << biv;
					ba << bia;
                }
            }

            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
            Eigen::VectorXd x_v = qr.solve(bv);
			Eigen:: VectorXd x_a = qr.solve(ba);
		}
	}
}

void FlowFit_sg()
{
	// TODO второй шаг FlowFit, поиск поля давлений
	// Примечание: воспользоваться уравнением на странице 9
	// laplas(p) + grad(u * div(u)) = 0
}

int main()
{
	vector<particle> partilceList();
	GetInput("Perm_Tracs_Maks.csv", partilceList);
	FlowFit_fg();
	FlowFit_sg();

	//cin.get();
	return 0;
}

