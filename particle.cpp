#include <vector>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#define PI = 3.1415
#define FCROSS = 0.1
using namespace std;

class particle //Частица. Хранит идентификатор и набор положений в пространстве-времени
{
    private:
    int _id; // ID
    vector<double>* _t; // Время, ниже пространственные координаты по x, y и z
    vector<double>* _x; 
    vector<double>* _y;
    vector<double>* _z;
    double min_x;
    double min_y;
    double min_z;
    double max_x;
    double max_y;
    double max_z;
    vector<vector<double>>* _traceFitWeights; // набор весовых коэфицентов для определения траектории движения частицы

    public:
    particle(int id){ // Конструктор TODO КК и деструктор
        this->_id = id;
        _t = new vector<double>();
        _x = new vector<double>();
        _y = new vector<double>();
        _z = new vector<double>();
        min_x = 100;
        min_y = 100;
        min_z = 100;
        max_x = -100;
        max_y = -100;
        max_z = -100;
        _traceFitWeights = new vector<vector<double>>();
        for (int i = 0; i < 3; i++)
            _traceFitWeights->push_back(new vector(0.0, 0.0));
    }
    
    double _bSpline(double t, double step) // Базовый кубический сплайн
    {
        double result = 0;
        if (t < step) result = 4/6 - pow(abs(t), 2.0) + 0.5 * pow(abs(t), 3.0);
        else if (t < step * 2) result = pow(2 - abs(t), 3.0) / 6;
        return result;
    }
    
    double _d_bSpline(double t, double step) // Первая производная _bSpline
    {
        double result = 0;
        if (t < step) result = -2 * abs(t) + 1.5 * pow(abs(t), 2.0);
        else if (t < step * 2) result = -2 + 2 * abs(t) - 0.5 * pow(abs(t), 2.0);
        return result;
    }
    double _dd_bSpline(double t, double step)  // Вторая производная _bSpline
    {
        double result = 0;
        if (t < step) result = -2 + 3 * abs(t);
        else if (t < step * 2) result = 2 - abs(t);
        return result;
    }
    double _ddd_bSpline(double t, double step)  // Третья производная _bSpline
    {
        double result = 0;
        if (t < step) result = 3.0;
        else if (t < step * 2) result = -1;
        return result;
    }

    double GetPosition(int coord, double t) // Возвращает пространственную координату в момент времени t
    {
        double result = 0;
        int i = 0;
        for (;i < _t-> length() && _t[i] != t, i++) {}
        switch(coord)
        case 1: result = _x[i];
        case 2: result = _y[i];
        case 3: result = _z[i];
        return result;
    }

    double GetMin(int coord)
    {
        double result = 0;
        switch(coord)
        case 1: result = min_x;
        case 2: result = min_y;
        case 3: result = min_z;
        return result;
    }

    double GetMax(int coord)
    {
        double result = 0;
        switch(coord)
        case 1: result = max_x;
        case 2: result = max_y;
        case 3: result = max_z;
        return result;
    }

    double Distance(double t, double x, double y, double z) // Возвращает расстояние от частицы до заданной точки
    {
        double result = 0;
        double curr_x = GetPosition(1, t);
        double curr_y = GetPosition(2, t);
        double curr_z = GetPosition(3, t);
        result = sqrt((x - curr_x)*(x - curr_x) + (y - curr_y)*(y - curr_y) + (z - curr_z)*(z - curr_z));
        return result;
    }
    
    double GetVelocity(int coord, double t) // Возвращает скорость точки вдоль координаты в момент времени t
    {
        double result = 0;
        for (int i = 0; i < _traceFitWeights[coord - 1].length(); i++)
        {
            result += _traceFitWeights[coord - 1][i] * _d_bSpline(t - _t[i]);
        }

        return result;
    }

    double GetAcceleration(int coord, double t) // Возвращает ускорение точки вдоль координаты в момент времени t
    {
        double result = 0;
        for (int i = 0; i < _traceFitWeights[coord - 1].length(); i++)
        {
            result += _traceFitWeights[coord - 1][i] * _dd_bSpline(t - _t[i]);
        }

        return result;
    }

    int GetID() // Получение ID
    {
        return _id;
    }

    void addPosition(double t, double x, double y, double z) // Добавление нового положения частицы
    {
        _t->push_back(t);
        _x->push_back(x);
        _y->push_back(y);
        _z->push_back(z);
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
        if (z < min_z) min_z = z;
        if (z > max_z) max_z = z;
        for (auto iter{_traceFitWeights->begin()}; iter != _traceFitWeights->end(); ++iter)
            *iter.push_back(0.0);
    }
    
    void print() // Вывод хранящихся начальных данных, тест корректности считывания
    {
        cout << _id << endl;
        for (auto iter{_positions->begin()}; iter != _positions->end(); ++iter)
        {
            cout << *iter[0] << " " << *iter[1] << " " << *iter[2] << " " << *iter[3] << endl;
        }
    }
    
    void TraceFit() // Поиск коэфицентов в траектории частицы
    {
        // Траектория представляется в виде линейной комбинации базовых сплайнов с весовыми коэфицентами
        // Для поиска коэфицентов используется метод наименьших квадратов: искомые коэфиценты являются экстремумом функции
        // F(c) = sum_(i=1)^(n) abs(p_c(t_i) - y_i)^2 + sum_(i=1)^(n-1) abs(lambda*p'''_c((t_i + t_i+1) / 2))
        // Из условий dF(c)/dc_i = 0 формируется и решается СЛАУ, находится искомый вектор с

        double step = *_t[1] - *_t[0];
        for (int coord = 1; i < 4; i++)
        {
            vector<double> weights(_traceFitWeights[coord-1]);
            double lambda = 1 / pow(PI * FCROSS, 3);
            Eigen::MatrixXd A(weights.length(), weights.length());
            Eigen::VectorXd b(weights.length());
            for (int i = 0; i < weights.length(); i++)
            {
                for (int j = 0; j < weights.length(); j++)
                {
                    double aij = 0;
                    for (int k = 0; k < _t->length(); k++)
                    {
                        aij += _bSpline(_t[k] - _t[j], step) * _bSpline(_t[k] - _t[i]);
                        if (k < _t->length() - 1) aij += lambda * _ddd_bSpline(_t[k] + step/2 - _t[j], step)
                        * _ddd_bSpline(_t[k] + step/2 - _t[i]);
                    }
                    double bi = 0;
                    for (int k = 0; k < _t->length(); k++)
                    {
                        bi += _bSpline(_t[k] - _t[i]) * GetPosition(coord, _t[k]);
                    }
                    A << aij;
                    b << bi;
                }
            }

            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
            Eigen::VectorXd x = qr.solve(b);
            for (int i = 0; i < weights.length(); i++)
                _traceFitWeights[coord-1][i] = x(i);
        }
        
    }
};