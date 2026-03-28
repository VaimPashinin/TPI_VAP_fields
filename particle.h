#include <vector>

struct _position
{
    _position(double t, double x, double y, double z){}
    double t;
    double x;
    double y;
    double z;
};

class particle
{
    int _id;
    std::vector<_position> _positions;

    public:
    particle(int id) {}
    void addPosition(double t, double x, double y, double z) {}
    void print() {}
};
