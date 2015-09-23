#include "SimG4CMS/Forward/interface/TotemHit.h"

TotemHit::TotemHit(int det1, float eloss1, float pabs1, int pid1, int trackid1,
                   int parent1, float x1, float y1, float z1, float vx1, float vy1,
                   float vz1, float px1, float py1, float pz1, float vpx1, float vpy1, float vpz1) :
        _det(det1),
        _eloss(eloss1),
        _pabs(pabs1),
        _pid(pid1),
        _trackid(trackid1),
        _parent(parent1),
        _x(x1),
        _y(y1),
        _z(z1),
        _vx(vx1),
        _vy(vy1),
        _vz(vz1),
        _px(px1),
        _py(py1),
        _pz(pz1),
        _vpx(vpx1),
        _vpy(vpy1),
        _vpz(vpz1) {
}

std::ostream &operator<<(std::ostream &os, const TotemHit &hit) {
    os << "TotemHit: x=" << hit.x() << " y=" << hit.y() << " z=" << hit.z() << std::endl;
    return os;
}

void TotemHit::Print() const {
    std::cout << (*this);
}

