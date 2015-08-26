/*
Hit class for analysis in TotemTest
*/

#ifndef TOTEMHIT_H
#define TOTEMHIT_H

#include<iostream>

class TotemHit {
public:
    TotemHit(int, float, float, int, int,
             int, float, float, float, float,
             float, float, float, float, float, float, float, float);

    TotemHit() { }

    int det() { return _det; }

    int pid() { return _pid; }

    int trackid() { return _trackid; }

    int parent() { return _parent; }

    float eloss() { return _eloss; }

    float pabs() { return _pabs; }

    float x() const { return _x; }

    float y() const { return _y; }

    float z() const { return _z; }

    float vx() { return _vx; }

    float vy() { return _vy; }

    float vz() { return _vz; }

    float px() { return _px; }

    float py() { return _py; }

    float pz() { return _pz; }

    float vpx() { return _vpx; }

    float vpy() { return _vpy; }

    float vpz() { return _vpz; }

    void Print() const;

private:
    int _det;
    float _eloss;
    float _pabs;
    int _pid;
    int _trackid;
    int _parent;
    float _x;
    float _y;
    float _z;
    float _vx;
    float _vy;
    float _vz;
    float _px;
    float _py;
    float _pz;
    float _vpx;
    float _vpy;
    float _vpz;

};

#endif
