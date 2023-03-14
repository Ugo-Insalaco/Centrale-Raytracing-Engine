#ifndef VECTOR_H
#define VECTOR_H
class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0);
	double& operator[](int i);
	double operator[](int i) const;
	Vector& operator+=(const Vector& v);
	double norm2 () const;
    double normalize();
    void print();
	int getLongestAxis();
	double coord[3];
};
#endif