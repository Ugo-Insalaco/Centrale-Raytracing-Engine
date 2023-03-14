using namespace std;

static inline double sqr(double x) { return x * x; }

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	double& operator[](int i) { return coord[i]; }
	double operator[](int i) const { return coord[i]; }

	Vector& operator+=(const Vector& v) {
		coord[0] += v[0];
		coord[1] += v[1];
		coord[2] += v[2];
		return *this;
	}

	double norm2() const {
		return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
	}

    double normalize() {
        double n = sqrt(this->norm2());
        coord[0]/= n;
        coord[1]/= n;
        coord[2]/= n;
        return n;
    }

    void print(){
        cout << "x: " << coord[0] << " y: " << coord[1] << " z: " << coord[2] << endl;
    }

	int getLongestAxis(){
		if(coord[0] >= coord[1] && coord[0] >= coord[2]) {
			return 0;
		}
		if(coord[1] >= coord[0] && coord[1] >= coord[2]) {
			return 1;
		}
		return 2;
	}
	double coord[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, Vector& b){
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(
    Vector& a, double b){
    return Vector(a[0]/b, a[1]/b, a[2]/b);
}

double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b){
    return Vector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

Vector computeBarycenter(const Vector& a, const Vector& b, const Vector& c){
	return 0.33*Vector(a[0]+b[0]+c[0],a[1]+b[1]+c[1],a[2]+b[2]+c[2]);
}

void translate(Vector& vector, Vector translation){
	vector[0] += translation[0];
	vector[1] += translation[1];
	vector[2] += translation[2];
}

void rotateVector(Vector& vector, Vector center, Vector rotation){
		float alpha = rotation[0];
        float beta = rotation[1];
        float gamma = rotation[2];
        translate(vector, -1*center);
		vector[0] = vector[0];
		vector[1] = cos(alpha)*vector[1]-sin(alpha)*vector[2];
		vector[2] = sin(alpha)*vector[1]+cos(alpha)*vector[2];

		vector[0] = cos(beta)*vector[0]-sin(beta)*vector[2];
		vector[1] = vector[1];
		vector[2] = sin(beta)*vector[0]+cos(beta)*vector[2];

		vector[0] = cos(gamma)*vector[0]-sin(gamma)*vector[1];
		vector[1] = sin(gamma)*vector[0]+cos(gamma)*vector[1];
		vector[2] = vector[2]; 
        translate(vector, center);
	}