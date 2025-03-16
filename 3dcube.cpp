#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

struct Matrix3x3 {
    double m[3][3];
};

struct Vec {
    double x;
    double y;
    double z;

    Vec(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}

    void applyMatrix(const Matrix3x3& mat) {
        double newX = mat.m[0][0] * x + mat.m[0][1] * y + mat.m[0][2] * z;
        double newY = mat.m[1][0] * x + mat.m[1][1] * y + mat.m[1][2] * z;
        double newZ = mat.m[2][0] * x + mat.m[2][1] * y + mat.m[2][2] * z;

        x = newX;
        y = newY;
        z = newZ;
    }
};

bool isOnCubeFace(double x, double y, double z, double a, double epsilon) {
    return 
    (fabs(fabs(x) - a) <= epsilon && y >= -a - epsilon && y <= a + epsilon && z >= -a - epsilon && z <= a + epsilon) ||
    (fabs(fabs(y) - a) <= epsilon && x >= -a - epsilon && x <= a + epsilon && z >= -a - epsilon && z <= a + epsilon) ||
    (fabs(fabs(z) - a) <= epsilon && x >= -a - epsilon && x <= a + epsilon && y >= -a - epsilon && y <= a + epsilon);
}

void generate_cube(vector<Vec>& vectors, int edge_points) {
    const int points_per_edge = edge_points * 2 + 1;
    const double step = 2.0 / (points_per_edge - 1);
    const double a = 1.0;

    for (int i = 0; i < points_per_edge; ++i) {
        double x = -a + i * step;
        for (int j = 0; j < points_per_edge; ++j) {
            double y = -a + j * step;
            for (int k = 0; k < points_per_edge; ++k) {
                double z = -a + k * step;
                if (isOnCubeFace(x, y, z, a, 1e-4)) {
                    vectors.emplace_back(x, y, z);
                }
            }
        }
    }
}

Matrix3x3 createRotationMatrix(double alpha, double beta, double gamma) {
    double radAlpha = alpha * M_PI / 180.0;
    double radBeta = beta * M_PI / 180.0;
    double radGamma = gamma * M_PI / 180.0;

    double cosA = cos(radAlpha);
    double sinA = sin(radAlpha);
    double cosB = cos(radBeta);
    double sinB = sin(radBeta);
    double cosG = cos(radGamma);
    double sinG = sin(radGamma);

    Matrix3x3 mat = {{
        {cosB * cosG, -cosB * sinG, sinB},
        {sinA * sinB * cosG + cosA * sinG, -sinA * sinB * sinG + cosA * cosG, -sinA * cosB},
        {-cosA * sinB * cosG + sinA * sinG, cosA * sinB * sinG + sinA * cosG, cosA * cosB}
    }};

    return mat;
}

void rotateVectors(double alpha, double beta, double gamma, vector<Vec>& vectors) {
    Matrix3x3 mat = createRotationMatrix(alpha, beta, gamma);
    for (Vec& v : vectors) {
        v.applyMatrix(mat);
    }
}

bool compareByZ(const Vec& v1, const Vec& v2) {
    return v1.z > v2.z;
}

void sortVectors(vector<Vec>& vectors) {
    sort(vectors.begin(), vectors.end(), compareByZ);
}