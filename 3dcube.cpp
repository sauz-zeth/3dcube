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

bool isOnTorusSurface(
    double x, double y, double z, 
    double R, double r, 
    double epsilon
) {
    // Альтернативная форма уравнения тора (без корня):
    // (x² + y² + z² + R² - r²)² = 4R²(x² + y²)
    double left = pow(x*x + y*y + z*z + R*R - r*r, 2);
    double right = 4 * R*R * (x*x + y*y);
    return fabs(left - right) <= epsilon;
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
                // if (isOnTorusSurface(x, y, z, a, 0.5, 1e-4)) {
                //     vectors.emplace_back(x, y, z);
                // }
            }
        }
    }
}

void generateEightPetalRose(std::vector<Vec>& points, int resolution) {
    for (int i = 0; i < resolution; ++i) {
        double theta = 2 * M_PI * i / resolution; // Азимутальный угол
        for (int j = 0; j < resolution; ++j) {
            double phi = M_PI * j / resolution; // Полярный угол

            // Радиус как функция углов
            double r = 1 + sin(4 * theta) * sin(4 * phi);

            // Преобразование в декартовы координаты
            double x = r * sin(phi) * cos(theta);
            double y = r * sin(phi) * sin(theta);
            double z = r * cos(phi);

            points.emplace_back(x, y, z);
        }
    }
}

void generateBoysSurface(std::vector<Vec>& points, int resolution) {
    
    for (int i = 0; i < resolution; ++i) {
        double u = M_PI * i / resolution; // Параметр u
        for (int j = 0; j < resolution; ++j) {
            double v = 2 * M_PI * j / resolution; // Параметр v

            // Вычисление координат
            double denominator = 2 - sqrt(2) * sin(3 * u) * sin(2 * v);
            double x = (sqrt(2) * pow(cos(u), 2) * cos(2 * v) + cos(v) * sin(2 * u)) / denominator;
            double y = (sqrt(2) * pow(cos(u), 2) * sin(2 * v) - sin(v) * sin(2 * u)) / denominator;
            double z = (3 * pow(cos(u), 2)) / denominator;

            points.emplace_back(x, y, z);
        }
    }
}

void generateHyperbolicParaboloid(std::vector<Vec>& points, int resolution) {
    double step = 2.0 / (resolution - 1); // Шаг между точками
    for (int i = 0; i < resolution; ++i) {
        double x = -1.0 + i * step; // От -1 до 1
        for (int j = 0; j < resolution; ++j) {
            double y = -1.0 + j * step; // От -1 до 1
            double z = x * x - y * y;  // Уравнение гиперболического параболоида
            points.emplace_back(x, y, z);
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