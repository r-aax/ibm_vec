// Математические примитивы.

#ifndef MTH_H
#define MTH_H

// Среднее значение.
#define MTH_AVG(A, B) (0.5 * ((A) + (B)))

// Квадрат разницы между значениями.
#define MTH_DIFF2(A, B) (((A) - (B)) * ((A) - (B)))

// Квадрат рассмояния между точками.
#define MTH_DIST2(X1, Y1, Z1, X2, Y2, Z2) (MTH_DIFF2(X1, X2) + MTH_DIFF2(Y1, Y2) + MTH_DIFF2(Z1, Z2))

// Расстояние до точки.
double
dist_to_point(double x,
              double y,
              double z,
              double cx,
              double cy,
              double cz);

// Расстояние до сферы.
double
dist_to_sphere(double x,
               double y,
               double z,
               double sx,
               double sy,
               double sz,
               double r);

// Печать матрицы 4*4.
void
m4x4_print(double (&m)[4][4]);

// Печать двух матриц 4*4.
void
m4x4_print_duplex(double (&a)[4][4],
                  double (&b)[4][4]);

// Инициализация вектора.
void
m4x4_init_vec(double (&v)[4],
              double v0,
              double v1,
              double v2,
              double v3);

// Инициализация матрицы.
void
m4x4_init_mat(double (&m)[4][4],
              double m00,
              double m01,
              double m02,
              double m03,
              double m10,
              double m11,
              double m12,
              double m13,
              double m20,
              double m21,
              double m22,
              double m23,
              double m30,
              double m31,
              double m32,
              double m33);

// Скалярное произведение векторов.
double
m4x4_scalar_product(double (&a)[4],
                    double (&b)[4]);

// Умножение матрицы на вектор-столбец.
void
m4x4_mul_on_vec(double (&m)[4][4],
                double (&v)[4],
                double (&r)[4]);

// Перемножение матриц 4*4.
void
m4x4_mul(double (&a)[4][4],
         double (&b)[4][4],
         double (&c)[4][4]);

// Инициализация единичной матрицы.
void
m4x4_init_E(double (&m)[4][4]);

// Копирование матрицы.
void
m4x4_copy(double (&src)[4][4],
          double (&dst)[4][4]);

// Умножение строки на число.
void
m4x4_div_line(double (&m)[4][4],
              int i,
              double k);

// Перестановка в матрице 4*4 двух строк.
void
m4x4_swap_lines(double (&m)[4][4],
                int i1,
                int i2);

// Прибаление к одной строки другой строки, умноженной на число.
void
m4x4_add_to_1_line_2k(double (&m)[4][4],
                      int i1,
                      int i2,
                      double k);

// Номер ведущей строки по данному столбцу.
int
m4x4_lead_line(double (&m)[4][4],
               int j);

// Инвертирование матрицы 4*4.
bool
m4x4_invert(double (&a)[4][4],
            double (&b)[4][4]);

#endif
