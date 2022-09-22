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
float
dist_to_point(float x,
              float y,
              float z,
              float cx,
              float cy,
              float cz);

// Расстояние до сферы.
float
dist_to_sphere(float x,
               float y,
               float z,
               float sx,
               float sy,
               float sz,
               float r);

// Печать вектора длины 4.
void
m4x4_print_vec(float (&v)[4]);

// Печать матрицы 4*4.
void
m4x4_print(float (&m)[4][4]);

// Печать двух матриц 4*4.
void
m4x4_print_duplex(float (&a)[4][4],
                  float (&b)[4][4]);

// Инициализация вектора.
void
m4x4_init_vec(float (&v)[4],
              float v0,
              float v1,
              float v2,
              float v3);

// Инициализация матрицы.
void
m4x4_init_mat(float (&m)[4][4],
              float m00,
              float m01,
              float m02,
              float m03,
              float m10,
              float m11,
              float m12,
              float m13,
              float m20,
              float m21,
              float m22,
              float m23,
              float m30,
              float m31,
              float m32,
              float m33);

// Скалярное произведение векторов.
float
m4x4_scalar_product(float (&a)[4],
                    float (&b)[4]);

// Умножение вектор-строки на матрицу.
void
m4x4_mul_vec_mat(float (&v)[4],
                 float (&m)[4][4],
                 float (&r)[4]);

// Умножение матрицы на вектор-столбец.
void
m4x4_mul_mat_vec(float (&m)[4][4],
                 float (&v)[4],
                 float (&r)[4]);

// Перемножение матриц 4*4.
void
m4x4_mul(float (&a)[4][4],
         float (&b)[4][4],
         float (&c)[4][4]);

// Инициализация единичной матрицы.
void
m4x4_init_E(float (&m)[4][4]);

// Копирование матрицы.
void
m4x4_copy(float (&src)[4][4],
          float (&dst)[4][4]);

// Умножение строки на число.
void
m4x4_div_line(float (&m)[4][4],
              int i,
              float k);

// Перестановка в матрице 4*4 двух строк.
void
m4x4_swap_lines(float (&m)[4][4],
                int i1,
                int i2);

// Прибаление к одной строки другой строки, умноженной на число.
void
m4x4_add_to_1_line_2k(float (&m)[4][4],
                      int i1,
                      int i2,
                      float k);

// Номер ведущей строки по данному столбцу.
int
m4x4_lead_line(float (&m)[4][4],
               int j);

// Инвертирование матрицы 4*4.
bool
m4x4_invert(float (&a)[4][4],
            float (&b)[4][4]);

#endif
