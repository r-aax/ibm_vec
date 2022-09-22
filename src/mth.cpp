// Математические примитивы.

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <omp.h>

#include "mth.h"
#include "conf.h"

using namespace std;

// Расстояние до точки.
float
dist_to_point(float x,
              float y,
              float z,
              float cx,
              float cy,
              float cz)
{
    return sqrt(MTH_DIST2(x, y, z, cx, cy, cz));
}

// Расстояние до сферы.
float
dist_to_sphere(float x,
               float y,
               float z,
               float sx,
               float sy,
               float sz,
               float r)
{
    float d = MTH_DIST2(x, y, z, sx, sy, sz);

    d = sqrt(d);

    return abs(d - r);
}

// Печать вектора длины 4.
void
m4x4_print_vec(float (&v)[4])
{
    cout << "Vector 4 : " << endl;
    cout << "--------------------------------" << endl;

    for (int j = 0; j < 4; j++)
    {
        cout << setw(10) << v[j] << " ";
    }

    cout << endl;

    cout << "--------------------------------" << endl;
}

// Печать матрицы 4*4.
void
m4x4_print(float (&m)[4][4])
{
    cout << "Matrix 4*4 : " << endl;
    cout << "--------------------------------" << endl;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << setw(10) << m[i][j] << " ";
        }

        cout << endl;
    }

    cout << "--------------------------------" << endl;
}

// Печать двух матриц 4*4.
void
m4x4_print_duplex(float (&a)[4][4],
                  float (&b)[4][4])
{
    cout << "Matrixes 4*4 duplex : " << endl;
    cout << "----------------------------------------------------" << endl;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << setw(10) << a[i][j] << " ";
        }

        cout << "\t\t\t";

        for (int j = 0; j < 4; j++)
        {
            cout << setw(10) << b[i][j] << " ";
        }

        cout << endl;
    }

    cout << "----------------------------------------------------" << endl;
}

// Инициализация вектора.
void
m4x4_init_vec(float (&v)[4],
              float v0,
              float v1,
              float v2,
              float v3)
{
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
}

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
              float m33)
{
    m[0][0] = m00;
    m[0][1] = m01;
    m[0][2] = m02;
    m[0][3] = m03;
    m[1][0] = m10;
    m[1][1] = m11;
    m[1][2] = m12;
    m[1][3] = m13;
    m[2][0] = m20;
    m[2][1] = m21;
    m[2][2] = m22;
    m[2][3] = m23;
    m[3][0] = m30;
    m[3][1] = m31;
    m[3][2] = m32;
    m[3][3] = m33;
}

// Скалярное произведение векторов.
float
m4x4_scalar_product(float (&a)[4],
                    float (&b)[4])
{
    float r = 0.0;

    for (int i = 0; i < 4; i++)
    {
        r += a[i] * b[i];
    }

    return r;
}

// Умножение вектор-строки на матрицу.
void
m4x4_mul_vec_mat(float (&v)[4],
                 float (&m)[4][4],
                 float (&r)[4])
{
    for (int j = 0; j < 4; j++)
    {
        float f = 0.0;

        for (int k = 0; k < 4; k++)
        {
            f += v[k] * m[k][j];
        }

        r[j] = f;
    }
}

// Умножение матрицы на вектор-столбец.
void
m4x4_mul_mat_vec(float (&m)[4][4],
                 float (&v)[4],
                 float (&r)[4])
{
    for (int i = 0; i < 4; i++)
    {
        float f = 0.0;

        for (int k = 0; k < 4; k++)
        {
            f += m[i][k] * v[k];
        }

        r[i] = f;
    }
}

// Перемножение матриц 4*4.
void
m4x4_mul(float (&a)[4][4],
         float (&b)[4][4],
         float (&c)[4][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            float f = 0.0;

            for (int k = 0; k < 4; k++)
            {
                f += a[i][k] * b[k][j];
            }

            c[i][j] = f;
        }
    }
}

// Инициализация единичной матрицы.
void
m4x4_init_E(float (&m)[4][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            m[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

// Копирование матрицы.
void
m4x4_copy(float (&src)[4][4],
          float (&dst)[4][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            dst[i][j] = src[i][j];
        }
    }
}

// Умножение строки на число.
void
m4x4_div_line(float (&m)[4][4],
              int i,
              float k)
{
    for (int j = 0; j < 4; j++)
    {
        m[i][j] /= k;
    }
}

// Перестановка в матрице 4*4 двух строк.
void
m4x4_swap_lines(float (&m)[4][4],
                int i1,
                int i2)
{
    if (i1 != i2)
    {
        for (int j = 0; j < 4; j++)
        {
            float f = m[i1][j];

            m[i1][j] = m[i2][j];
            m[i2][j] = f;
        }
    }
}

// Прибаление к одной строки другой строки, умноженной на число.
void
m4x4_add_to_1_line_2k(float (&m)[4][4],
                      int i1,
                      int i2,
                      float k)
{
    for (int j = 0; j < 4; j++)
    {
        m[i1][j] += m[i2][j] * k;
    }
}

// Номер ведущей строки по данному столбцу.
// Ищется только среди строк с номерами, не меньшими, чем j.
int
m4x4_lead_line(float (&m)[4][4],
               int j)
{
    int lead_i = j;
    float lead_v = abs(m[j][j]);

    for (int i = j + 1; i < 4; i++)
    {
        float v = abs(m[i][j]);

        if (v > lead_v)
        {
            lead_i = i;
            lead_v = v;
        }
    }

    return lead_i;
}

// Инвертирование матрицы 4*4.
bool
m4x4_invert(float (&a)[4][4],
            float (&b)[4][4])
{
    float t[4][4];

    m4x4_copy(a, t);
    m4x4_init_E(b);

#if MTH_DEBUG_PRINT == 1
    cout << "INI" << endl;
    m4x4_print_duplex(t, b);
#endif

    for (int si = 0; si < 4; si++)
    {
        int lead_i = m4x4_lead_line(t, si);

        // Ведущий элемент на первое место.
        m4x4_swap_lines(t, si, lead_i);
        m4x4_swap_lines(b, si, lead_i);

#if MTH_DEBUG_PRINT == 1
        cout << "SWP: " << si << ", " << lead_i << endl;
        m4x4_print_duplex(t, b);
#endif

        float tsisi = t[si][si];

        if (abs(tsisi) < 1.0e-10)
        {
            return false;
        }

        // Зануляем все, что ниже и выше.
        for (int ti = 0; ti < 4; ti++)
        {
            if (ti != si)
            {
                float k = -t[ti][si] / tsisi;

                m4x4_add_to_1_line_2k(t, ti, si, k);
                m4x4_add_to_1_line_2k(b, ti, si, k);

#if MTH_DEBUG_PRINT == 1
                cout << "FMA: to " << ti << " add " << si << " mul " << k << endl;
                m4x4_print_duplex(t, b);
#endif

            }
        }

        // Саму строку нормируем.
        m4x4_div_line(t, si, tsisi);
        m4x4_div_line(b, si, tsisi);

#if MTH_DEBUG_PRINT == 1
        cout << "DIV: " << si << " on " << tsisi << endl;
        m4x4_print_duplex(t, b);
#endif

    }

    return true;
}
