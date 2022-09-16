// Математические примитивы.

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <fstream>
#include <math.h>

#include "mth.h"

using namespace std;

// Расстояние до точки.
double
dist_to_point(double x,
              double y,
              double z,
              double cx,
              double cy,
              double cz)
{
    return sqrt(MTH_DIST2(x, y, z, cx, cy, cz));
}

// Расстояние до сферы.
double
dist_to_sphere(double x,
               double y,
               double z,
               double sx,
               double sy,
               double sz,
               double r)
{
    double d = MTH_DIST2(x, y, z, sx, sy, sz);

    d = sqrt(d);

    return abs(d - r);
}

// Печать матрицы 4*4.
void
m4x4_print(double (&m)[4][4])
{
    cout << "Matrix 4*4 : " << endl;
    cout << "--------------------------------" << endl;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << m[i][j] << " ";
        }

        cout << endl;
    }

    cout << "--------------------------------" << endl;
}

// Печать двух матриц 4*4.
void
m4x4_print_duplex(double (&a)[4][4],
                  double (&b)[4][4])
{
    cout << "Matrixes 4*4 duplex : " << endl;
    cout << "----------------------------------------------------" << endl;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << a[i][j] << " ";
        }

        cout << "\t\t\t";

        for (int j = 0; j < 4; j++)
        {
            cout << b[i][j] << " ";
        }

        cout << endl;
    }

    cout << "----------------------------------------------------" << endl;
}

// Перемножение матриц 4*4.
void
m4x4_mul(double (&a)[4][4],
         double (&b)[4][4],
         double (&c)[4][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            double f = 0.0;

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
m4x4_init_E(double (&m)[4][4])
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
m4x4_copy(double (&src)[4][4],
          double (&dst)[4][4])
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
m4x4_div_line(double (&m)[4][4],
              int i,
              double k)
{
    for (int j = 0; j < 4; j++)
    {
        m[i][j] /= k;
    }
}

// Перестановка в матрице 4*4 двух строк.
void
m4x4_swap_lines(double (&m)[4][4],
                int i1,
                int i2)
{
    if (i1 != i2)
    {
        for (int j = 0; j < 4; j++)
        {
            double f = m[i1][j];

            m[i1][j] = m[i2][j];
            m[i2][j] = f;
        }
    }
}

// Прибаление к одной строки другой строки, умноженной на число.
void
m4x4_add_to_1_line_2k(double (&m)[4][4],
                      int i1,
                      int i2,
                      double k)
{
    for (int j = 0; j < 4; j++)
    {
        m[i1][j] += m[i2][j] * k;
    }
}

// Номер ведущей строки по данному столбцу.
// Ищется только среди строк с номерами, не меньшими, чем j.
int
m4x4_lead_line(double (&m)[4][4],
               int j)
{
    int lead_i = j;
    double lead_v = abs(m[j][j]);

    for (int i = j + 1; i < 4; i++)
    {
        double v = abs(m[i][j]);

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
m4x4_invert(double (&a)[4][4],
            double (&b)[4][4])
{
    double t[4][4];

    m4x4_copy(a, t);
    m4x4_init_E(b);

    for (int si = 0; si < 4; si++)
    {
        int lead_i = m4x4_lead_line(t, si);

        // Ведущий элемент на первое место.
        m4x4_swap_lines(t, si, lead_i);
        m4x4_swap_lines(b, si, lead_i);

        double tsisi = t[si][si];

        if (tsisi < 1.0e-10)
        {
            return false;
        }

        // Зануляем все, что ниже и выше.
        for (int ti = 0; ti < 4; ti++)
        {
            if (ti != si)
            {
                double k = -t[ti][si] / tsisi;

                m4x4_add_to_1_line_2k(t, ti, si, k);
                m4x4_add_to_1_line_2k(b, ti, si, k);
            }
        }

        // Саму строку нормируем.
        m4x4_div_line(t, si, tsisi);
        m4x4_div_line(b, si, tsisi);
    }

    return true;
}
