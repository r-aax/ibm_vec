#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <fstream>
#include <math.h>

#include "conf.h"
#include "mth.h"

using namespace std;

// Реализация численного метода.
//
// 1. Разностная схема
//    Н.С. Смирнова. Сравнение схем с расщеплением потока для численного решения уравнений Эйлера сжимаемого газа.
//    (ТРУДЫ МФТИ, 2018, Том 10, № 1).

// Количество ячеек.
#define CELLS_COUNT (NX * NY * NZ)

// Линеаризация индекса.
#define LIN(IX, IY, IZ) ((IZ * NY + IY) * NX + IX)

// Определение циклов.
#define LOOP1 for (int i = 0; i < CELLS_COUNT; i++)
#define LOOP3 for (int iz = 0; iz < NZ; iz++) for (int iy = 0; iy < NY; iy++) for (int ix = 0; ix < NX; ix++)

// Данные расчетной области.
// Примитивные данные:
//   - плотность,
//   - составляющие скорости по x, y, z,
//   - давление.
// Сохраняющиеся данные:
//   - плотность импульса
//   - плотность полной энергии
float r[CELLS_COUNT];
float u[CELLS_COUNT];
float v[CELLS_COUNT];
float w[CELLS_COUNT];
float p[CELLS_COUNT];
//
float ru[CELLS_COUNT];
float rv[CELLS_COUNT];
float rw[CELLS_COUNT];
float E[CELLS_COUNT];
//
float fp_r[CELLS_COUNT];
float fp_ru[CELLS_COUNT];
float fp_rv[CELLS_COUNT];
float fp_rw[CELLS_COUNT];
float fp_E[CELLS_COUNT];
float fn_r[CELLS_COUNT];
float fn_ru[CELLS_COUNT];
float fn_rv[CELLS_COUNT];
float fn_rw[CELLS_COUNT];
float fn_E[CELLS_COUNT];

// Вспомогательные данные расчетной области.
// Тип ячейки:
//   0 - обычная расчетная ячейка.
//   1 - граничная ячейка.
//   2 - фиктивная ячейка.
//   3 - внутренняя ячейка, которая не принимает участие в расчетах.
#define KIND_COMMON 0
#define KIND_BORDER 1
#define KIND_GHOST 2
#define KIND_INNER 3
int kind[CELLS_COUNT];

// Данные сфер для обтекания.
float sph_x[SPHERES_COUNT];
float sph_y[SPHERES_COUNT];
float sph_z[SPHERES_COUNT];
float sph_r[SPHERES_COUNT];

// Инициализация сферы.
void
sphere_init_xy(int i,
               float sphere_x,
               float sphere_y,
               float sphere_r)
{
    sph_x[i] = sphere_x;
    sph_y[i] = sphere_y;
    sph_z[i] = (static_cast<float>(NZ) / 2.0);
    sph_r[i] = sphere_r;
}

// Инициализация расчетной области.
void
calc_area_init()
{
    int left_cell_count = 1;

    // Н.С. Смирнова. Сравнение схем с расщеплением потока для численного решения уравнений Эйлера сжимаемого газа.
    // (ТРУДЫ МФТИ, 2018, Том 10, № 1).
    // Таблица 1. Тест 1.
    LOOP3
    {
        int i = LIN(ix, iy, iz);

        if (ix < left_cell_count)
        {
            r[i] = 1.0;
            u[i] = 0.75;
            v[i] = 0.0;
            w[i] = 0.0;
            p[i] = 1.0;
        }
        else
        {
            r[i] = 0.125;
            u[i] = 0.0;
            v[i] = 0.0;
            w[i] = 0.0;
            p[i] = 0.1;
        }
    }

    // Инициализация сфер.
    sphere_init_xy(0,  2.0, 2.0, 1.05);
    sphere_init_xy(1,  4.0, 4.0, 0.75);
    sphere_init_xy(2,  5.0, 1.0, 0.75);
    sphere_init_xy(3,  7.0, 2.5, 1.15);
    sphere_init_xy(4, 10.0, 0.0, 1.05);
}

// Определение типов ячеек.
void
calc_area_define_cells_kinds()
{
    LOOP1 kind[i] = KIND_COMMON;

    LOOP3
    {
        int i = LIN(ix, iy, iz);
        float xl = DH * static_cast<float>(ix);
        float xr = DH * static_cast<float>(ix + 1);
        float yl = DH * static_cast<float>(iy);
        float yr = DH * static_cast<float>(iy + 1);
        float zl = DH * static_cast<float>(iz);
        float zr = DH * static_cast<float>(iz + 1);

        for (int si = 0; si < SPHERES_COUNT; si++)
        {
            float dxl2 = MTH_DIFF2(xl, sph_x[si]);
            float dxr2 = MTH_DIFF2(xr, sph_x[si]);
            float dyl2 = MTH_DIFF2(yl, sph_y[si]);
            float dyr2 = MTH_DIFF2(yr, sph_y[si]);
            float dzl2 = MTH_DIFF2(zl, sph_z[si]);
            float dzr2 = MTH_DIFF2(zr, sph_z[si]);
            float r2 = sph_r[si] * sph_r[si];
            int lll_in = static_cast<int>(dxl2 + dyl2 + dzl2 <= r2);
            int llr_in = static_cast<int>(dxl2 + dyl2 + dzr2 <= r2);
            int lrl_in = static_cast<int>(dxl2 + dyr2 + dzl2 <= r2);
            int lrr_in = static_cast<int>(dxl2 + dyr2 + dzr2 <= r2);
            int rll_in = static_cast<int>(dxr2 + dyl2 + dzl2 <= r2);
            int rlr_in = static_cast<int>(dxr2 + dyl2 + dzr2 <= r2);
            int rrl_in = static_cast<int>(dxr2 + dyr2 + dzl2 <= r2);
            int rrr_in = static_cast<int>(dxr2 + dyr2 + dzr2 <= r2);
            int in_cnt = lll_in + llr_in + lrl_in + lrr_in + rll_in + rlr_in + rrl_in + rrr_in;

            if (in_cnt == 0)
            {
                // Все ячейки снаружи окружности.
                // Ничего не делаем.
                ; 
            }
            else if (in_cnt == 8)
            {
                // Все ячейки внутри окружности.
                kind[i] = KIND_INNER;

                // Конфликт может быть только с одной сферой.
                break;
            }
            else
            {
                float xc = MTH_AVG(xl, xr);
                float yc = MTH_AVG(yl, yr);
                float zc = MTH_AVG(zl, zr);

                if (MTH_DIST2(xc, yc, zc, sph_x[si], sph_y[si], sph_z[si]) <= r2)
                {
                    // Центр ячейки внутри сферы.
                    kind[i] = KIND_GHOST;
                }
                else
                {
                    // Центр ячейки снаружи сферы.
                    kind[i] = KIND_BORDER;
                }

                // Конфликт может быть только с одной сферой.
                break;
            }
        }
    }

    // Еще один проход по внутренним ячейкам.
    // Если соседом внутренней ячейки является граничная или обычная ячейка,
    // то внутренняя ячейка становится фиктивной.
    LOOP3
    {
        int i = LIN(ix, iy, iz);

        if (kind[i] == KIND_INNER)
        {
            if (ix > 0)
            {
                int bi = LIN(ix - 1, iy, iz);

                if ((kind[bi] == KIND_COMMON) || (kind[bi] == KIND_BORDER))
                {
                    kind[i] = KIND_GHOST;

                    continue;
                }
            }

            if (ix < NX - 1)
            {
                int bi = LIN(ix + 1, iy, iz);

                if ((kind[bi] == KIND_COMMON) || (kind[bi] == KIND_BORDER))
                {
                    kind[i] = KIND_GHOST;

                    continue;
                }
            }

            if (iy > 0)
            {
                int bi = LIN(ix, iy - 1, iz);

                if ((kind[bi] == KIND_COMMON) || (kind[bi] == KIND_BORDER))
                {
                    kind[i] = KIND_GHOST;

                    continue;
                }
            }

            if (iy < NY - 1)
            {
                int bi = LIN(ix, iy + 1, iz);

                if ((kind[bi] == KIND_COMMON) || (kind[bi] == KIND_BORDER))
                {
                    kind[i] = KIND_GHOST;

                    continue;
                }
            }

            if (iz > 0)
            {
                int bi = LIN(ix, iy, iz - 1);

                if ((kind[bi] == KIND_COMMON) || (kind[bi] == KIND_BORDER))
                {
                    kind[i] = KIND_GHOST;

                    continue;
                }
            }

            if (iz < NZ - 1)
            {
                int bi = LIN(ix, iy, iz + 1);

                if ((kind[bi] == KIND_COMMON) || (kind[bi] == KIND_BORDER))
                {
                    kind[i] = KIND_GHOST;

                    continue;
                }
            }
        }
    }
}

// Экспорт в ParaView.
void
calc_area_paraview_export(int i)
{
    ostringstream ss;
    ss << "out/export_" << setfill('0') << setw(10) << i << ".dat";
    ofstream f(ss.str());

    f << "TITLE=\"[" << NX << " * " << NY << " * " << NZ << "] calc area\"" << endl;
    f << "VARIABLES=\"X\", \"Y\", \"Z\", \"Rho\", \"U\", \"V\", \"W\", \"P\", \"Kind\"" << endl;
    f << "ZONE T=\"single zone\"" << endl;
    f << "NODES=" << (8 * CELLS_COUNT) << endl;
    f << "ELEMENTS=" << CELLS_COUNT << endl;
    f << "DATAPACKING=BLOCK" << endl;
    f << "ZONETYPE=FEBRICK" << endl;
    f << "VARLOCATION=([4-9]=CELLCENTERED)" << endl;

#define LX (ix * DH)
#define HX ((ix + 1) * DH)
#define LY (iy * DH)
#define HY ((iy + 1) * DH)
#define LZ (iz * DH)
#define HZ ((iz + 1) * DH)

    LOOP3 f << LX << " " << HX << " " << LX << " " << HX << " " << LX << " " << HX << " " << LX << " " << HX << " "; f << endl;
    LOOP3 f << LY << " " << LY << " " << HY << " " << HY << " " << LY << " " << LY << " " << HY << " " << HY << " "; f << endl;
    LOOP3 f << LZ << " " << LZ << " " << LZ << " " << LZ << " " << HZ << " " << HZ << " " << HZ << " " << HZ << " "; f << endl;
    LOOP1 f << r[i] << " "; f << endl;
    LOOP1 f << u[i] << " "; f << endl;
    LOOP1 f << v[i] << " "; f << endl;
    LOOP1 f << w[i] << " "; f << endl;
    LOOP1 f << p[i] << " "; f << endl;
    LOOP1 f << kind[i] << " "; f << endl;

#undef LX
#undef HX
#undef LY
#undef HY
#undef LZ
#undef HZ

    LOOP1
    {
        int s = 8 * i;

        // Нумерация точек в элементе начинается с единицы.
        f << (s + 1) << " " << (s + 2) << " " << (s + 4) << " " << (s + 3) << " "
          << (s + 5) << " " << (s + 6) << " " << (s + 8) << " " << (s + 7) << endl;
    }

    f.close();
}

// Перевод примитивных величин в консервативные.
void
d_to_u()
{
    LOOP1
    {
        ru[i] = r[i] * u[i];
        rv[i] = r[i] * v[i];
        rw[i] = r[i] * w[i];

        // E = rho * (V^2/2 + e) = rho * (V^2/2 + p/((GAMMA - 1.0) * rho))
        //   = rho * V^2/2 + p/(GAMMA - 1)
        E[i] = 0.5 * r[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]) + p[i] / (GAMMA - 1.0);
    }
}

// Перевод консервативных величин в примитивные.
void
u_to_d()
{
    LOOP1
    {
        u[i] = ru[i] / r[i];
        v[i] = rv[i] / r[i];
        w[i] = rw[i] / r[i];

        // E = rho * V^2/2 + p/(GAMMA - 1.0)
        // p = (E - rho * V^2/2) * (GAMMA - 1.0)
        p[i] = (E[i] - 0.5 * r[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i])) * (GAMMA - 1.0);
    }
}

// Вычисление f/g/h pos/neg.
void
calc_f()
{
    LOOP1
    {
        float a = sqrt(GAMMA * p[i] / r[i]);
        float l1 = u[i] - a;
        float l2 = u[i];
        float l5 = u[i] + a;
        float lp1 = 0.5 * (l1 + abs(l1));
        float lp2 = 0.5 * (l2 + abs(l2));
        float lp5 = 0.5 * (l5 + abs(l5));
        float ln1 = 0.5 * (l1 - abs(l1));
        float ln2 = 0.5 * (l2 - abs(l2));
        float ln5 = 0.5 * (l5 - abs(l5));
        float k = 0.5 * r[i] / GAMMA;
        float V2 = u[i] * u[i] + v[i] * v[i] + w[i] * w[i];
        float H = 0.5 * V2 + a * a / (GAMMA - 1.0);

        fp_r[i] = k * (lp1 + 2.0 * (GAMMA - 1.0) * lp2 + lp5);
        fp_ru[i] = k * ((u[i] - a) * lp1 + 2.0 * (GAMMA - 1.0) * u[i] * lp2 + (u[i] + a) * lp5);
        fp_rv[i] = k * (v[i] * lp1 + 2.0 * (GAMMA - 1.0) * v[i] * lp2 + v[i] * lp5);
        fp_rw[i] = k * (w[i] * lp1 + 2.0 * (GAMMA - 1.0) * w[i] * lp2 + w[i] * lp5);
        fp_E[i] = k * ((H - u[i] * a) * lp1 + (GAMMA - 1.0) * V2 * lp2 + (H + u[i] * a) * lp5);

        fn_r[i] = k * (ln1 + 2.0 * (GAMMA - 1.0) * ln2 + ln5);
        fn_ru[i] = k * ((u[i] - a) * ln1 + 2.0 * (GAMMA - 1.0) * u[i] * ln2 + (u[i] + a) * ln5);
        fn_rv[i] = k * (v[i] * ln1 + 2.0 * (GAMMA - 1.0) * v[i] * ln2 + v[i] * ln5);
        fn_rw[i] = k * (w[i] * ln1 + 2.0 * (GAMMA - 1.0) * w[i] * ln2 + w[i] * ln5);
        fn_E[i] = k * ((H - u[i] * a) * ln1 + (GAMMA - 1.0) * V2 * ln2 + (H + u[i] * a) * ln5);
    }
}

// Вычисление f/g/h pos/neg.
void
calc_fgh()
{
    calc_f();
}

// Вычисление потоков.
void
calc_flows()
{
    calc_fgh();

    LOOP3
    {
        //
        // Направление X.
        //

        int i = LIN(ix, iy, iz);
        int li = i;
        int ri = i;

        if (ix > 0)
        {
            li = LIN(ix - 1, iy, iz);
        }

        if (ix < NX - 1)
        {
            ri = LIN(ix + 1, iy, iz);
        }

        r[i] -= (DT / DH) * (fp_r[i] + fn_r[ri] - fp_r[li] - fn_r[i]);
        ru[i] -= (DT / DH) * (fp_ru[i] + fn_ru[ri] - fp_ru[li] - fn_ru[i]);
        rv[i] -= (DT / DH) * (fp_rv[i] + fn_rv[ri] - fp_rv[li] - fn_rv[i]);
        rw[i] -= (DT / DH) * (fp_rw[i] + fn_rw[ri] - fp_rw[li] - fn_rw[i]);
        E[i] -= (DT / DH) * (fp_E[i] + fn_E[ri] - fp_E[li] - fn_E[i]);
    }
}

// Шаг вычислений.
void
step()
{
    d_to_u();
    calc_flows();
    u_to_d();
}

// Точка входа.
int
main()
{
    calc_area_init();
    calc_area_define_cells_kinds();
    calc_area_paraview_export(0);

    for (int i = 0; i < TIME_STEPS; i++)
    {
        cout << ".... step " << i << " of " << TIME_STEPS << endl;

        step();
        calc_area_paraview_export(i + 1);
    }
}
