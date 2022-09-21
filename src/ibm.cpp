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
#define LIN(IX, IY, IZ) (((IZ) * NY + (IY)) * NX + (IX))

// Получение координат из линеаризованного номера.
#define UNLINX(I) ((I) % (NX))
#define UNLINY(I) (((I) / (NX)) % (NY))
#define UNLINZ(I) ((I) / ((NX) * (NY)))

// Перевод индекса в координату центра.
#define CCORD(I) ((static_cast<double>(I) + 0.5) * DH)

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
double r[CELLS_COUNT];
double u[CELLS_COUNT];
double v[CELLS_COUNT];
double w[CELLS_COUNT];
double p[CELLS_COUNT];
//
double ru[CELLS_COUNT];
double rv[CELLS_COUNT];
double rw[CELLS_COUNT];
double E[CELLS_COUNT];
//
double fp_r[CELLS_COUNT];
double fp_ru[CELLS_COUNT];
double fp_rv[CELLS_COUNT];
double fp_rw[CELLS_COUNT];
double fp_E[CELLS_COUNT];
double fn_r[CELLS_COUNT];
double fn_ru[CELLS_COUNT];
double fn_rv[CELLS_COUNT];
double fn_rw[CELLS_COUNT];
double fn_E[CELLS_COUNT];
//
double gp_r[CELLS_COUNT];
double gp_ru[CELLS_COUNT];
double gp_rv[CELLS_COUNT];
double gp_rw[CELLS_COUNT];
double gp_E[CELLS_COUNT];
double gn_r[CELLS_COUNT];
double gn_ru[CELLS_COUNT];
double gn_rv[CELLS_COUNT];
double gn_rw[CELLS_COUNT];
double gn_E[CELLS_COUNT];
//
double hp_r[CELLS_COUNT];
double hp_ru[CELLS_COUNT];
double hp_rv[CELLS_COUNT];
double hp_rw[CELLS_COUNT];
double hp_E[CELLS_COUNT];
double hn_r[CELLS_COUNT];
double hn_ru[CELLS_COUNT];
double hn_rv[CELLS_COUNT];
double hn_rw[CELLS_COUNT];
double hn_E[CELLS_COUNT];

// Координаты ближайшей точки на обтекаемой формы.
double p0_x[CELLS_COUNT];
double p0_y[CELLS_COUNT];
double p0_z[CELLS_COUNT];

// Координаты нормали в точке поверхности.
double p0_normal_x[CELLS_COUNT];
double p0_normal_y[CELLS_COUNT];
double p0_normal_z[CELLS_COUNT];

// Три точки для шаблона.
int t1[CELLS_COUNT];
int t2[CELLS_COUNT];
int t3[CELLS_COUNT];

// Матрицы.
double mat_0p123[CELLS_COUNT][4][4];
double mat_g123[CELLS_COUNT][4][4];
double mat_ee[CELLS_COUNT][4][4];
double vec_1xyz[CELLS_COUNT][4];
double vec_1x0y0z0[CELLS_COUNT][4];
double vec_d[CELLS_COUNT][4];
double vec_base[CELLS_COUNT][4];

// Данные для экспорта.
int export_ids[CELLS_COUNT];
int export_count;

// Вспомогательные данные расчетной области.
// Тип ячейки:
//   0 - обычная расчетная ячейка.
//   1 - граничная ячейка.
//   2 - фиктивная ячейка.
//   3 - внутренняя ячейка, которая не принимает участие в расчетах.
#define KIND_NO     0
#define KIND_COMMON 1
#define KIND_BORDER 2
#define KIND_GHOST  3
#define KIND_INNER  4
int kind[CELLS_COUNT];

// Данные сфер для обтекания.
double sph_x[SPHERES_COUNT];
double sph_y[SPHERES_COUNT];
double sph_z[SPHERES_COUNT];
double sph_r[SPHERES_COUNT];

// Инициализация сферы.
void
sphere_init_xy(int i,
               double sphere_x,
               double sphere_y,
               double sphere_r)
{
    sph_x[i] = sphere_x;
    sph_y[i] = sphere_y;
    sph_z[i] = (static_cast<double>(NZ) / 2.0) * DH;
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

        if ((ix < left_cell_count) && (iy < left_cell_count))
        {
            r[i] = 1.0;
            u[i] = 0.75;
            v[i] = 0.5;
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
        double xl = DH * static_cast<double>(ix);
        double xr = DH * static_cast<double>(ix + 1);
        double yl = DH * static_cast<double>(iy);
        double yr = DH * static_cast<double>(iy + 1);
        double zl = DH * static_cast<double>(iz);
        double zr = DH * static_cast<double>(iz + 1);

        for (int si = 0; si < SPHERES_COUNT; si++)
        {
            double dxl2 = MTH_DIFF2(xl, sph_x[si]);
            double dxr2 = MTH_DIFF2(xr, sph_x[si]);
            double dyl2 = MTH_DIFF2(yl, sph_y[si]);
            double dyr2 = MTH_DIFF2(yr, sph_y[si]);
            double dzl2 = MTH_DIFF2(zl, sph_z[si]);
            double dzr2 = MTH_DIFF2(zr, sph_z[si]);
            double r2 = sph_r[si] * sph_r[si];
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
                double xc = MTH_AVG(xl, xr);
                double yc = MTH_AVG(yl, yr);
                double zc = MTH_AVG(zl, zr);

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

    // Превращаем все граничные ячейки в обычные.
    LOOP1
    {
        if (kind[i] == KIND_BORDER)
        {
            kind[i] = KIND_COMMON;
        }
    }

    // Инициализация данных экспорта

    export_count = 0;

    LOOP1
    {
        if (kind[i] != KIND_INNER)
        {
            export_ids[export_count] = i;
            export_count++;
        }
    }
}

// Вычисление для фиктивных ячеек следующих величин:
// координаты ближайшей точки на сфере,
// нормаль в этой точке.
void
calc_nearest_sphere_points_and_normals()
{
    LOOP3
    {
        int i = LIN(ix, iy, iz);

        if (kind[i] == KIND_GHOST)
        {
            double x = (static_cast<double>(ix) + 0.5) * DH;
            double y = (static_cast<double>(iy) + 0.5) * DH;
            double z = (static_cast<double>(iz) + 0.5) * DH;
            double sx = 0.0;
            double sy = 0.0;
            double sz = 0.0;
            double sr = 0.0;
            double d = 0.0;

            int min_si = 0;
            double min_d = (static_cast<double>(NX + NY + NZ)) * DH; // заведомо большое расстояние

            for (int si = 0; si < SPHERES_COUNT; si++)
            {
                sx = sph_x[si];
                sy = sph_y[si];
                sz = sph_z[si];
                sr = sph_r[si];
                d = dist_to_sphere(x, y, z, sx, sy, sz, sr);

                if (d < min_d)
                {
                    min_si = si;
                    min_d = d;
                }
            }

            // Рассматривается сфера с индексом si.
            sx = sph_x[min_si];
            sy = sph_y[min_si];
            sz = sph_z[min_si];
            sr = sph_r[min_si];
            d = dist_to_sphere(x, y, z, sx, sy, sz, sr);

            double dist_to_center = dist_to_point(x, y, z, sx, sy, sz);

            if (dist_to_center > sr)
            {

#if IBM_DEBUG_PRINT == 1
                // Точка снаружи сферы.
                // Это странно, если фиктивная ячейка находится внутри сферы.
                cout << "Internal error : ghost cell is outside sphere." << endl;
                exit(1);
#endif

            }
            else
            {
                // Точка внутри сферы.

                // Инициализация координат точек на сфере.
                p0_x[i] = x + (sx - x) * (-d / dist_to_center);
                p0_y[i] = y + (sy - y) * (-d / dist_to_center);
                p0_z[i] = z + (sz - z) * (-d / dist_to_center);

                // Инициализация нормалей.
                p0_normal_x[i] = (x - sx) / sr;
                p0_normal_y[i] = (y - sy) / sr;
                p0_normal_z[i] = (z - sz) / sr;
            }
        }
        else
        {
            p0_x[i] = 0.0;
            p0_y[i] = 0.0;
            p0_z[i] = 0.0;

            p0_normal_x[i] = 0.0;
            p0_normal_y[i] = 0.0;
            p0_normal_z[i] = 0.0;
        }
    }
}

// Расширенное получение типа.
int
ext_kind(int ix,
         int iy,
         int iz)
{
    if ((ix < 0) || (ix >= NX) || (iy < 0) || (iy >= NY) || (iz < 0) || (iz >= NZ))
    {
        return KIND_NO;
    }

    return kind[LIN(ix, iy, iz)];
}

// Поиск шаблона.
void
find_templates_points(int ix,
                      int iy,
                      int iz,
                      int *tt[])
{
    // Выполняем полный перебор по трем точкам.

    for (int lx1 = ix - 1; lx1 <= ix + 1; lx1++)
    {
        for (int ly1 = iy - 1; ly1 <= iy + 1; ly1++)
        {
            for (int lz1 = iz - 1; lz1 <= iz + 1; lz1++)
            {
                // Проверка выхода за границы области.
                if ((lx1 < 0) || (lx1 >= NX) || (ly1 < 0) || (ly1 >= NY) || (lz1 < 0) || (lz1 >= NZ))
                {
                    continue;
                }

                int bi1 = LIN(lx1, ly1, lz1);

                // Проверка типа ячейки.
                if (kind[bi1] != KIND_COMMON)
                {
                    continue;
                }

                for (int lx2 = ix - 1; lx2 <= ix + 1; lx2++)
                {
                    for (int ly2 = iy - 1; ly2 <= iy + 1; ly2++)
                    {
                        for (int lz2 = iz - 1; lz2 <= iz + 1; lz2++)
                        {
                            // Проверка выхода за границы области.
                            if ((lx2 < 0) || (lx2 >= NX) || (ly2 < 0) || (ly2 >= NY) || (lz2 < 0) || (lz2 >= NZ))
                            {
                                continue;
                            }

                            int bi2 = LIN(lx2, ly2, lz2);

                            // Проверка типа ячейки.
                            if (kind[bi2] != KIND_COMMON)
                            {
                                continue;
                            }

                            for (int lx3 = ix - 1; lx3 <= ix + 1; lx3++)
                            {
                                for (int ly3 = iy - 1; ly3 <= iy + 1; ly3++)
                                {
                                    for (int lz3 = iz - 1; lz3 <= iz + 1; lz3++)
                                    {
                                        // Проверка выхода за границы области.
                                        if ((lx3 < 0) || (lx3 >= NX) || (ly3 < 0) || (ly3 >= NY) || (lz3 < 0) || (lz3 >= NZ))
                                        {
                                            continue;
                                        }

                                        int bi3 = LIN(lx3, ly3, lz3);

                                        // Проверка типа ячейки.
                                        if (kind[bi3] != KIND_COMMON)
                                        {
                                            continue;
                                        }

                                        // Проверяем, лежат ли точки на одной прямой.
                                        if ((lx3 - lx2 == lx2 - lx1) && (ly3 - ly2 == ly2 - ly1) && (lz3 - lz2 == lz2 - lz1))
                                        {
                                            continue;
                                        }

                                        // Проверяем компланарность векторов из фиктивной ячейки к ячейкам шаблонов.
                                        int ax = lx1 - ix;
                                        int ay = ly1 - iy;
                                        int az = lz1 - iz;
                                        int bx = lx2 - ix;
                                        int by = ly2 - iy;
                                        int bz = lz2 - iz;
                                        int cx = lx3 - ix;
                                        int cy = ly3 - iy;
                                        int cz = lz3 - iz;
                                        int det = ax * (by * cz - bz * cy)
                                                  - ay * (bx * cz - bz * cx)
                                                  + az * (bx * cy - by * cx);

                                        if (det == 0)
                                        {
                                            continue;
                                        }

                                        // Все проверки выполнены, можно сохранять точки.
                                        *tt[0] = bi1;
                                        *tt[1] = bi2;
                                        *tt[2] = bi3;

                                        return;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

#if IBM_DEBUG_PRINT == 1
    // Если все перебрали и не нашли шаблов, то выходим с ошибкой.
    cout << "Error : find_templates_points : template is not found." << endl;
    cout << "ix = " << ix << ", iy = " << iy << ", iz = " << iz << endl;
    cout << "Neighborhood:" << endl;
    cout << ext_kind(ix - 1, iy + 1, iz - 1) << " " << ext_kind(ix, iy + 1, iz - 1) << " " << ext_kind(ix + 1, iy + 1, iz - 1) << endl;
    cout << ext_kind(ix - 1, iy, iz - 1) << " " << ext_kind(ix, iy, iz - 1) << " " << ext_kind(ix + 1, iy, iz - 1) << endl;
    cout << ext_kind(ix - 1, iy - 1, iz - 1) << " " << ext_kind(ix, iy - 1, iz - 1) << " " << ext_kind(ix + 1, iy - 1, iz - 1) << endl;
    cout << "------------------" << endl;
    cout << "       " << ext_kind(ix - 1, iy + 1, iz) << " " << ext_kind(ix, iy + 1, iz) << " " << ext_kind(ix + 1, iy + 1, iz) << endl;
    cout << "       " << ext_kind(ix - 1, iy, iz) << " " << ext_kind(ix, iy, iz) << " " << ext_kind(ix + 1, iy, iz) << endl;
    cout << "       " << ext_kind(ix - 1, iy - 1, iz) << " " << ext_kind(ix, iy - 1, iz) << " " << ext_kind(ix + 1, iy - 1, iz) << endl;
    cout << "------------------" << endl;
    cout << "              " << ext_kind(ix - 1, iy + 1, iz + 1) << " " << ext_kind(ix, iy + 1, iz + 1) << " " << ext_kind(ix + 1, iy + 1, iz + 1) << endl;
    cout << "              " << ext_kind(ix - 1, iy, iz + 1) << " " << ext_kind(ix, iy, iz + 1) << " " << ext_kind(ix + 1, iy, iz + 1) << endl;
    cout << "              " << ext_kind(ix - 1, iy - 1, iz + 1) << " " << ext_kind(ix, iy - 1, iz + 1) << " " << ext_kind(ix + 1, iy - 1, iz + 1) << endl;
    exit(1);
#endif

}

// Определение расчетных шаблонов.
void
define_templates()
{
    LOOP3
    {
        int i = LIN(ix, iy, iz);

        if (kind[i] == KIND_GHOST)
        {
            int *tt[] = { &t1[i], &t2[i], &t3[i] };

            find_templates_points(ix, iy, iz, tt);
        }
        else
        {
            t1[i] = 0;
            t2[i] = 0;
            t3[i] = 0;
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
    f << "VARIABLES=\"X\", \"Y\", \"Z\", \"Id\", "
      << "\"Rho\", \"U\", \"V\", \"W\", \"P\", "
      << "\"Kind\", \"P0X\", \"P0Y\", \"P0Z\", \"P0NormalX\", \"P0NormalY\", \"P0NormalZ\", "
      << "\"T1\", \"T2\", \"T3\"" << endl;
    f << "ZONE T=\"single zone\"" << endl;
    f << "NODES=" << (8 * export_count) << endl;
    f << "ELEMENTS=" << export_count << endl;
    f << "DATAPACKING=BLOCK" << endl;
    f << "ZONETYPE=FEBRICK" << endl;
    f << "VARLOCATION=([4-19]=CELLCENTERED)" << endl;

#define LX (ix * DH)
#define HX ((ix + 1) * DH)
#define LY (iy * DH)
#define HY ((iy + 1) * DH)
#define LZ (iz * DH)
#define HZ ((iz + 1) * DH)

    LOOP3
        if (kind[LIN(ix, iy, iz)] != KIND_INNER)
            f << LX << " " << HX << " " << LX << " " << HX << " " << LX << " " << HX << " " << LX << " " << HX << " ";
    f << endl;

    LOOP3
        if (kind[LIN(ix, iy, iz)] != KIND_INNER)
            f << LY << " " << LY << " " << HY << " " << HY << " " << LY << " " << LY << " " << HY << " " << HY << " ";
    f << endl;

    LOOP3
        if (kind[LIN(ix, iy, iz)] != KIND_INNER)
            f << LZ << " " << LZ << " " << LZ << " " << LZ << " " << HZ << " " << HZ << " " << HZ << " " << HZ << " ";
    f << endl;

    for(int i = 0; i < export_count; i++) f << i << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << r[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << u[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << v[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << w[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << p[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << kind[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << p0_x[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << p0_y[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << p0_z[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << p0_normal_x[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << p0_normal_y[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << p0_normal_z[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << t1[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << t2[export_ids[i]] << " "; f << endl;
    for(int i = 0; i < export_count; i++) f << t3[export_ids[i]] << " "; f << endl;

#undef LX
#undef HX
#undef LY
#undef HY
#undef LZ
#undef HZ

    for (int i = 0; i < export_count; i++)
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
        // В потоках участвуют обычные и фиктивные ячейки.
        if ((kind[i] == KIND_COMMON) || (kind[i] == KIND_GHOST))
        {
            ru[i] = r[i] * u[i];
            rv[i] = r[i] * v[i];
            rw[i] = r[i] * w[i];

            // E = rho * (V^2/2 + e) = rho * (V^2/2 + p/((GAMMA - 1.0) * rho))
            //   = rho * V^2/2 + p/(GAMMA - 1)
            E[i] = 0.5 * r[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]) + p[i] / (GAMMA - 1.0);
        }
    }
}

// Перевод консервативных величин в примитивные.
void
u_to_d()
{
    LOOP1
    {
        // В основные величины пересчитываем только обычные ячейки.
        if (kind[i] == KIND_COMMON)
        {
            u[i] = ru[i] / r[i];
            v[i] = rv[i] / r[i];
            w[i] = rw[i] / r[i];

            // E = rho * V^2/2 + p/(GAMMA - 1.0)
            // p = (E - rho * V^2/2) * (GAMMA - 1.0)
            p[i] = (E[i] - 0.5 * r[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i])) * (GAMMA - 1.0);

#if IBM_DEBUG_PRINT == 1
            // Плотность, давление и энергия это положительные величины.
            if ((r[i] <= 0.0) || (p[i] <= 0.0) || (E[i] <= 0.0))
            {
                cout << "Error : u_to_d : negative density, pressure of energy." << endl;
                cout << "r / p / E = " << r[i] << " / " << p[i] << " / " << E[i] << endl;
                exit(1);
            }

            if (isnan(r[i]) || isnan(u[i]) || isnan(v[i]) || isnan(w[i]) || isnan(p[i]))
            {
                cout << "Error : u_to_d : not a number." << endl;
                cout << "r / ru / rv / rw / E = " << r[i] << " / " << ru[i] << " / " << rv[i] << " / " << rw[i] << " / " << E[i] << endl;
                cout << "r / u / v / w / p = " << r[i] << " / " << u[i] << " / " << v[i] << " / " << w[i] << " / " << p[i] << endl;
                exit(1);
            }
#endif

        }
    }
}

// Вычисление f pos/neg.
void
calc_f()
{
    LOOP1
    {
        if ((kind[i] == KIND_COMMON) || (kind[i] == KIND_GHOST))
        {
            double a = sqrt(GAMMA * p[i] / r[i]);
            double l1 = u[i] - a;
            double l2 = u[i];
            double l5 = u[i] + a;
            double lp1 = 0.5 * (l1 + abs(l1));
            double lp2 = 0.5 * (l2 + abs(l2));
            double lp5 = 0.5 * (l5 + abs(l5));
            double ln1 = 0.5 * (l1 - abs(l1));
            double ln2 = 0.5 * (l2 - abs(l2));
            double ln5 = 0.5 * (l5 - abs(l5));
            double k = 0.5 * r[i] / GAMMA;
            double V2 = u[i] * u[i] + v[i] * v[i] + w[i] * w[i];
            double H = 0.5 * V2 + a * a / (GAMMA - 1.0);

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

#if IBM_DEBUG_PRINT == 1
            if (isnan(fp_r[i]) || isnan(fn_r[i]))
            {
                cout << "Error : calc_f : fp_r / fn_r is not a number." << endl;
                cout << "fp_r / fn_r = " << fp_r[i] << ", " << fn_r[i] << endl;
                exit(1);
            }
#endif

        }
    }
}

// Вычисление g pos/neg.
void
calc_g()
{
    LOOP1
    {
        if ((kind[i] == KIND_COMMON) || (kind[i] == KIND_GHOST))
        {
            double a = sqrt(GAMMA * p[i] / r[i]);
            double l1 = v[i] - a;
            double l2 = v[i];
            double l5 = v[i] + a;
            double lp1 = 0.5 * (l1 + abs(l1));
            double lp2 = 0.5 * (l2 + abs(l2));
            double lp5 = 0.5 * (l5 + abs(l5));
            double ln1 = 0.5 * (l1 - abs(l1));
            double ln2 = 0.5 * (l2 - abs(l2));
            double ln5 = 0.5 * (l5 - abs(l5));
            double k = 0.5 * r[i] / GAMMA;
            double V2 = u[i] * u[i] + v[i] * v[i] + w[i] * w[i];
            double H = 0.5 * V2 + a * a / (GAMMA - 1.0);

            gp_r[i] = k * (lp1 + 2.0 * (GAMMA - 1.0) * lp2 + lp5);
            gp_ru[i] = k * (u[i] * lp1 + 2.0 * (GAMMA - 1.0) * u[i] * lp2 + u[i] * lp5);
            gp_rv[i] = k * ((v[i] - a) * lp1 + 2.0 * (GAMMA - 1.0) * v[i] * lp2 + (v[i] + a) * lp5);
            gp_rw[i] = k * (w[i] * lp1 + 2.0 * (GAMMA - 1.0) * w[i] * lp2 + w[i] * lp5);
            gp_E[i] = k * ((H - v[i] * a) * lp1 + (GAMMA - 1.0) * V2 * lp2 + (H + v[i] * a) * lp5);

            gn_r[i] = k * (ln1 + 2.0 * (GAMMA - 1.0) * ln2 + ln5);
            gn_ru[i] = k * (u[i] * ln1 + 2.0 * (GAMMA - 1.0) * u[i] * ln2 + u[i] * ln5);
            gn_rv[i] = k * ((v[i] - a) * ln1 + 2.0 * (GAMMA - 1.0) * v[i] * ln2 + (v[i] + a) * ln5);
            gn_rw[i] = k * (w[i] * ln1 + 2.0 * (GAMMA - 1.0) * w[i] * ln2 + w[i] * ln5);
            gn_E[i] = k * ((H - v[i] * a) * ln1 + (GAMMA - 1.0) * V2 * ln2 + (H + v[i] * a) * ln5);

#if IBM_DEBUG_PRINT == 1
            if (isnan(gp_r[i]) || isnan(gn_r[i]))
            {
                cout << "Error : calc_g : gp_r / gn_r is not a number." << endl;
                cout << "gp_r / gn_r = " << gp_r[i] << ", " << gn_r[i] << endl;
                exit(1);
            }
#endif

        }
    }
}

// Вычисление h pos/neg.
void
calc_h()
{
    LOOP1
    {
        if ((kind[i] == KIND_COMMON) || (kind[i] == KIND_GHOST))
        {
            double a = sqrt(GAMMA * p[i] / r[i]);
            double l1 = w[i] - a;
            double l2 = w[i];
            double l5 = w[i] + a;
            double lp1 = 0.5 * (l1 + abs(l1));
            double lp2 = 0.5 * (l2 + abs(l2));
            double lp5 = 0.5 * (l5 + abs(l5));
            double ln1 = 0.5 * (l1 - abs(l1));
            double ln2 = 0.5 * (l2 - abs(l2));
            double ln5 = 0.5 * (l5 - abs(l5));
            double k = 0.5 * r[i] / GAMMA;
            double V2 = u[i] * u[i] + v[i] * v[i] + w[i] * w[i];
            double H = 0.5 * V2 + a * a / (GAMMA - 1.0);

            hp_r[i] = k * (lp1 + 2.0 * (GAMMA - 1.0) * lp2 + lp5);
            hp_ru[i] = k * (u[i] * lp1 + 2.0 * (GAMMA - 1.0) * u[i] * lp2 + u[i] * lp5);
            hp_rv[i] = k * (v[i] * lp1 + 2.0 * (GAMMA - 1.0) * v[i] * lp2 + v[i] * lp5);
            hp_rw[i] = k * ((w[i] - a) * lp1 + 2.0 * (GAMMA - 1.0) * w[i] * lp2 + (w[i] + a) * lp5);
            hp_E[i] = k * ((H - w[i] * a) * lp1 + (GAMMA - 1.0) * V2 * lp2 + (H + w[i] * a) * lp5);

            hn_r[i] = k * (ln1 + 2.0 * (GAMMA - 1.0) * ln2 + ln5);
            hn_ru[i] = k * (u[i] * ln1 + 2.0 * (GAMMA - 1.0) * u[i] * ln2 + u[i] * ln5);
            hn_rv[i] = k * (v[i] * ln1 + 2.0 * (GAMMA - 1.0) * v[i] * ln2 + v[i] * ln5);
            hn_rw[i] = k * ((w[i] - a) * ln1 + 2.0 * (GAMMA - 1.0) * w[i] * ln2 + (w[i] + a) * ln5);
            hn_E[i] = k * ((H - w[i] * a) * ln1 + (GAMMA - 1.0) * V2 * ln2 + (H + w[i] * a) * ln5);

#if IBM_DEBUG_PRINT == 1
            if (isnan(hp_r[i]) || isnan(hn_r[i]))
            {
                cout << "Error : calc_h : hp_r / hn_r is not a number." << endl;
                cout << "hp_r / hn_r = " << hp_r[i] << ", " << hn_r[i] << endl;
                exit(1);
            }
#endif

        }
    }
}

// Вычисление f/g/h pos/neg.
void
calc_fgh()
{
    calc_f();
    calc_g();
    calc_h();
}

// Вычисление потоков.
void
calc_flows()
{
    calc_fgh();

    LOOP3
    {
        int i = LIN(ix, iy, iz);

        // Вычисления производятся только для COMMON ячеек.
        if (kind[i] == KIND_COMMON)
        {
            //
            // Направление X.
            //

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

            //
            // Направление Y.
            //

            li = i;
            ri = i;

            if (iy > 0)
            {
                li = LIN(ix, iy - 1, iz);
            }

            if (iy < NY - 1)
            {
                ri = LIN(ix, iy + 1, iz);
            }

            r[i] -= (DT / DH) * (gp_r[i] + gn_r[ri] - gp_r[li] - gn_r[i]);
            ru[i] -= (DT / DH) * (gp_ru[i] + gn_ru[ri] - gp_ru[li] - gn_ru[i]);
            rv[i] -= (DT / DH) * (gp_rv[i] + gn_rv[ri] - gp_rv[li] - gn_rv[i]);
            rw[i] -= (DT / DH) * (gp_rw[i] + gn_rw[ri] - gp_rw[li] - gn_rw[i]);
            E[i] -= (DT / DH) * (gp_E[i] + gn_E[ri] - gp_E[li] - gn_E[i]);

            //
            // Направление Z.
            //

            li = i;
            ri = i;

            if (iz > 0)
            {
                li = LIN(ix, iy, iz - 1);
            }

            if (iz < NZ - 1)
            {
                ri = LIN(ix, iy, iz + 1);
            }

            r[i] -= (DT / DH) * (hp_r[i] + hn_r[ri] - hp_r[li] - hn_r[i]);
            ru[i] -= (DT / DH) * (hp_ru[i] + hn_ru[ri] - hp_ru[li] - hn_ru[i]);
            rv[i] -= (DT / DH) * (hp_rv[i] + hn_rv[ri] - hp_rv[li] - hn_rv[i]);
            rw[i] -= (DT / DH) * (hp_rw[i] + hn_rw[ri] - hp_rw[li] - hn_rw[i]);
            E[i] -= (DT / DH) * (hp_E[i] + hn_E[ri] - hp_E[li] - hn_E[i]);

#if IBM_DEBUG_PRINT == 1
            if (isnan(r[i]))
            {
                cout << "Error : calc_flows : density is not a number.";
                cout << "r = " << r[i] << endl;
                exit(1);
            }
#endif

        }
    }
}

// Подготовка к аппроксимации значений в фиктивных ячейках.
void
pre_approximate_values()
{
    LOOP3
    {
        int i = LIN(ix, iy, iz);

        if (kind[i] == KIND_GHOST)
        {
            // Требуется выполнить следующую аппроксимацию:
            // 1. Плотность - скалярная величина.
            // 2. Давление - скалярная величина.
            // 3. Скорость - векторная величина.

            int tmpl1 = t1[i];
            int tmpl2 = t2[i];
            int tmpl3 = t3[i];
            int tmpl1x = UNLINX(tmpl1);
            int tmpl1y = UNLINY(tmpl1);
            int tmpl1z = UNLINZ(tmpl1);
            int tmpl2x = UNLINX(tmpl2);
            int tmpl2y = UNLINY(tmpl2);
            int tmpl2z = UNLINZ(tmpl2);
            int tmpl3x = UNLINX(tmpl3);
            int tmpl3y = UNLINY(tmpl3);
            int tmpl3z = UNLINZ(tmpl3);
            double x1 = CCORD(tmpl1x);
            double y1 = CCORD(tmpl1y);
            double z1 = CCORD(tmpl1z);
            double x2 = CCORD(tmpl2x);
            double y2 = CCORD(tmpl2y);
            double z2 = CCORD(tmpl2z);
            double x3 = CCORD(tmpl3x);
            double y3 = CCORD(tmpl3y);
            double z3 = CCORD(tmpl3z);

            m4x4_init_vec(vec_1xyz[i],
                          1.0, CCORD(ix), CCORD(iy), CCORD(iz));

            m4x4_init_mat(mat_0p123[i],
                          0.0, p0_normal_x[i], p0_normal_y[i], p0_normal_z[i],
                          1.0, x1, y1, z1,
                          1.0, x2, y2, z2,
                          1.0, x3, y3, z3);

            if (!m4x4_invert(mat_0p123[i], mat_0p123[i]))
            {

#if IBM_DEBUG_PRINT == 1
                cout << "Error : can not invert matrix B<0'123>." << endl;
                m4x4_print(mat_0p123[i]);
                exit(1);
#endif

            }

            m4x4_init_mat(mat_g123[i],
                          1.0, CCORD(ix), CCORD(iy), CCORD(iz),
                          1.0, x1, y1, z1,
                          1.0, x2, y2, z2,
                          1.0, x3, y3, z3);

            if (!m4x4_invert(mat_g123[i], mat_g123[i]))
            {

#if IBM_DEBUG_PRINT == 1
                cout << "Error : can not invert matrix B<G123>." << endl;
                m4x4_print(mat_g123[i]);
                exit(1);
#endif

            }

            m4x4_init_vec(vec_1x0y0z0[i],
                          1.0, p0_x[i], p0_y[i], p0_z[i]);
            m4x4_mul_vec_mat(vec_1x0y0z0[i], mat_g123[i], vec_d[i]);
            m4x4_mul_vec_mat(vec_1xyz[i], mat_0p123[i], vec_base[i]);

            if ((abs(p0_normal_x[i]) <= abs(p0_normal_y[i]))
                && (abs(p0_normal_x[i]) <= abs(p0_normal_z[i])))
            {
                m4x4_init_mat(mat_ee[i],
                              p0_normal_x[i], p0_normal_y[i], p0_normal_z[i], 0.0,
                              -p0_normal_y[i], p0_normal_x[i], 0.0, 0.0,
                              0.0, -p0_normal_z[i], p0_normal_y[i], 0.0,
                              0.0, 0.0, 0.0, 1.0);
            }
            else
            {
                m4x4_init_mat(mat_ee[i],
                              p0_normal_x[i], p0_normal_y[i], p0_normal_z[i], 0.0,
                              -p0_normal_y[i], p0_normal_x[i], 0.0, 0.0,
                              -p0_normal_z[i], 0.0, p0_normal_x[i], 0.0,
                              0.0, 0.0, 0.0, 1.0);
            }

            if (!m4x4_invert(mat_ee[i], mat_ee[i]))
            {

#if IBM_DEBUG_PRINT == 1
                cout << "Error : can not invert matrix EE." << endl;
                m4x4_print(mat_ee[i]);
                exit(1);
#endif

            }
        }
    }
}

// Аппроксимация значений в фиктивных ячейках.
void
approximate_values()
{
    LOOP3
    {
        int i = LIN(ix, iy, iz);

        if (kind[i] == KIND_GHOST)
        {
            // Требуется выполнить следующую аппроксимацию:
            // 1. Плотность - скалярная величина.
            // 2. Давление - скалярная величина.
            // 3. Скорость - векторная величина.

            int tmpl1 = t1[i];
            int tmpl2 = t2[i];
            int tmpl3 = t3[i];
            int tmpl1x = UNLINX(tmpl1);
            int tmpl1y = UNLINY(tmpl1);
            int tmpl1z = UNLINZ(tmpl1);
            int tmpl2x = UNLINX(tmpl2);
            int tmpl2y = UNLINY(tmpl2);
            int tmpl2z = UNLINZ(tmpl2);
            int tmpl3x = UNLINX(tmpl3);
            int tmpl3y = UNLINY(tmpl3);
            int tmpl3z = UNLINZ(tmpl3);
            double vec_phi[4];
            double vec_a[4];

            //
            // Аппроксимация плотности.
            //

            m4x4_init_vec(vec_phi,
                          0.0, r[tmpl1], r[tmpl2], r[tmpl3]);
            m4x4_mul_mat_vec(mat_0p123[i], vec_phi, vec_a);
            r[i] = m4x4_scalar_product(vec_a, vec_1xyz[i]);

#if IBM_DEBUG_PRINT == 1
            if (isnan(r[i]) || (r[i] <= 0.0))
            {
                cout << "Error : approximate_values : wrong value of r." << endl;
                cout << "r = " << r[i] << endl;
                exit(0);
            }
#endif

            //
            // Аппроксимация давления.
            // Матрица B та же, надо только поменять phi.
            //

            m4x4_init_vec(vec_phi,
                          0.0, p[tmpl1], p[tmpl2], p[tmpl3]);
            m4x4_mul_mat_vec(mat_0p123[i], vec_phi, vec_a);
            p[i] = m4x4_scalar_product(vec_a, vec_1xyz[i]);

#if IBM_DEBUG_PRINT == 1
            if (isnan(p[i]) || (p[i] <= 0.0))
            {
                cout << "Error : approximate_values : wrong value of p." << endl;
                cout << "p = " << p[i] << endl;
                cout << "tmpl ps : " << p[tmpl1] << ", " << p[tmpl2] << ", " << p[tmpl3] << endl;
                exit(0);
            }
#endif

            //
            // Аппроксимация скорости.
            //

            double q;
            double t_xy;
            double t_xz;
            double t_yz;
            double vec_t_xy[4];
            double vec_t_xz[4];
            double vec_t_yz[4];
            double vec_q_ts[4];
            double vec_vg[4];

            // Явное вычисление q.
            q = - (vec_d[i][1] * (u[tmpl1] * p0_normal_x[i] + v[tmpl1] * p0_normal_y[i] + w[tmpl1] * p0_normal_z[i])
                   + vec_d[i][2] * (u[tmpl2] * p0_normal_x[i] + v[tmpl2] * p0_normal_y[i] + w[tmpl2] * p0_normal_z[i])
                   + vec_d[i][3] * (u[tmpl3] * p0_normal_x[i] + v[tmpl3] * p0_normal_y[i] + w[tmpl3] * p0_normal_z[i]))
                  / (vec_d[i][0]);

#if IBM_DEBUG_PRINT == 1
            if (isnan(q))
            {
                cout << "Error : Q is not a number." << endl;
                cout << "normal : " << p0_normal_x[i] << ", " << p0_normal_y[i] << ", " << p0_normal_z[i] << endl;
                cout << "tmpl1 velocity : " << u[tmpl1] << ", " << v[tmpl1] << ", " << w[tmpl1] << endl;
                cout << "tmpl2 velocity : " << u[tmpl2] << ", " << v[tmpl2] << ", " << w[tmpl2] << endl;
                cout << "tmpl3 velocity : " << u[tmpl3] << ", " << v[tmpl3] << ", " << w[tmpl3] << endl;
                cout << "vec_d :" << endl;
                m4x4_print_vec(vec_d[i]);
                cout << "Q = " << q << endl;
                exit(1);
            }
#endif

            // Явное вычисление t_xy, t_xz, t_yz.
            m4x4_init_vec(vec_t_xy,
                          0.0,
                          -u[tmpl1] * p0_normal_y[i] + v[tmpl1] * p0_normal_x[i],
                          -u[tmpl2] * p0_normal_y[i] + v[tmpl2] * p0_normal_x[i],
                          -u[tmpl3] * p0_normal_y[i] + v[tmpl3] * p0_normal_x[i]);
            m4x4_init_vec(vec_t_xz,
                          0.0,
                          -u[tmpl1] * p0_normal_z[i] + w[tmpl1] * p0_normal_x[i],
                          -u[tmpl2] * p0_normal_z[i] + w[tmpl2] * p0_normal_x[i],
                          -u[tmpl3] * p0_normal_z[i] + w[tmpl3] * p0_normal_x[i]);
            m4x4_init_vec(vec_t_yz,
                          0.0,
                          -v[tmpl1] * p0_normal_z[i] + w[tmpl1] * p0_normal_y[i],
                          -v[tmpl2] * p0_normal_z[i] + w[tmpl2] * p0_normal_y[i],
                          -v[tmpl3] * p0_normal_z[i] + w[tmpl3] * p0_normal_y[i]);
            t_xy = m4x4_scalar_product(vec_base[i], vec_t_xy);
            t_xz = m4x4_scalar_product(vec_base[i], vec_t_xz);
            t_yz = m4x4_scalar_product(vec_base[i], vec_t_yz);

            if ((abs(p0_normal_x[i]) <= abs(p0_normal_y[i]))
                && (abs(p0_normal_x[i]) <= abs(p0_normal_z[i])))
            {
                m4x4_init_vec(vec_q_ts,
                              q, t_xy, t_yz, 1.0);
            }
            else
            {
                m4x4_init_vec(vec_q_ts,
                              q, t_xy, t_xz, 1.0);
            }

            m4x4_mul_mat_vec(mat_ee[i], vec_q_ts, vec_vg);

#if IBM_DEBUG_PRINT == 1
            if (isnan(vec_vg[0]) || isnan(vec_vg[1]) || isnan(vec_vg[2]))
            {
                cout << "Error : not a number in one or more components of velocity vector." << endl;
                cout << "mat_ee_inv : " << endl;
                m4x4_print(mat_ee[i]);
                cout << "vec_ts_q :" << endl;
                m4x4_print_vec(vec_q_ts);
                cout << "vec_vg :" << endl;
                m4x4_print_vec(vec_vg);
                exit(1);
            }
#endif

            u[i] = vec_vg[0];
            v[i] = vec_vg[1];
            w[i] = vec_vg[2];
        }
    }
}

// Шаг вычислений.
void
step()
{
    approximate_values();
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
    calc_nearest_sphere_points_and_normals();
    define_templates();
    pre_approximate_values();
    calc_area_paraview_export(0);

    for (int i = 0; i < TIME_STEPS; i++)
    {
        if ((i + 1) % EXPORT_DISCRETION == 0)
        {
            cout << ".... step " << (i + 1) << " of " << TIME_STEPS << endl;
        }

        step();

        if ((i + 1) % EXPORT_DISCRETION == 0)
        {
            calc_area_paraview_export(i + 1);
        }
    }
}
