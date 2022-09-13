#include <fstream>

#include "conf.h"

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

// Данные расчетной области.
// Примитивные данные:
//   - плотность,
//   - составляющие скорости по x, y, z,
//   - давление.
// Сохраняющиеся данные:
//   - плотность импульса
//   - плотность полной энергии
float rho[CELLS_COUNT];
float u[CELLS_COUNT];
float v[CELLS_COUNT];
float w[CELLS_COUNT];
float p[CELLS_COUNT];
//
float rho_u[CELLS_COUNT];
float rho_v[CELLS_COUNT];
float rho_w[CELLS_COUNT];
float E[CELLS_COUNT];

// Инициализация расчетной области.
void
calc_area_init()
{
    for (int i = 0; i < CELLS_COUNT; i++)
    {
        rho[i] = 0.5;
        u[i] = 1.5;
        v[i] = 2.5;
        w[i] = 3.5;
        p[i] = 4.5;
        rho_u[i] = 5.5;
        rho_v[i] = 6.5;
        rho_w[i] = 7.5;
        E[i] = 8.5;
    }
}

// Экспорт в ParaView.
void
calc_area_paraview_export()
{
    ofstream f("export.dat");

    f << "TITLE=\"[" << NX << " * " << NY << " * " << NZ << "] calc area\"" << endl;
    f << "VARIABLES=\"X\", \"Y\", \"Z\", \"Rho\"" << endl;
    f << "ZONE T=\"single zone\"" << endl;
    f << "NODES=" << (8 * CELLS_COUNT) << endl;
    f << "ELEMENTS=" << CELLS_COUNT << endl;
    f << "DATAPACKING=BLOCK" << endl;
    f << "ZONETYPE=FEBRICK" << endl;
    f << "VARLOCATION=([4-4]=CELLCENTERED)" << endl;

#define LOOP for (int ix = 0; ix < NX; ix++) for (int iy = 0; iy < NY; iy++) for (int iz = 0; iz < NZ; iz++)
#define LX (ix * DH)
#define HX ((ix + 1) * DH)
#define LY (iy * DH)
#define HY ((iy + 1) * DH)
#define LZ (iz * DH)
#define HZ ((iz + 1) * DH)

    LOOP f << LX << " " << HX << " " << LX << " " << HX << " " << LX << " " << HX << " " << LX << " " << HX << " "; f << endl;
    LOOP f << LY << " " << LY << " " << HY << " " << HY << " " << LY << " " << LY << " " << HY << " " << HY << " "; f << endl;
    LOOP f << LZ << " " << LZ << " " << LZ << " " << LZ << " " << HZ << " " << HZ << " " << HZ << " " << HZ << " "; f << endl;
    LOOP f << rho[LIN(ix, iy, iz)] << " "; f << endl;

#undef LOOP
#undef LX
#undef HX
#undef LY
#undef HY
#undef LZ
#undef HZ

    for (int i = 0; i < CELLS_COUNT; i++)
    {
        int s = 8 * i;

        // Нумерация точек в элементе начинается с единицы.
        f << (s + 1) << " " << (s + 2) << " " << (s + 4) << " " << (s + 3) << " "
          << (s + 5) << " " << (s + 6) << " " << (s + 8) << " " << (s + 7) << endl;
    }

    f.close();
}

// Точка входа.
int
main()
{
    calc_area_init();
    calc_area_paraview_export();
}
