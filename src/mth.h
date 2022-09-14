// Математические примитивы.

#ifndef MTH_H
#define MTH_H

// Среднее значение.
#define MTH_AVG(A, B) (0.5 * ((A) + (B)))

// Квадрат разницы между значениями.
#define MTH_DIFF2(A, B) (((A) - (B)) * ((A) - (B)))

// Квадрат рассмояния между точками.
#define MTH_DIST2(X1, Y1, Z1, X2, Y2, Z2) (MTH_DIFF2(X1, X2) + MTH_DIFF2(Y1, Y2) + MTH_DIFF2(Z1, Z2))

#endif
