// Настройки.

#ifndef CONF_H
#define CONF_H

// Размеры расчетной области.
// Область 10 * 5.
#define NX 112
#define NY 50
#define NZ 2

// Шаг по пространству и по времени.
#define DH 0.1  // м
#define DT 0.01 // с

// Количество шагов по времени.
#define TIME_STEPS 500

// Показатель адиабаты.
#define GAMMA 1.4

// Количество шаров для обтекания.
#define SPHERES_COUNT 5

// Дискретность экспорта.
#define EXPORT_DISCRETION 500

// Опции отладки.
#define MTH_DEBUG_PRINT 0
#define IBM_DEBUG_PRINT 0
#define INTEL_RUN 0
#define USE_AVX512 0
#define NONPROB_BRANCH_LOCALIZATION 0

#endif
