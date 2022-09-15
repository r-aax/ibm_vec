// Настройки.

#ifndef CONF_H
#define CONF_H

// Размеры расчетной области.
// Область 10 * 5.
#define NX 100
#define NY 50
#define NZ 1

// Шаг по пространству и по времени.
#define DH 0.1  // м
#define DT 0.01 // с

// Количество шагов по времени.
#define TIME_STEPS 500

// Показатель адиабаты.
#define GAMMA 1.4

// Количество шаров для обтекания.
#define SPHERES_COUNT 5

#endif
