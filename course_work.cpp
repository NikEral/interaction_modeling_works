/*
Курсовое задание
Задание:
    1.  Создать класс нуклона, который будет хранить информацию о заряде, массе и положении в декартовой системе координат
    2.  Создать класс ядра, который содержит в себе набор нуклонов, заряд, массу, радиус(R), и другие параметры ядра:
         a. r0
         b. V0
         c. a
         d. beta2
         e. beta4
    3.  Случайным образом согласно потенциалу Вудса-Саксона(по варианту) сгенерировать нуклоны ядра.
    4.  Пересечь 2 таких ядра(каждый генерируется отдельно!) с расстоянием между их центрами случайное от 0 до 2R
    5.  Посчитать число нуклонов обоих ядер суммарно, попавших в область пересечения
    6.  Повторить процесс 10 000 раз
    7.  Построить график распределения получившегося числа нуклонов в области перекрытия(двух ядер по отдельности и суммарно)
    8.  Профитировать функцией, которая опишет данное распределение
    9.  Нарисовать график с функцией фитирования, на рисунке также должна быть написана использованная функция.
    10. Сохранить гистограмму в root файл, канвас в png


    №   Ядро    A   Z   N   V₀ (МэВ)    r₀ (фм) a (фм)  β₂      β₄    Тип деформации
    7   ⁵⁶Fe    56  26  30  51.0        1.25    0.535   0.00    0.00    Сферическое
*/
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// #define A 56
// #define Z 26
// #define N 30
// #define V0 51.0 /MeV
// #define r0 1.25 /fm
// #define a 0.535 /fm
// #define beta2 0.0
// #define beta4 0.0

class Nucleon {
public:
    int charge;
    double mass;
    double position[3];   
    Nucleon(int ch, double m, double x, double y, double z) : charge(ch), mass(m) {
        position[0] = x;
        position[1] = y;
        position[2] = z;
    }


};

class Nucleus {
public:
    vector<Nucleon> nucleons;
    int charge;
    double mass;
    double R;
    double r0, V0, a, beta2, beta4;
    Nucleus(int ch, double m, double r0, double V0, double a, double b2, double b4)
        : charge(ch), mass(m), r0(r0), V0(V0), a(a), beta2(b2), beta4(b4) {
        R = r0 * pow(mass, 1.0/3.0);
        // Генерация нуклонов 
        int max_attempts = 100000; 
        for(int i = 0; i < max_attempts; i++){
            if (nucleons.size() >= mass) break;
            
        }
        }
};

class Ferrum : public Nucleus {
public:
    Ferrum() : Nucleus(26, 56.0, 1.25, 51.0, 0.535, 0.0, 0.0) {
        
        for (int i = 0; i < 56; ++i) {
            int nucleon_charge = (i < 26) ? 1 : 0; // Протоны и нейтроны
            double nucleon_mass = 1.0; // Упрощенное значение массы
            double x = static_cast<double>(rand()) / RAND_MAX * R;
            double y = static_cast<double>(rand()) / RAND_MAX * R;
            double z = static_cast<double>(rand()) / RAND_MAX * R;
            nucleons.emplace_back(nucleon_charge, nucleon_mass, x, y, z);
        }
    }
};
int main() {
    // Nucleus Ferrum(26, 56.0, 4.8, 1.25, 51.0, 0.535, 0.0, 0.0);

    return 0;
}