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
#define ATTEMPTS 10000
// #define A 56
// #define Z 26
// #define N 30
// #define r0 1.25 /fm
// #define V0 51.0 /MeV
// #define a 0.535 /fm
// #define beta2 0.0
// #define beta4 0.0

#include <root/TCanvas.h>
#include <root/TF1.h>
#include <root/TH1F.h>
#include <root/TRandom.h>
#include <root/TMath.h>

class Nucleon {
public:
    Int_t charge;
    Double_t mass;
    Double_t position[3];   
    Nucleon(Int_t ch, Double_t m, Double_t x, Double_t y, Double_t z) : charge(ch), mass(m) {
        position[0] = x;
        position[1] = y;
        position[2] = z;
    }
};

class Nucleus {
public:
    vector<Nucleon> nucleons;
    Int_t charge;
    Double_t mass;
    Double_t R;
    Double_t r0, V0, a, beta2, beta4;
    Double_t position[3];
    Nucleus(int ch, Double_t m, Double_t r0, Double_t V0, Double_t a, Double_t b2, Double_t b4, Double_t x=0, Double_t y=0, Double_t z=0)
        : charge(ch), mass(m), r0(r0), V0(V0), a(a), beta2(b2), beta4(b4), position{x, y, z} {
        R = r0 * pow(mass, 1.0/3.0);
        // Генерация нуклонов 
        Int_t max_attempts = 100000;
        TF1* ws_density = new TF1("ws_density", "[0]/(1+exp((x-[1])/[2]))", 0, 2*R);

        Double_t rho0 = 1;
        ws_density->SetParameters(rho0, R, a);
        Double_t integral = ws_density->Integral(0, R);
        rho0 /= integral;
        ws_density->SetParameters(rho0, R, a);

        for(Int_t i = 0; i < mass; i++){
            if (nucleons.size() >= mass) break;
            Double_t R = ws_density->GetRandom();
            Double_t theta = gRandom->Uniform(0, TMath::Pi());
            Double_t phi = gRandom->Uniform(0, 2 * TMath::Pi());
            Double_t x = R * TMath::Sin(theta) * TMath::Cos(phi);
            Double_t y = R * TMath::Sin(theta) * TMath::Sin(phi);
            Double_t z = R * TMath::Cos(theta);
            Int_t nucleon_charge = (i < charge) ? 1 : 0; // Протоны и нейтроны
            nucleons.push_back(Nucleon(nucleon_charge, 1.0, x+position[0], y+position[1], z+position[2]));
        }
    }
};

class Ferrum : public Nucleus {
public:
    Ferrum(Double_t x = 0, Double_t y = 0, Double_t z = 0) : Nucleus(26, 56.0, 1.25, 51.0, 0.535, 0.0, 0.0, x, y, z) {}
};
int main() {
    TCanvas* canvas = new TCanvas("canvas", "Nucleus Overlap", 800, 600);
    TH1F* hist_overlap_n1 = new TH1F("hist_overlap_n1", "Overlap Nucleons Nucleus 1;Number of Overlap Nucleons;Counts", 30, 0, 56);
    TH1F* hist_overlap_n2 = new TH1F("hist_overlap_n2", "Overlap Nucleons Nucleus 2;Number of Overlap Nucleons;Counts", 30, 0, 56);
    TH1F* hist_total_overlap = new TH1F("hist_total_overlap", "Total Overlap Nucleons;Number of Overlap Nucleons;Counts", 60, 0, 112);
    const Double_t r0 = 1.25;
    const Double_t mass = 56.0;
    const Double_t R = r0 * pow(mass, 1.0/3.0);
    for(Int_t i = 0; i < ATTEMPTS; i++){
        Double_t distance = gRandom->Uniform(0, 2*R);
        Ferrum nucleus1;
        Ferrum nucleus2(distance, 0, 0);
        Int_t count_overlap_n1 = 0;
        for(auto n : nucleus1.nucleons){
        /*
        n1 in nucleus2 -> (x - x2)^2 + (y - y2)^2 + (z - z2)^2 <= R^2
        */
            if(pow(n.position[0] - nucleus2.position[0], 2) +
               pow(n.position[1] - nucleus2.position[1], 2) +
               pow(n.position[2] - nucleus2.position[2], 2) <= pow(R, 2)){
                count_overlap_n1++;
            }
        }
        hist_overlap_n1->Fill(count_overlap_n1);
        Int_t count_overlap_n2 = 0;
        for(auto n : nucleus2.nucleons){
            if(pow(n.position[0] - nucleus1.position[0], 2) +
               pow(n.position[1] - nucleus1.position[1], 2) +
               pow(n.position[2] - nucleus1.position[2], 2) <= pow(R, 2)){
                count_overlap_n2++;
            }
        }
        hist_overlap_n2->Fill(count_overlap_n2);
        hist_total_overlap->Fill(count_overlap_n1 + count_overlap_n2);
    }
    hist_total_overlap->SetLineColor(kGreen);
    hist_total_overlap->Draw("hist");
    hist_overlap_n1->SetLineColor(kRed);
    hist_overlap_n1->Draw("hist same");
    hist_overlap_n2->SetLineColor(kBlue);
    hist_overlap_n2->Draw("hist same");
    canvas->Print("overlap_nucleons.png");

    return 0;
}