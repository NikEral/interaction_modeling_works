/*
1.  Задать TLorentzVector частицы с массой лямбда гиперона (1.115 ГэВ) и компонентами импульса, случайно взятыми из равномерного распределения 0.2-1 ГэВ/c
2.  Перейти в систему покоя лямбды и создать TLorentzVector-ы для двух дочерних частиц (протон и отрицательно заряженный пион) так, чтобы выполнялись законы 
    сохранения энергии-импульса. 
3.  Вернуться обратно в лабораторную систему и записать в соответствующие гистограммы для каждого компонента импульса каждой частицы значения этих дочерних 
    частиц(должно получиться 6 гистограмм с px,py,pz для протонов и пионов)
4.  Провести пункты 1-3 10 000 раз, нарисовать гистограммы.
*/

#include <root/TLorentzVector.h>
/*
/usr/include/root/TLorentzVector.h:27:10: fatal error: Math/Vector4D.h: No such file or directory
   27 | #include "Math/Vector4D.h"
      |          ^~~~~~~~~~~~~~~~~
compilation terminated.
*/

#include <root/TRandom.h>
#include <root/TH1D.h>
#include <root/TCanvas.h>
#define ENTRIES 10000
#define MASS_LAMBDA 1.115 // GeV/c^2
#define MASS_PROTON 0.938 // GeV/c^2
#define MASS_PION 0.139 // GeV/c^2

int main() {
    TH1F *h_px_proton = new TH1F("h_px_proton", "Px Proton;Px (GeV/c);Counts", 100, -1, 1);
    TH1F *h_py_proton = new TH1F("h_py_proton", "Py Proton;Py (GeV/c);Counts", 100, -1, 1);
    TH1F *h_pz_proton = new TH1F("h_pz_proton", "Pz Proton;Pz (GeV/c);Counts", 100, -1, 1);
    
    TH1F *h_px_pion = new TH1F("h_px_pion", "Px Pion;Px (GeV/c);Counts", 100, -1, 1);
    TH1F *h_py_pion = new TH1F("h_py_pion", "Py Pion;Py (GeV/c);Counts", 100, -1, 1);
    TH1F *h_pz_pion = new TH1F("h_pz_pion", "Pz Pion;Pz (GeV/c);Counts", 100, -1, 1);

    TLorentzVector* lambda = new TLorentzVector();
    for (int i = 0; i < ENTRIES; i++) {
        Double_t x = gRandom->Uniform(0.2, 1.0);
        Double_t y = gRandom->Uniform(0.2, 1.0);
        Double_t z = gRandom->Uniform(0.2, 1.0);
        lambda->SetXYZM(x, y, z, MASS_LAMBDA);
        TVector3 beta_lambda = lambda->BoostVector();
        Double_t p_star = TMath::Sqrt((TMath::Power(MASS_LAMBDA, 2) - TMath::Power(MASS_PROTON + MASS_PION, 2))  *
                                      (TMath::Power(MASS_LAMBDA, 2) - TMath::Power(MASS_PROTON - MASS_PION, 2))) /
                                      (2 * MASS_LAMBDA);
        Double_t theta = gRandom->Uniform(0, TMath::Pi());
        Double_t phi = gRandom->Uniform(0, 2 * TMath::Pi());
        
        Double_t px_star = p_star * TMath::Sin(theta) * TMath::Cos(phi);
        Double_t py_star = p_star * TMath::Sin(theta) * TMath::Sin(phi);
        Double_t pz_star = p_star * TMath::Cos(theta);

        TLorentzVector proton;
        proton.SetXYZM(px_star, py_star, pz_star, MASS_PROTON);
        TLorentzVector pion;
        pion.SetXYZM(-px_star, -py_star, -pz_star, MASS_PION);

        proton.Boost(beta_lambda);
        pion.Boost(beta_lambda);

        h_px_proton->Fill(proton.Px());
        h_py_proton->Fill(proton.Py());
        h_pz_proton->Fill(proton.Pz());
        
        h_px_pion->Fill(pion.Px());
        h_py_pion->Fill(pion.Py());
        h_pz_pion->Fill(pion.Pz());
    }
    TCanvas *canvas = new TCanvas("canvas", "Lambda Decay", 1200, 800);
    canvas->Divide(3, 2);
    
    canvas->cd(1);
    h_px_proton->Draw();
    
    canvas->cd(2);
    h_py_proton->Draw();
    
    canvas->cd(3);
    h_pz_proton->Draw();
    
    canvas->cd(4);
    h_px_pion->Draw();
    
    canvas->cd(5);
    h_py_pion->Draw();
    
    canvas->cd(6);
    h_pz_pion->Draw();
    canvas->Print("hw3/hw3.png");
}