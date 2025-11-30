/*
Создать дерево TTree с полями: 
Температура(распределено по гауссу с параметрами(1,0,1) +4)
Давление(распределено по гауссу с параметрами(1,0,5)+750)
Время(в формате числа минут, каждое следующее значение +1 минута)
Дождь(Значение 0 или 1, определяется случайно)
Число вхождений: 1440
Написать функцию, которая будет выводить все значения в момент времени формата часы:минуты.
*/
#include <iostream>
#include <root/TTree.h>
#include <root/TRandom.h>
#include <root/TFile.h>
#define ENTRIES 1440

using namespace std;

void format_print(TTree *tree, const char* time_str) {
    // Тут хорошо было бы написать парсер времени, но времени нет, поэтому сделаем костыль
    Int_t hours = atoi(time_str);
    Int_t minutes = atoi(time_str + 3);
    Int_t total_minutes = hours * 60 + minutes;
    Double_t temperature;
    Double_t pressure;
    Int_t time;
    Bool_t rain;

    tree->SetBranchAddress("temperature", &temperature);
    tree->SetBranchAddress("pressure", &pressure);
    tree->SetBranchAddress("time", &time);
    tree->SetBranchAddress("rain", &rain);      
    tree->GetEntry(total_minutes);
    cout << "Time: " << time_str << endl;
    cout << "Temperature: " << temperature << endl;
    cout << "Pressure: " << pressure << endl;
    cout << "Rain: " << (rain ? "Yes" : "No") << endl;
}

int main(){
    TFile *file = new TFile("hw2/hw2.root", "recreate");
    Double_t temperature;
    Double_t pressure;
    Int_t time;
    Bool_t rain;
    TTree* tree = new TTree("tree", "tree");
    tree->Branch ("temperature", &temperature, "temperature/D");
    tree->Branch("pressure", &pressure, "pressure/D");
    tree->Branch("time", &time, "time/I");
    tree->Branch("rain", &rain, "rain/O");

    for (int i = 0; i < ENTRIES; i++) {
        temperature = gRandom->Gaus(1, 1) + 4;
        pressure = gRandom->Gaus(1, 5) + 750;
        time = i;
        rain = gRandom->Integer(2);
        tree->Fill();
    }
    format_print(tree, "04:20");
    file->Write();
    file->Close();
    // tree->Print();
}