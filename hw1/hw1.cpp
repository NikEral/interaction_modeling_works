/*
1. Создать гистограмму со статистикой 100 000 вхождений, распределенной по функции: sin([0]x)+cos([1]*x)* x[2]
2. Нарисовать и вывести на рисунок значения среднего, среднеквадратичного отклонения, максимального и минимального значения, и интеграла в пределах (4;8)
3. Варианты значений параметров [0]-[2] (распределите поровну между собой)
    1. 1,4,2
    2. 3,7,1
    3. 2,3,-1
    4. 6,7,-2 */

#include <iostream>
#include <root/TH1F.h>
#include <root/TF1.h>
#include <root/TCanvas.h>
#include <root/TPaveText.h>
#define ENTRIES 1000000
#define LEFT_LIMIT 0
#define RIGHT_LIMIT 10
#define PARAM_0 1
#define PARAM_1 4
#define PARAM_2 2

int main() {
    TH1F* hist = new TH1F("hist", "Histogram;X;Entries", 100, LEFT_LIMIT, RIGHT_LIMIT);
    TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
    TF1* func = new TF1("func", "sin([0]*x) + cos([1]*x)*x*[2]", LEFT_LIMIT, RIGHT_LIMIT);
    func->SetParameters(PARAM_0, PARAM_1, PARAM_2);
    hist->FillRandom("func", ENTRIES);
    hist->Draw("hist");
    double mean = hist->GetMean();
    double std_dev = hist->GetStdDev();
    double max = hist->GetMaximum();
    double min = hist->GetMinimum();
    
    int bin4 = hist->FindBin(4.0);
    int bin8 = hist->FindBin(8.0);
    double integral = hist->Integral(bin4, bin8);
    
    TPaveText *stats = new TPaveText(0.7, 0.7, 1, 1, "NDC");
    stats->SetFillColor(0);
    stats->SetTextSize(0.03);
    stats->AddText(Form("Mean: %.4f", mean));
    stats->AddText(Form("StdDev: %.4f", std_dev));
    stats->AddText(Form("Max: %.4f", max));
    stats->AddText(Form("Min: %.4f", min));
    stats->AddText(Form("Integral (4-8): %.4f", integral));
    stats->Draw();
    canvas->Update();
    canvas->Print("hw1/hw1_hist.png");
}