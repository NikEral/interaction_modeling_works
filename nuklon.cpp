#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TF1.h"
#include "TBox.h"
#include "TPaveText.h"
class Nucleon {
private:
    	int charge; 
    	double mass;  
	double x, y, z;

public:
	Nucleon(int ch = 0, double m = 0.0, double x_pos = 0.0, double y_pos = 0.0, double z_pos = 0.0) : charge(ch), mass(m), x(x_pos), y(y_pos), z(z_pos) {}
    	int getCharge() const { return charge; }
    	double getMass() const { return mass; }
    	double getX() const { return x; }
    	double getY() const { return y; }
    	double getZ() const { return z; }
	double getRadius() const {
		return std::sqrt(x*x + y*y + z*z);
	}
	void setPosition(double x_pos, double y_pos, double z_pos) {
        	x = x_pos;
        	y = y_pos;
        	z = z_pos;
    	}
};
class Nucleus {
private:
	std::vector<Nucleon> nucleons;
    	int total_charge;
    	double total_mass;
    	double radius;

    	double r0;
    	double V0;
    	double a;
    	double beta2;
    	double beta4;

	double ro0;
public:
	Nucleus(int A = 0, double r0_val = 0.0, double V0_val = 0.0, double a_val = 0.0, double beta2_val = 0.0, double beta4_val = 0.0) : total_charge(0), total_mass(0), r0(r0_val), V0(V0_val), a(a_val), beta2(beta2_val), beta4(beta4_val), ro0(1.0) {
		radius = r0 * std::pow(A, 1.0/3.0);	
}
	

	const std::vector<Nucleon>& getNucleons() const { return nucleons; }
    	int getTotalCharge() const { return total_charge; }
    	double getTotalMass() const { return total_mass; }
    	double getRadius() const { return radius; }
    	int getTotalNucleons() const { return nucleons.size(); }
    	double getR0() const { return r0; }
    	double getV0() const { return V0; }
    	double getA() const { return a; }


	void generation(int protons, int neutrons) {
		nucleons.clear();
		total_charge = protons;
        	total_mass = protons * 938.27 + neutrons * 939.57;

		std::random_device rd;
        	std::mt19937 gen(rd());
        	std::uniform_real_distribution<> dis_phi(0, 2 * M_PI); //равномерное распределение для азимутального
        	std::uniform_real_distribution<> dis_cos(-1, 1); //равномерное распределение для косинуса зенитного

		for (int i = 0; i < protons; ++i) {
            		double r = VoodsSaxon(gen);
            		double phi = dis_phi(gen);
            		double cos_theta = dis_cos(gen);
            		double sin_theta = std::sqrt(1 - cos_theta * cos_theta);

            		double x = r * sin_theta * std::cos(phi);
            		double y = r * sin_theta * std::sin(phi);
            		double z = r * cos_theta;

            		nucleons.emplace_back(1, 938.27, x, y, z);
        	}
		for (int i = 0; i < neutrons; ++i) {
            		double r = VoodsSaxon(gen);
            		double phi = dis_phi(gen);
            		double cos_theta = dis_cos(gen);
            		double sin_theta = std::sqrt(1 - cos_theta * cos_theta);

            		double x = r * sin_theta * std::cos(phi);
            		double y = r * sin_theta * std::sin(phi);
            		double z = r * cos_theta;

            		nucleons.emplace_back(0, 939.57, x, y, z);
        	}
		
	}

	double VoodsSaxon(std::mt19937& gen) {
		std::uniform_real_distribution<> dis(0.0, 1.0);
        	double max_r = radius + 5 * a; // граница распределения
        
       		 while (true) {
            		double r = dis(gen) * max_r;
			double R = radius * (1 + beta2 * std::pow(r/radius, 2) + beta4 * std::pow(r/radius, 4));
            		double prob = ro0 / (1.0 + std::exp((r - R) / a));
            
            		if (dis(gen) < prob) {
                		return r;
            		}
        	}
	}

	void shift(double dx, double dy, double dz) {
        	for (auto& nucleon : nucleons) {
            		nucleon.setPosition(nucleon.getX() + dx, nucleon.getY() + dy, nucleon.getZ() + dz);
        	}
    	}
};



int count_overlap_nucleons(const Nucleus& nucleus1, const Nucleus& nucleus2, double distance, double& overlap1, double& overlap2) {
                overlap1 = 0;
                overlap2 = 0;
                int total_overlap = 0;

                double R1 = nucleus1.getRadius();
                double R2 = nucleus2.getRadius();

                for (const auto& nucleon : nucleus1.getNucleons()) {
                        double x = nucleon.getX();
                        double dist_to_center2 = std::sqrt(std::pow(x - distance, 2) + std::pow(nucleon.getY(), 2) + std::pow(nucleon.getZ(), 2));
                        if (dist_to_center2 <= R2) {
                                overlap1++;
                                total_overlap++;
                        }
                }

                for (const auto& nucleon : nucleus2.getNucleons()) {
                        double x = nucleon.getX();
                        double dist_to_center1 = std::sqrt(std::pow(x - distance, 2) + std::pow(nucleon.getY(), 2) + std::pow(nucleon.getZ(), 2));
                        if (dist_to_center1 <= R1) {
                                overlap2++;
                                total_overlap++;
                        }
                }

                return total_overlap;
}


double doubleGaussian(double* x, double* par) {
                return par[0] * std::exp(-0.5 * std::pow((x[0] - par[1])/par[2], 2)) + par[3] * std::exp(-0.5 * std::pow((x[0] - par[4])/par[5], 2));
}



void nuklon() {
	gROOT->SetBatch(kTRUE);
	int A = 32;
    	int protons = 16;
    	int neutrons = 16;
	double r0 = 1.21;
    	double V0 = 53.5;
    	double a = 0.52;   
    	double beta2 = 0.0;
    	double beta4 = 0.0;
	TH1D* hist_total = new TH1D("histAll", "Total overlap nucleons;Number of nucleons;Counts", 100, 0, protons + neutrons);
    	TH1D* hist1 = new TH1D("hist1", "Nucleus 1;Number of nucleons;Counts", 100, 0, protons + neutrons);
    	TH1D* hist2 = new TH1D("hist2", "Nucleus 2;Number of nucleons;Counts", 100, 0, protons + neutrons);
	int n = 10000;
	std::random_device rd;
    	std::mt19937 gen(rd());

	
	for (int i = 0; i < n; ++i) {

        	Nucleus nucleus1(A, r0, V0, a, beta2, beta4);
        	Nucleus nucleus2(A, r0, V0, a, beta2, beta4);

        	nucleus1.generation(protons, neutrons);
        	nucleus2.generation(protons, neutrons);
		
		std::uniform_real_distribution<> dis_distance(0.0, 2 * nucleus1.getRadius());
        	double distance = dis_distance(gen);

        	nucleus2.shift(distance, 0, 0);

        	double overlap1, overlap2;
        	int total_overlap = count_overlap_nucleons(nucleus1, nucleus2, distance, overlap1, overlap2);

        	hist_total->Fill(total_overlap);
        	hist1->Fill(overlap1);
        	hist2->Fill(overlap2);
    	}
	
	TF1* fit = new TF1("fit", doubleGaussian, 0, protons + neutrons, 6);
	double mean = hist_total->GetMean();
    	double rms = hist_total->GetRMS();

    	fit->SetParameters(hist_total->GetMaximum(), mean, rms/2, hist_total->GetMaximum()/2, mean + rms, rms/2);

    	fit->SetParNames("A1", "mean1", "sigma1", "A2", "mean2", "sigma2");

    	hist_total->Fit("fit", "Q");

	TCanvas* canvas = new TCanvas("canvas", "Voods-Saxon", 1920, 1080);
	hist_total->SetLineColor(kBlue);
    	hist_total->SetLineWidth(2);
	hist_total->SetFillColor(kBlue);
    	hist_total->SetFillStyle(3004);

    	hist1->SetLineColor(kGreen);
    	hist1->SetLineWidth(2);

    	hist2->SetLineColor(kRed);
    	hist2->SetLineWidth(2);

	hist_total->Draw();
	hist1->Draw("same");
    	hist2->Draw("same");
	fit->Draw("same");

	TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    	legend->AddEntry(hist_total, "Summ", "f");
    	legend->AddEntry(hist1, "Core 1", "l");
    	legend->AddEntry(hist2, "Core 2", "l");
    	legend->AddEntry(fit, "Double Gauss", "l");
    	legend->Draw();

	TFile* output = new TFile("Result.root", "RECREATE");
    	hist_total->Write();
    	hist1->Write();
    	hist2->Write();
    	fit->Write();
    	canvas->Write();
    	output->Close();

	canvas->SaveAs("Result.png");
}



























