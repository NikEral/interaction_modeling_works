void script(Int_t nbins, Double_t rlimit, Double_t llimit, Int_t entries, const char* filename){
    TFile* file = new TFile(filename, "recreate");
    TH1F* hist = new TH1F("hist", "Y;X;Counts", nbins, llimit, rlimit);
    hist->FillRandom("gaus");
    file->Write();
    file->Close();
}
