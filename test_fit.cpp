#include "fit_double_gaussian.cpp"

void test_fit()
{
	auto f = new TFile("exp8_4S_mva.root");
	auto lambda0 = (TTree*)f->Get("lambda0");

	TH1F * hist = new TH1F("hist", "hist", 200, 1.1, 1.13);
	lambda0->Draw("M>>hist", "mva>=0.001");
	
	double_gaussian_result res = fit_double_gaussian(hist, true, "mva");
	cout << "Status of the fit: " << res.status << endl;
}

