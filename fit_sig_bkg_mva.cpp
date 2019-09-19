#include "fit_double_gaussian.cpp"

// Find the numbers of signals and backgrounds for all the mva cuts of interest
// Default number of bins: 200
// Arguments:
//      filename    - write "mva nsig nbkg" to this file
//      append      - append to file if true
void fit_sig_bkg_mva(TString filename, bool append = true)
{
    auto input = new TFile("exp8_4S_mva.root");
    auto lambda0 = (TTree*)input->Get("lambda0");
    auto hist = new TH1F("hist", "hist", 200, 1.1, 1.13);

    std::ofstream fout;

    if(append)
    {
        fout.open(filename.Data(), std::ios_base::app);
    }
    else
    {
        fout.open(filename.Data());
    }

    for(double mva = 0.; mva <= .01; mva += 0.01)
    {
        TString cut = TString::Format("mva >= %.04f", mva);
        lambda0->Draw("M>>hist", cut, "goff");
        TString filename = TString::Format("mva_%03d.png", int(1000 * mva));

        double_gaussian_result res = fit_double_gaussian(hist, true, cut, filename);
        while(res.status != 0) {
            cout << "Fit does not seem right. Do it again for mva = " << mva << endl;
            res = fit_double_gaussian(hist, true, cut, filename);
        }
        cout << mva << ", " << res.nsig <<  ", " << res.nbkg << endl;
        fout << mva << ", " << res.nsig <<  ", " << res.nbkg << endl;
    }
    fout.close();
}