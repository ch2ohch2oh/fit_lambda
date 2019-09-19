// Double Gaussian fitter
// Dazhi W. @ September 2019
//
// This script is a double Gaussian fitter. The signal distribution is 
// assumed to be a double Gussian (two Gaussians with identical means).
// The background is taken to be second order polynomial.
// The initial values of the 3 coeeficients of the polynomial are assigned
// uniform random numbers from -1 to +1 so that the user do not have to
// change the parameter every time the fit fails.
//
// The signal and background numbers are guessed based on sideband 
// subtraction.
// 

using namespace RooFit;

// An wrapper to return the result of the fit
class double_gaussian_result
{
public:
    double_gaussian_result(RooFitResult * res)
    {
        this->fit_result = res;
    }

    RooFitResult * fit_result;
    double nsig;
    double nbkg;
    int status;
};

// Estimate the number of signals based on sideband subtraction
// Arguments:
//      hist  - histogram containing the data
//      lower - lower bound of the signal peak
//      upper - upper bound of the signal peak
double estimate_sig(TH1F* hist, double lower, double upper)
{
    // ROOT calculates integral based on bin index. Silly!
    int bin_min = 0;
    int bin_max = hist->GetSize();
    int bin_lower = hist->FindBin(lower);
    int bin_upper = hist->FindBin(upper);

    //cout << x_min << x_max << endl;
    double sideband_density = (hist->Integral(bin_min, bin_lower) + hist->Integral(bin_upper, bin_max)) / 
        (bin_max - bin_upper + bin_lower - bin_min);
    return hist->Integral(bin_lower, bin_upper) - sideband_density * (bin_upper - bin_lower);
}

// Fit a double Gaussian + second order polynomial background
// Arguments:
//      hist        - data histogram     
//      make_plot   - make plot or not
//      title       - title on the plot 
//      filename    - filename for the output figure
double_gaussian_result fit_double_gaussian(TH1F* hist, bool make_plot = false, TString title = "Mass", TString filename = "")
{
    // This is the range of the fit
    double mass_min = hist->GetXaxis()->GetXmin();
    double mass_max = hist->GetXaxis()->GetXmax();
    
    // Peak position
    double mass_mean = 1.1156;
    // Constrain the mass peak winthin this mass_mean += mass_windows/2
    double mass_window = 0.01;
    // Mass window for estimation of signal number via sideband subtraction
    double signal_window = 0.02;

    double estimated_total = hist->GetEntries();
    double estimated_sig = estimate_sig(hist, mass_mean - signal_window / 2, mass_mean + signal_window / 2);
    double estimated_bkg = estimated_total - estimated_sig;
    // Make it safe
    estimated_total *= 1.5;


    // Narrow Gaussian
    double sigma1_min = 0.0001;
    double sigma1_max = 0.005;
    double sigma1_default = 0.001;
    
    // Wide Gaussian
    double s2overs1_min = 1;
    double s2overs1_max = 10;
    double s2overs1_default = 5;

    // Double Gaussian signal
    RooRealVar mass("mass", "Mass", mass_min, mass_max);
    RooRealVar sigma1("sigma1", "Sigma of the narrow Gaussian", sigma1_default, sigma1_min, sigma1_max);
    RooRealVar s2overs1("s2overs1", "Ratio of the two sigmas", s2overs1_default, s2overs1_min, s2overs1_max);
    RooFormulaVar sigma2("sigma2", "Sigma of the wide Gaussian" , "@0 * @1", RooArgList(sigma1, s2overs1));
    RooRealVar mean("mean", "Center of the mass", mass_mean, mass_mean - mass_window / 2, mass_mean + mass_window / 2);
    RooGaussian g1("g1", "Narrow Gaussian", mass, mean, sigma1);
    RooGaussian g2("g2", "Wide Gaussian", mass, mean, sigma2);
    
    // Signal yield
    RooRealVar nsig1("nsig1", "Yield for wide Gaussian", 0.5 * estimated_sig, 0, estimated_total);
    RooRealVar nsig2("nsig2", "Yield for the narrow Gaussian", 0.5 * estimated_sig, 0, estimated_total);

    // Second order polynomial background
    time_t timer;
    auto ra = new TRandom(time(&timer));
    RooRealVar c0("c0", "c0", ra->Uniform(-1, 1), -10, 10);
    RooRealVar c1("c1", "c1", ra->Uniform(-1 ,1), -10, 10);
    RooRealVar c2("c2", "c2", ra->Uniform(-1, 1), -10, 10);
    RooPolynomial bkg("bkg", "Polynomial background", mass, RooArgList(c0, c1, c2));
    
    // Background yield
    RooRealVar nbkg("nbkg", "Background yield", estimated_bkg, 0, estimated_total);

    RooAddPdf model("model", "Signal + Background", RooArgList(g1, g2, bkg), RooArgList(nsig1, nsig2, nbkg));

    RooDataHist data("data", "Data", mass, hist);

    RooFitResult * fit_res = model.fitTo(data, Extended(kTRUE), Save());

    double total_sig = nsig1.getVal() + nsig2.getVal();

    if(make_plot)
    {
        TCanvas * can = new TCanvas("can", "Mass");
        RooPlot * frame = mass.frame(Title(title));
        data.plotOn(frame, Binning(100));
        model.plotOn(frame);
        model.plotOn(frame, LineColor(kRed));
        model.plotOn(frame, Components("bkg"));
        model.paramOn(frame, Format("NEU",AutoPrecision(1)));
        frame->Draw();

        // Make model_paramBox transparent
        // Have to use GetPrimitive otherwise error
        auto model_paramBox = (TPaveText*)can->GetPrimitive("model_paramBox");
        model_paramBox->SetFillStyle(0);
        model_paramBox->SetBorderSize(0);
        model_paramBox->SetTextSizePixels(13);
        model_paramBox->AddText(TString::Format("total sig = %.0f", total_sig));
        model_paramBox->AddText(TString::Format("sig / bkg = %.2f", total_sig / nbkg.getVal()));

        if(filename != "")
        {
            can->SaveAs(filename);
        }
    }

    double_gaussian_result res(fit_res);
    res.nsig = total_sig;
    res.nbkg = nbkg.getVal();
    res.status = fit_res->status();

    return res;
}