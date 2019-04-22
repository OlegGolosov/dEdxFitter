enum eParticles
{
  kPimin = 0,
  kKmin,
  kEmin,
  kPiplus,
  kKplus,
  kPplus,
  kEplus,
  kNparticles
};

TString histFilePath = "input/dedx_na49.root";
//TString histFilePath = "input/dedx_na61.root";
TString histName = "hTrackdEdx_allTPC";
TString cutsFilePath = "input/dedx_cuts.root";
double pScale = 20.;
vector <int> cuts = {kPimin, kPiplus, kPplus};
vector <int> particles = {kPimin, kKmin, kEmin, kPiplus, kKplus, kPplus, kEplus}; 

int nIterations = 1;

struct particleInfo 
{
  public:
    TString name;
    int charge;
    float mass;
};

unordered_map <int, particleInfo> particleList = 
{
  {kPimin,  {"pimin", -1, 0.140}},
  {kKmin,   {"kmin",  -1, 0.495}},
  {kEmin,   {"emin",  -1, 0.000511}},
  {kPiplus, {"piplus", 1, 0.140}},
  {kKplus,  {"kplus",  1, 0.495}},
  {kPplus,  {"pplus",  1, 0.938}},
  {kEplus,  {"eplus",  1, 0.000511}}
};

TH2* cutTH2 (TH2 *hFull, TCutG *cut);
Double_t bb_impl (Double_t *x, Double_t *par); // from analysis note
Double_t bb_impl_ (Double_t *x, Double_t *par); // from na49root

void fit_bb ()
{   
  int nParameters = 8;
  unordered_map <int, double*> parameters;
  unordered_map <int, TH2D*> hdEdxCut;
  unordered_map <int, TF1*> fitFunc;
  unordered_map <int, TString> funcName;
  unordered_map <int, float> fitXmin;
  unordered_map <int, float> fitXmax;
  
  int nParticles = particles.size ();
  int nCuts = cuts.size ();

  TCutG* graphCut;
  float xmin, xmax;
  TString cutName;
  int charge;
  float mass;
  int particleIndex;
  
  TFile *histFile = new TFile (histFilePath, "read");
  TFile *cutsFile = new TFile (cutsFilePath, "read");
  TH2D *hdEdx = (TH2D*) histFile -> Get (histName);
  double gXmin = hdEdx -> GetXaxis () -> GetXmin ();
  double gXmax = hdEdx -> GetXaxis () -> GetXmax ();
  
  for (auto particle : particles) 
  {
    parameters.emplace (particle, new double [nParameters]);
    funcName.emplace (particle, Form ("bb_%s", particleList.at(particle).name.Data()));
    parameters.at (particle) [0] = pScale;
    parameters.at (particle) [1] = particleList.at(particle).mass;
    { // na61?
//      parameters.at (particle) [2] = 0.42;
//      parameters.at (particle) [3] = 13.85;
//      parameters.at (particle) [4] = 2.28;
//      parameters.at (particle) [5] = 0.26;
    }
    { // na49
      parameters.at (particle) [2] = 1.613702;
      parameters.at (particle) [3] = 10.407406;
      parameters.at (particle) [4] = 2.463701;
      parameters.at (particle) [5] = 0.164279;
    }
    parameters.at (particle) [6] = 1.;
    parameters.at (particle) [7] = 1.;
    if (particleList.at(particle).charge > 0) {xmin = 0.; xmax = gXmax;}
    else {xmin = gXmin; xmax = 0.;}
    fitFunc.emplace (particle, new TF1 (funcName.at(particle), bb_impl_, xmin, xmax, 8));
  }
  
  TCanvas *c = new TCanvas ("fit");
  c -> Divide (2, 1);
  c -> cd (1);
  gPad -> SetLogz();
  
  for (auto particle : cuts)
  {
    cutName = particleList.at(particle).name.Data();
    graphCut = (TCutG*) cutsFile -> Get (cutName);
    hdEdxCut.emplace (particle, (TH2D*) cutTH2 (hdEdx, graphCut));
    fitXmin.emplace (particle, hdEdxCut.at (particle) -> GetXaxis () -> GetBinCenter (hdEdxCut.at (particle) -> FindFirstBinAbove (0, 1)));
    fitXmax.emplace (particle, hdEdxCut.at (particle) -> GetXaxis () -> GetBinCenter (hdEdxCut.at (particle) -> FindLastBinAbove (0, 1)));
    if (particleList.at(particle).charge > 0) {xmin = 0.; xmax = gXmax;}
    else {xmin = gXmin; xmax = 0.;}
    cout << xmin << "\t" << xmax << endl;
    
    fitFunc.at (particle) -> SetParameters (parameters.at (particle));
    fitFunc.at (particle) -> SetParNames ("pScale", "mass", "a1", "a2", "b", "c", "AB", "delta");
    fitFunc.at (particle) -> FixParameter (0, parameters.at (particle) [0]); // pScale
    fitFunc.at (particle) -> FixParameter (1, parameters.at (particle) [1]); // mass
    fitFunc.at (particle) -> FixParameter (2, parameters.at (particle) [2]); // a1
    fitFunc.at (particle) -> FixParameter (3, parameters.at (particle) [3]); // a2
    fitFunc.at (particle) -> FixParameter (4, parameters.at (particle) [4]); // b
    fitFunc.at (particle) -> FixParameter (5, parameters.at (particle) [5]); // c
    fitFunc.at (particle) -> FixParameter (6, 999.); // do not fit scaleFactor
    fitFunc.at (particle) -> FixParameter (7, 999.); // do not fit delta
    
//    for (int i = 0; i < nIterations; i++)
//    {
      hdEdxCut.at (particle) -> Fit (funcName.at(particle), "", "colz same", fitXmin.at (particle), fitXmax.at (particle));
//    }
    // iterate1
/*    fitFunc [particle] -> FixParameter (5, 10.);
    fitFunc [particle] -> FixParameter (6, 2.);
    fitFunc [particle] -> FixParameter (7, 0.);
    hdEdxCut [particle] -> Fit (funcName, "0");
    
    fitFunc [particle] -> FixParameter (4, fitFunc [particle] -> GetParameter (4));
    fitFunc [particle] -> ReleaseParameter (5)е;
    hdEdxCut [particle] -> Fit (funcName, "0");
    
    fitFunc [particle] -> FixParameter (5, fitFunc [particle] -> GetParameter (5));
    fitFunc [particle] -> ReleaseParameter (6);
    hdEdxCut [particle] -> Fit (funcName, "0");
    
    fitFunc [particle] -> FixParameter (6, fitFunc [particle] -> GetParameter (6));
    fitFunc [particle] -> ReleaseParameter (7);
    hdEdxCut [particle] -> Fit (funcName, "0");
*/
  }
  
  double *parameters_avg = new double [nParameters];
  
  parameters_avg [0] = pScale;
  for (int par = 2; par < nParameters; par++)
    parameters_avg [par] = 0.;
    
  for (auto particle : cuts)
  {
    for (int par = 2; par < nParameters; par++)
      parameters_avg [par] += (fitFunc.at (particle) -> GetParameter (par) / nCuts);
  }
  
  c -> cd (2);
  gPad -> SetLogz ();
  hdEdx -> Draw ("colz");
  for (auto particle : particles)
  { 
    fitFunc.at (particle) -> SetParameters (parameters_avg);
//    fitFunc.at (particle) -> SetParameters (parameters.at(particle));
    fitFunc.at (particle) -> SetParameter (1, particleList.at(particle).mass);
    fitFunc.at (particle) -> Draw ("same"); 
  }
}


Double_t bb_impl_ (Double_t *x, Double_t *par)
{
//  Returns expected dE/dx values, according to Bethe-Bloch parametrization
//  (Called by GetRelRise())
//  ** Code provided by Christof Roland **
 
  double Ptot = exp (fabs (x [0])) / par [0];
  double mass = par [1];

  // Parameters for Bethe Bloch-parametrization     
  double fParaC =  1.597546;// 1.5624;
  double fParaD =  9.8;
  double fParaE =  2.38;
  double fParaF =  0.2; 
  
  Float_t bg, dedx;
  Float_t c, d, e, f, x1, x2, p0, dfq, d1, d2, d3;
  
  bg = Ptot / mass;
  
  c =  fParaC;
  d =  fParaD;
  e =  fParaE;
  f =  fParaF;
  
  { // besides this part everything is taken from na49root
    c = par [2];
    d = par [3];
    e = par [4];
    f = par [5];
  }
  
  x1 = pow(10,(e - 1.0/3.0 * sqrt(2.0*log(10.0)/(3.0 * f))));
  x2 = pow(10,(e + 2.0/3.0 * sqrt(2.0*log(10.0)/(3.0 * f))));
  
  p0 = c/(d + 2.0*log(10.0)*e - 1.0);
  
  if (bg<x1)
    dfq = 0;
  else
  {
    dfq = -2.0 * log(bg) + 2.0 * log(10.0) * e;
    if(bg<x2)
    {
      d1 = 2.0/3.0*sqrt(2.0*log(10.0)/(3.0 * f));
      d2 = log(bg)/log(10.0);
      d3 = pow(( e + d1 - d2),3);
      dfq -= f*d3;
    }
  }
  
  dedx = p0*( (1+ bg*bg)/(bg*bg) * ( d + log(1+(bg*bg)) + dfq)-1.0 );
  
  return dedx;
}

Double_t bb_impl (Double_t *x, Double_t *par)
{
  double p = exp (fabs (x [0])) / par [0];
  double m = par [1];
  double E = sqrt (p * p + m * m);
  double v = p / m;
  double beta = p / E;
  double beta2 = (p * p) / (E * E);
  double gamma = sqrt (1. - beta2);
  double a1 = par [2];
  double a2 = par [3];
  double b = par [4];
  double c = par [5];
  double scaleFactor = par [6];
//  double scaleFactor = a1 / (a2 + 2.0 * log (10.0) * b - 1.0);
  double delta;
  double d = 3.;
  
  if (fabs (par [7] - 999.) > 0.1) // fit delta
    delta = par [7];
  else if (beta * gamma < a1) 
    delta = 0.;
  else if (beta * gamma < a2) 
    delta = 2. * (log (beta * gamma) - b) + c * pow (log (a2) - log (beta * gamma), d);
  else delta = 2. * (log (beta * gamma) - b);

//  double bb = par [1] / beta2 * (log (par [2] * v * beta2 / (1. - beta2)) - beta2); 
  double bb = scaleFactor / beta2 * (log (beta2 / (1. - beta2)) - beta2 - delta);
  return bb;
}

TH2* cutTH2 (TH2 *hFull, TCutG *cut)
{
  TString name = hFull -> GetName ();
  name += "_cut";
  TH2 *hCut = (TH2*) hFull -> Clone (name);
  int nBinsX = hCut -> GetNbinsX ();
  int nBinsY = hCut -> GetNbinsY ();
  double x, y;
  
  for (int i = 1; i <= nBinsX; i++) 
  {
    for (int j = 1; j <= nBinsY; j++) 
    {
      x = hCut -> GetXaxis () -> GetBinCenter (i);
      y = hCut -> GetYaxis () -> GetBinCenter (j);
      if (!cut -> IsInside (x, y)) 
      {
        hCut -> SetBinContent (i, j, 0.);
        hCut -> SetBinError (i, j, 0.);
      }
    }
//    hCut -> Draw ("colz");
//    gPad -> SetLogz ();
//    cut -> SetLineColor (kRed);
//    cut -> Draw ("same");
  }
  return hCut;
}
