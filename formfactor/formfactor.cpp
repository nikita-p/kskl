#include <iostream>
#include <TComplex.h>
#include <TF1.h>

using namespace std;

class MDVM{

    static const double ALPHA;
    static const double C;
    
    static const double m_phi;
    static const double m_rho;
    static const double m_omg;
    
    static const double m_k0;
    static const double m_pi;
    static const double m_pi0;
    
    static const double w0_phi;
    static const double w0_rho;
    static const double w0_omg;
    
    static double BETA(double s);
    //бета
    static double PV2(double s, double M, double Mn); 
    //PV2: фазовый объём (почти: он типа нормирован) распада на 2 одинаковые частицы, s: энергия**2, M: масса векторного мезона, W0: ширина распада,  Mn: масса частицы в распаде
    static double PV3(double s, double m1, double m2, double m3);
    //PV3: фазовый объём распада на 3 частицы
    
    static double WRhoX(double s, double W0, double MX){
        return W0 * PV2(s, MX, m_pi);    }
    static double WOmgX(double s, double W0, double MX){
        return W0 * PV3(s,m_pi,m_pi,m_pi0) / PV3(MX*MX,m_pi,m_pi,m_pi0); }
    static double WPhiX(double s, double W0, double MX){
        return W0 * PV2(s, MX, m_k0);  }
    
    
    static TComplex BW(double, double, double, double(*WX)(double, double, double));
    //функция Брейта-Вигнера
    static TComplex BW_Rho(double s){
        return BW(s, m_rho, w0_rho, WRhoX); }
    static TComplex BW_Omg(double s){
        return BW(s, m_omg, w0_omg, WOmgX); }
    static TComplex BW_Phi(double s){
        return BW(s, m_phi, w0_phi, WPhiX); }
        
    static TComplex BW_Rho1(double s){
        return BW(s, 1465, 400, WRhoX); }
    static TComplex BW_Omg1(double s){
        return BW(s, 1420, 200, WOmgX); }
    static TComplex BW_Rho2(double s){
        return BW(s, 1570, 144, WRhoX); }
    static TComplex BW_Omg2(double s){
        return BW(s, 1650, 315, WOmgX); }
    static TComplex BW_Phi1(double s){
        return BW(s, 1680, 150, WPhiX); }
    static TComplex BW_Rho3(double s){
        return BW(s, 1720, 250, WRhoX); }
    static TComplex BW_Rho4(double s){ //из pdg видно, что это какая-то хренотень
        return BW(s, 1880, 160, WRhoX); }
    static TComplex BW_Rho5(double s){ //из pdg видно, что это какая-то полнейшая хренотень
        return BW(s, 2155, 320, WRhoX); }
    static TComplex BW_Phi2(double s){
        return BW(s, 2188, 83, WPhiX);  }
    
public:
    
    static TComplex F0(double* x, double* par, bool mode);
    //формфактор, нулевое приближение
    //mode: 0 - short/long; 1 - charged;
    static TComplex F1(double* x, double* par, bool mode);
    //формфактор с учётом omega(1400)
    static double Cross_Section(double* x, double* par, bool mode);
    static double Cross_Section_Neutral(double* x, double* par){
        return Cross_Section(x, par, 0);    }
    static double Cross_Section_Charged(double* x, double* par){
        return Cross_Section(x, par, 1);    }
    
    static TF1* Cross_Section(bool mode);
    
};

/*___________________________________________________________*/

const double MDVM::C     = 0.389379292E12; //(MeV)^2 * nb
const double MDVM::ALPHA = 7.297E-3;

const double MDVM::m_rho = 775.26;
const double MDVM::m_omg = 782.65;
const double MDVM::m_phi = 1019.464;//1;

const double MDVM::m_k0  = 497.6;
const double MDVM::m_pi  = 139.57;
const double MDVM::m_pi0 = 135.;

const double MDVM::w0_phi = 4.247;//9;
const double MDVM::w0_rho = 149.1;
const double MDVM::w0_omg = 8.49;

/*__________________________________________________________*/

    double MDVM::BETA(double s){ //бета
        double E = TMath::Sqrt(s)/2.;
        double P = TMath::Sqrt(E*E - m_k0*m_k0);
        return P/E;
    }

    double MDVM::PV2(double s, double M, double Mn){
        if(s/4. < Mn*Mn )
            return 0;
        double w = (M*M/s)*pow((s/4. - Mn*Mn)/((M*M)/4. - Mn*Mn),3/2.);
        return w;
    }
    
    double MDVM::PV3(double s, double m1, double m2, double m3){
        double pv = (pow(TMath::Pi(),3)/2.)*( pow(m1*m2*m3,1/2.)*pow(sqrt(s) - m1 - m2 - m3,2)/pow(m1 + m2 + m3,3/2.) );
        return pv;
    }

    TComplex MDVM::BW(double s, double MX, double WX0, double(*WX)(double, double, double)){
        TComplex I(0, 1);
        TComplex bw = pow( MX,2 )/( pow(MX,2) - s - I*sqrt(s)*WX(s,WX0,MX) );
        return bw;
    }
   
    TComplex MDVM::F0(double* x, double* par, bool mode){
        double n = 1.005;//1.027;
        double s = TMath::Power(x[0]*1E3, 2);
        double CR = par[0];
        double CO = par[1];
        double CP = par[2];
        
        double KR = mode ? CR/2. : -CR/2.;
        double KO = CO/6.;
        double KP = mode ? CP/3. : n*CP/3.;
        
        TComplex F = KR*BW_Rho(s) + KO*BW_Omg(s) + KP*BW_Phi(s);
        return F;
    }
    
    TComplex MDVM::F1(double* x, double* par, bool mode){
        double s = TMath::Power(x[0]*1E3, 2);
        double n = 1.027;
        double CR = par[3];
        double CR2 = par[5];
        double CR3 = par[6];
        //double CR4 = par[7];
        double CR5 = 1 - par[0] - par[3] - par[5] - par[6];// - par[7];
        double CO = par[4];
        double CO2 = par[7];//1 - par[1] - par[4];
        double CO3 = 1 - par[1] - par[4] - par[7];
        double CP = par[8];//1 - par[2];
        double CP2 = 1 - par[2] - par[8];
        
        double KR = mode ? CR/2. : -CR/2.;
        double KR2 = mode ? CR2/2. : -CR2/2.;
        double KR3 = mode ? CR3/2. : -CR3/2.;
        //double KR4 = mode ? CR4/2. : -CR4/2.;
        double KR5 = mode ? CR5/2. : -CR5/2.;
        double KO = CO/6.;
        double KO2 = CO2/6.;
        double KP = mode ? CP/3. : n*CP/3.;
        double KP2 = mode ? CP2/3. : n*CP2/3.;
        double KO3 = CO3/6.;
        
        TComplex F1 = F0(x, par, mode) + KR*BW_Rho1(s) + KO*BW_Omg1(s) + KP*BW_Phi1(s) + KR2*BW_Rho2(s) + KR3*BW_Rho3(s) + /*KR4*BW_Rho4(s) +*/ (KR5+KO3)*BW_Rho5(s) + KO2*BW_Omg2(s) + KP2*BW_Phi2(s);
        return F1;
    }

    double MDVM::Cross_Section(double* x, double* par, bool mode){
        double s = TMath::Power(x[0]*1E3, 2);
        if(x[0]<0.4976*2)
            return 0;
        double fabs = TComplex::Abs( F1(x, par, mode) );
        double cs = (TMath::Pi()/3.)*TMath::Power(ALPHA,2)*C*TMath::Power(BETA(s),3)*TMath::Power(fabs,2)/s;
        return cs;
    }
   
    TF1* MDVM::Cross_Section(bool mode){
        const int Npars = 9;
        TF1* fcs_c = new TF1("Cross section", (mode ? Cross_Section_Charged : Cross_Section_Neutral), 0.98, 2.1, Npars);
        fcs_c->SetParNames("C_{#rho}", "C_{#omega}", "C_{#phi}", "C_{#rho(1450)}",  "C_{#omega(1420)}", "C_{#rho(1570)}", "C_{#rho(1720)}", "C_{#omega(1650)}", "C_{#phi(1680)}");
        fcs_c->SetParameters(1.123, 1.027, 1.101, -0.1, -0.1, 0.1, -0.1, 0.1, -0.1);
        return fcs_c;
    }
    
    
    
