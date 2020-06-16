#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

//jet class definitons
class jets
{
double pt, eta, phi, m;
int btag;

public:
jets() = default;
jets(double, double, double, double, int);

void setAll(double, double, double, double, int);
void setPt(double);
void setEta(double);
void setPhi(double);
void setM(double);
void setBTag(int);


double getPt();
double getEta();
double getPhi();
double getM();
int    getBTag();
TLorentzVector get4P();
};

// constructor and method definitions
inline jets::jets(double Pt, double Eta, double Phi, double M, int BTag)
{
        pt   = Pt;
        eta  = Eta;
        phi  = Phi;
        m    = M;
        btag = BTag;
}

void jets::setAll(double Pt, double Eta, double Phi, double M, int BTag)
{
        pt   = Pt;
        eta  = Eta;
        phi  = Phi;
        m    = M;
        btag = BTag;
}

void jets::setPt(double Pt) {
        pt = Pt;
}
void jets::setEta(double Eta) {
        eta = Eta;
}
void jets::setPhi(double Phi) {
        phi = Phi;
}
void jets::setM(double M) {
        m = M;
}
void jets::setBTag(int BTag) {
        btag = BTag;
}

double jets::getPt() {
        return pt;
}
double jets::getEta() {
        return eta;
}
double jets::getPhi() {
        return phi;
}
double jets::getM() {
        return m;
}
int jets::getBTag() {
        return btag;
}

TLorentzVector jets::get4P()
{
        TLorentzVector dummy4P;
        dummy4P.SetPtEtaPhiM(pt, eta, phi, m);
        return dummy4P;
}

//lepton class definitons
class leptons
{
double pt, eta, phi, m, q;
int pdgid;

public:
leptons() = default;
leptons(double, double, double, double, double, int);

void setAll(double, double, double, double, double, int);
void setPt(double);
void setEta(double);
void setPhi(double);
void setM(double);
void setQ(double);
void setPDGID(int);

double getPt();
double getEta();
double getPhi();
double getM();
double getQ();
int    getPDGID();
TLorentzVector get4P();
};

// constructor and method definitions
inline leptons::leptons(double Pt, double Eta, double Phi, double M, double Q, int PDGID)
{
        pt = Pt;
        eta = Eta;
        phi = Phi;
        m = M;
        q = Q;
        pdgid = PDGID;
}

void leptons::setAll(double Pt, double Eta, double Phi, double M, double Q, int PDGID)
{
        pt = Pt;
        eta = Eta;
        phi = Phi;
        m = M;
        q = Q;
        pdgid = PDGID;
}

void leptons::setPt(double Pt) {
        pt = Pt;
}
void leptons::setEta(double Eta) {
        eta = Eta;
}
void leptons::setPhi(double Phi) {
        phi = Phi;
}
void leptons::setM(double M) {
        m = M;
}
void leptons::setQ(double Q) {
        q = Q;
}
void leptons::setPDGID(int PDGID) {
        pdgid = PDGID;
}

double leptons::getPt() {
        return pt;
}
double leptons::getEta() {
        return eta;
}
double leptons::getPhi() {
        return phi;
}
double leptons::getM() {
        return m;
}
double leptons::getQ() {
        return q;
}
int leptons::getPDGID() {
        return pdgid;
}

TLorentzVector leptons::get4P()
{
        TLorentzVector dummy4P;
        dummy4P.SetPtEtaPhiM(pt, eta, phi, m);
        return dummy4P;
}

//candidate class definitons
class candidate
{
double pt, eta, phi, m;

public:
candidate() = default;

candidate(double, double, double, double);
candidate(jets, jets);
candidate(leptons, leptons);
candidate(candidate, jets);

void setAll(double, double, double, double);
void setPt(double);
void setEta(double);
void setPhi(double);
void setM(double);


double getPt();
double getEta();
double getPhi();
double getM();

TLorentzVector get4P();
};

// constructor and method definitions
inline candidate::candidate(double Pt, double Eta, double Phi, double M)
{
        pt   = Pt;
        eta  = Eta;
        phi  = Phi;
        m    = M;
}

inline candidate::candidate(jets daughter1, jets daughter2)
{
        TLorentzVector candidate_4P;
        candidate_4P = daughter1.get4P() + daughter2.get4P();

        pt  = candidate_4P.Pt();
        eta = candidate_4P.Eta();
        phi = candidate_4P.Phi();
        m   = candidate_4P.M();
}

inline candidate::candidate(leptons daughter1, leptons daughter2)
{
        TLorentzVector candidate_4P;
        candidate_4P = daughter1.get4P() + daughter2.get4P();

        pt  = candidate_4P.Pt();
        eta = candidate_4P.Eta();
        phi = candidate_4P.Phi();
        m   = candidate_4P.M();
}

inline candidate::candidate(candidate daughter1, jets daughter2)
{
        TLorentzVector candidate_4P;
        candidate_4P = daughter1.get4P() + daughter2.get4P();

        pt  = candidate_4P.Pt();
        eta = candidate_4P.Eta();
        phi = candidate_4P.Phi();
        m   = candidate_4P.M();
}


void candidate::setAll(double Pt, double Eta, double Phi, double M)
{
        pt   = Pt;
        eta  = Eta;
        phi  = Phi;
        m    = M;
}

void candidate::setPt(double Pt) {
        pt = Pt;
}
void candidate::setEta(double Eta) {
        eta = Eta;
}
void candidate::setPhi(double Phi) {
        phi = Phi;
}
void candidate::setM(double M) {
        m = M;
}

double candidate::getPt() {
        return pt;
}
double candidate::getEta() {
        return eta;
}
double candidate::getPhi() {
        return phi;
}
double candidate::getM() {
        return m;
}

TLorentzVector candidate::get4P()
{
        TLorentzVector dummy4P;
        dummy4P.SetPtEtaPhiM(pt, eta, phi, m);
        return dummy4P;
}


//////////////////////////////function definitions///////////////////////////////////////////////

double get_event_shape_variable(std::vector<jets> v_jets, std::string option = "aplanarity")
{
        if (v_jets.size() > 0)
        {
                //construct the sphericity tensor
                double Sxx(0), Sxy(0), Sxz(0);
                double Syx(0), Syy(0), Syz(0);
                double Szx(0), Szy(0), Szz(0);

                double sump2(0), sumE(0);

                for (unsigned int i_jet = 0; i_jet < v_jets.size(); ++i_jet)
                {
                        Sxx += v_jets.at(i_jet).get4P().Px() * v_jets.at(i_jet).get4P().Px();
                        Sxy += v_jets.at(i_jet).get4P().Px() * v_jets.at(i_jet).get4P().Py();
                        Sxz += v_jets.at(i_jet).get4P().Px() * v_jets.at(i_jet).get4P().Pz();

                        Syy += v_jets.at(i_jet).get4P().Py() * v_jets.at(i_jet).get4P().Py();
                        Syz += v_jets.at(i_jet).get4P().Py() * v_jets.at(i_jet).get4P().Pz();

                        Szz += v_jets.at(i_jet).get4P().Pz() * v_jets.at(i_jet).get4P().Pz();

                        sump2 += pow(v_jets.at(i_jet).get4P().E(), 2) - pow(v_jets.at(i_jet).get4P().M(), 2);
                        sumE  += v_jets.at(i_jet).get4P().E();
                }

                if (sump2 > 0)
                {
                        Sxx = Sxx / sump2;
                        Sxy = Sxy / sump2;
                        Sxz = Sxz / sump2;

                        Syy = Syy / sump2;
                        Syz = Syz / sump2;

                        Szz = Szz / sump2;
                }

                TMatrixDSym S(3);
                S(0, 0) = Sxx;
                S(0, 1) = Sxy;
                S(0, 2) = Sxz;
                S(1, 0) = Sxy;
                S(1, 1) = Syy;
                S(1, 2) = Syz;
                S(2, 0) = Sxz;
                S(2, 1) = Syz;
                S(2, 2) = Szz;

                //---- compute sphericity -----------------------
                TMatrixDSymEigen TEigen(S);
                TVectorD eigenValues(TEigen.GetEigenValues());

                double aplanarity = 1.5 * eigenValues(2);
                double sphericity = 1.5 * (eigenValues(1) + eigenValues(2));

                if (option == "aplanarity")
                        return aplanarity;
                else if (option == "sphericity")
                        return sphericity;
        }
}

double get_sumE(std::vector<jets> v_jets)
{
        double sumE(0);

        for (auto& jet : v_jets)
                sumE += jet.get4P().E();

        return sumE;
}

double legendre_l(int l, double x)
{
        if (l == 0) return 1;
        else if (l == 1) return x;
        else if (l == 2) return 0.5 * (3 * x * x - 1);
        else return 0.5 * (5 * x * x * x - 3 * x);
}

double cos_omega_ij(jets j1, jets j2)
{
        //"calculates the 3D angle between j1 and j2"
        TVector3 p1 = j1.get4P().Vect();
        TVector3 p2 = j2.get4P().Vect();
        return p1.Dot(p2) / (p1.Mag() * p2.Mag());
}

std::vector<float> get_fox_wolfram_parameters(std::vector<jets> v_jets, int l_max)
{
        std::vector<float> fox_wolfram(l_max, 0.);

        //---- compute Fox-Wolfram Moment-----------------------
        double px1, py1, pz1, px2, py2, pz2, p_i, p_j;
        double E_vis2 = get_sumE(v_jets) * get_sumE(v_jets);
        double cos_omega;

        for (unsigned int i_jet = 0; i_jet < v_jets.size(); ++i_jet)
        {
                for (unsigned int j_jet = 0; j_jet < v_jets.size(); ++j_jet)
                {
                        if (i_jet <= j_jet)
                        {
                                cos_omega = cos_omega_ij(v_jets.at(i_jet), v_jets.at(j_jet));
                                p_i = v_jets.at(i_jet).get4P().P();
                                p_j = v_jets.at(j_jet).get4P().P();

                                for (int l = 0; l < l_max; ++l)
                                {
                                        fox_wolfram.at(l) += (p_i*p_j/E_vis2) * legendre_l(l, cos_omega);
                                }
                        }
                }
        }
        return fox_wolfram;
}

std::vector<candidate> create_W_candidates(std::vector<jets> v_light_jets)
{
        std::vector<candidate> v_W_candidates;
        for (unsigned int i_jet = 0; i_jet < v_light_jets.size(); ++i_jet)
        {
                for (unsigned int j_jet = 0; j_jet < v_light_jets.size(); ++j_jet)
                {
                        if (i_jet < j_jet)
                        {
                                candidate dummy_candidate(v_light_jets.at(i_jet), v_light_jets.at(j_jet));
                                v_W_candidates.push_back(dummy_candidate);
                        }
                }
        }

        return v_W_candidates;
}

std::vector<candidate> create_had_top_candidates(std::vector<jets> v_b_jets, std::vector<candidate> v_W_candidates, bool sort = 1)
{
        std::vector<candidate> v_had_top_candidates;

        for (unsigned int i_w_candidate = 0; i_w_candidate < v_W_candidates.size(); ++i_w_candidate)
        {
                for (unsigned int j_bjet = 0; j_bjet < v_b_jets.size(); ++j_bjet)
                {
                        candidate dummy_candidate(v_W_candidates.at(i_w_candidate), v_b_jets.at(j_bjet));
                        v_had_top_candidates.push_back(dummy_candidate);
                }
        }

        std::sort(v_had_top_candidates.begin(), v_had_top_candidates.end(), [](candidate first, candidate second) {
                return abs(first.getM() - 173.) < abs(second.getM() - 173.);
        });

        return v_had_top_candidates;
}

std::vector<jets> create_b_jet_collection(std::vector<jets> v_jets)
{
        std::vector<jets> v_b_jets;

        for (unsigned int i_jet = 0; i_jet < v_jets.size(); ++i_jet)
        {
                if (v_jets.at(i_jet).getBTag() == 1)
                {
                        v_b_jets.push_back(v_jets.at(i_jet));
                }
        }

        return v_b_jets;
}

std::vector<jets> create_light_jet_collection(std::vector<jets> v_jets)
{
        std::vector<jets> v_light_jets;

        for (unsigned int i_jet = 0; i_jet < v_jets.size(); ++i_jet)
        {
                if (v_jets.at(i_jet).getBTag() != 1)
                {
                        v_light_jets.push_back(v_jets.at(i_jet));
                }
        }

        return v_light_jets;
}

std::vector<float> create_jet_pt_vector(std::vector<jets> v_jets, int n)
{
        std::vector<float> jet_pt;

        for (int i_jet = 0; i_jet < n; ++i_jet)
        {
                jet_pt.push_back(v_jets.at(i_jet).getPt());
        }
        return jet_pt;
}
