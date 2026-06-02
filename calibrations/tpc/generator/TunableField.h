#ifndef TUNABLEFIELD_H
#define TUNABLEFIELD_H

#include <TVector3.h>
#include <string>
#include <vector>

class RosseggerReader;
class PadrowReader;
class TH3F;

/**
 * \class TunableField
 * \brief Implements a tunable electric field model for the sPHENIX TPC 
 * accounting for primary space charge and IBF.
 */
class TunableField
{
 public:
  TunableField();
  virtual ~TunableField();

  // Part A: Primary Space Charge Handling
  void setSpaceChargeFilename(const std::string& filename) { m_scFilename = filename; }
  void loadPrimarySpaceCharge();
  void makePrimaryInto2D();
  void computePrimaryField();

  // Part B: IBF Handling and Rossegger Lookup
  void loadRossegger();
  void loadIBFSpaceCharge();
  void makeIBFinto2D();
  void makeIBFRotationallyPeriodic();
  void loadPadrowBoundaries();
  void adjustIBFToActiveArea();
  void computePadrowIBFFields();

  void calculateFieldContributions();

  // Diagnostic Output
  void setDiagnosticFile(const std::string& filename);
  void saveDiagnosticSpaceCharge(const std::string& name, TH3F* h);

  // Part C: Total Field Calculation
  void setPrimaryScale(float scale) { m_primaryScale = scale; }
  void setIBFGain(float gain) { m_ibfGain = gain; }
  void setExternalField(float E) { m_externalE = E; }
  void setPadrowKnockout(int padrow, bool knockout);
  void generateTotalField();  

  TVector3 getFieldAt(const TVector3& pos);

 private:
  void clearPadrowFields();
  void rotateAndAverage(TH3F* h);
  void FillGuardBins(TH3F* h);
  void computeField(TH3F* sourceCharge, TH3F* fieldR, TH3F* fieldP, TH3F* fieldZ);

  TH3F* makeUnguardedStandardTH3F(const std::string& name, const std::string& title);
  TH3F* makeGuardedStandardTH3F(const std::string& name, const std::string& title);

  RosseggerReader* m_rosseggerReader = nullptr;
  
  TH3F* m_primaryCharge = nullptr;
  TH3F* m_ibfCharge = nullptr;
  std::vector<TH3F*> m_padrowCharges;

  TH3F* m_primaryFieldR = nullptr;
  TH3F* m_primaryFieldP = nullptr;
  TH3F* m_primaryFieldZ = nullptr;

  TH3F* m_totalIbfFieldR = nullptr;
  TH3F* m_totalIbfFieldP = nullptr;
  TH3F* m_totalIbfFieldZ = nullptr;

  std::vector<TH3F*> m_padrowFieldR;
  std::vector<TH3F*> m_padrowFieldP;
  std::vector<TH3F*> m_padrowFieldZ;

  TFile* m_diagFile = nullptr;

  float m_rmin = 20.0;
  float m_rmax = 78.0;
  float m_zmax = 105.5;
  int m_nr = 50;
  int m_nphi = 36;
  int m_nz = 50;

  std::string m_scFilename;
  std::vector<std::pair<float, float>> m_padrowBoundaries;

  float m_primaryScale = 1.0;
  float m_ibfGain = 1.0;
  float m_externalE = 400.0; // V/cm
  std::vector<bool> m_padrowKnockout;

  bool m_fieldIsReady = false;
  bool m_componentFieldsAreReady = false;
};
#endif // TUNABLEFIELD_H