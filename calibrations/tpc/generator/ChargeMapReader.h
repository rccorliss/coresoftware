class TH3F;
class TTree;
class TVector;

class MultiArray;


//since we are never in the position of adding very large numbers to very small, floats are sufficient precision here.
class ChargeMapReader
{
 public:
  ChargeMapReader();//calls the below with default values.
  ChargeMapReader(int _n0, float _rmin, float _rmax, int _n1, float _phimin, float _phimax, int _n2,float _zmin, float _zmax);
  ~ChargeMapReader();

 private:
  MultiArray<float> *charge=nullptr;
  TH3* hSourceCharge=nullptr;
  TH3* hChargeDensity=nullptr;
  bool chargeHistExists=true;
  bool chargeArrayExists=false;
  int nBins[]={1,1,1};//r,phi,z bins of the output fixed-width array
  float lowerBound[]={0,0,0};
  float upperBound[]={999,999,999};
  float binWidth[]={999,999,999};

  bool CanInterpolateAt(float r, float phi, float z);//checks whether it is okay to interpolate at this position in the charge density hist
  void RegenerateCharge();//internal function to revise the internal array whenever the bounds change etc.
  void RegenerateDensity();//internal function to rebuild the charge density map when the input map changes.

 public:
  ChargeMapReader();
  ~ChargeMapReader();
  bool ReadSourceCharge(const char* filename, const char* histname);
  bool ReadSourceCharge(TH3* sourceHist);
  bool SetOutputParameters(int _nr, float _rmin, float _rmax, int _nphi, float _phimin, float _phimax, int _nz,float _zmin, float _zmax);
  bool SetOutputBounds(float _rmin, float _rmax, float _phimin, float _phimax, float _zmin, float _zmax);
  bool SetOutputBins(int _nr, int _nphi, int _nz);
  float GetChargeInBin(int r, int phi, int z);
  float GetChargeAtPosition(float r, float phi, float z);
}
