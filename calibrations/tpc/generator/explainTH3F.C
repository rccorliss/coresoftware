#include <TFile.h>
#include <TH3.h>
#include <TKey.h>
#include <TClass.h>
#include <TAxis.h>
#include <iostream>
#include <string>

/**
 * \brief Helper macro to inspect the binning of the first TH3 found in a ROOT file.
 * \param filename Path to the .root file to inspect.
 */
void explainTH3F(const std::string& filename)
{
  TFile* f = TFile::Open(filename.c_str());
  if (!f || f->IsZombie())
  {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  TH3* h = nullptr;
  TIter next(f->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*) next()))
  {
    TClass* cl = TClass::GetClass(key->GetClassName());
    if (cl && cl->InheritsFrom(TH3::Class()))
    {
      h = (TH3*) key->ReadObj();
      if (h) break;
    }
  }

  if (!h)
  {
    std::cerr << "Error: No TH3 object found in " << filename << std::endl;
    f->Close();
    return;
  }

  std::cout << "Inspecting TH3: " << h->GetName() << " (" << h->GetTitle() << ")" << std::endl;

  TAxis* axes[3] = {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()};
  std::string dims[3] = {"X", "Y", "Z"};

  for (int i = 0; i < 3; ++i)
  {
    std::cout << "\n" << dims[i] << " Axis \"" << axes[i]->GetName() << " = \"" << axes[i]->GetTitle() << "\"; range = [" << axes[i]->GetXmin() << ", " << axes[i]->GetXmax() << "], bins = " << axes[i]->GetNbins() << " dx:";
    int nbins = axes[i]->GetNbins();
    for (int b = 1; b <= nbins; ++b)
    {
      std::cout <<  axes[i]->GetBinWidth(b) << ", ";
    }
    std::cout << std::endl;
    }
  }

  f->Close();
}
