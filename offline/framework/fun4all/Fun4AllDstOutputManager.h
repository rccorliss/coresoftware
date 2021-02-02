// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLDSTOUTPUTMANAGER_H
#define FUN4ALL_FUN4ALLDSTOUTPUTMANAGER_H

#include "Fun4AllOutputManager.h"

#include <set>
#include <string>

class PHNodeIOManager;
class PHCompositeNode;

class Fun4AllDstOutputManager : public Fun4AllOutputManager
{
 public:
  Fun4AllDstOutputManager(const std::string &myname = "DSTOUT", const std::string &filename = "dstout.root");
  virtual ~Fun4AllDstOutputManager();

  int AddNode(const std::string &nodename);
  int AddRunNode(const std::string &nodename);
  int StripNode(const std::string &nodename);
  int StripRunNode(const std::string &nodename);
  void SaveRunNode(const int i) {m_SaveRunNodeFlag = i;}
  int outfileopen(const std::string &fname);

  void Print(const std::string &what = "ALL") const;

  int Write(PHCompositeNode *startNode);
  int WriteNode(PHCompositeNode *thisNode);

 private:
  PHNodeIOManager *dstOut = nullptr;
  int m_SaveRunNodeFlag = 1;
  std::set<std::string> savenodes;
  std::set<std::string> saverunnodes;
  std::set<std::string> stripnodes;
  std::set<std::string> striprunnodes;
};

#endif
