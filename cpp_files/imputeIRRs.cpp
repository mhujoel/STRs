// static executable:
// g++ -O3 -Wall -static-libgcc -static-libstdc++ imputeIRRs.cpp -o imputeIRRs -I/n/groups/price/poru/HSPH_SVN/src/EAGLE -I/home/pl88/boost_1_58_0/install/include -L/n/groups/price/poru/external_software/libstdc++/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -L/n/groups/price/poru/external_software/zlib/zlib-1.2.11 -L/home/pl88/boost_1_58_0/install/lib -Wl,-Bstatic -lboost_iostreams -lz

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// files from https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz
#include "Types.hpp"
#include "FileUtils.cpp"
#include "StringUtils.cpp"

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage:" << endl;
    cerr << "- arg1 = IRR data to phase and impute (col1=ID ,col2=IRRs; no header)" << endl;
    cerr << "- arg2 = IBD neighbors (from extractNeighborsAll_RAP)" << endl;
    cerr << "- arg3 = output file (automatically gzipped if .gz)" << endl;
    cerr << endl;
    return 1;
  }

  // read IRR data; create lookup table of IDs (IDtoInd)
  vector <double> IRRs;
  const int MAX_ID = 10000000;
  vector <int> IDs, IDtoInd(MAX_ID, -1);
  {
    FileUtils::AutoGzIfstream finIRRs; finIRRs.openOrExit(argv[1]);
    int ID; double IRR;
    while (finIRRs >> ID >> IRR) {
      assert(0 <= ID && ID < MAX_ID);
      IDtoInd[ID] = IRRs.size();
      IDs.push_back(ID);
      IRRs.push_back(IRR);
    }
    finIRRs.close();
  }
  int N = IRRs.size();
  cout << "Read IRR data for " << N << " samples" << endl;
  
  // read haplotype neighbor data
  const unsigned int minNbr = 1, maxNbr = 10;
  vector < vector <int> > hapNbrs(2*N);
  vector <double> hapIRRs(2*N, NAN);
  {
    FileUtils::AutoGzIfstream finNeighbors; finNeighbors.openOrExit(argv[2]);
    string line; getline(finNeighbors, line);
    int ID, hap, nbrInd; double cMlen, cMedge; int IDnbr, hapNbr;
    while (finNeighbors >> ID >> hap >> nbrInd >> cMlen >> cMedge >> IDnbr >> hapNbr) {
      if (ID < 0 || IDnbr < 0) continue;
      assert(ID < MAX_ID);
      assert(hap == 1 || hap == 2);
      assert(IDnbr < MAX_ID);
      if (IDtoInd[ID]!=-1 && IDtoInd[IDnbr]!=-1 && hapNbrs[2*IDtoInd[ID]+hap-1].size() < maxNbr)
	hapNbrs[2*IDtoInd[ID]+hap-1].push_back(2*IDtoInd[IDnbr]+hapNbr-1);
    }
    finNeighbors.close();
  }
  
  // initialize hapIRRs to half of IRRs
  double meanIRRs = 0;
  {
    int NtoPhase = 0, totNbrs = 0;
    for (int i = 0; i < N; i++) {
      if (hapNbrs[2*i].size() >= minNbr && hapNbrs[2*i+1].size() >= minNbr) {
	hapIRRs[2*i] = hapIRRs[2*i+1] = IRRs[i]/2;
	NtoPhase++;
	totNbrs += hapNbrs[2*i].size() + hapNbrs[2*i+1].size();
	meanIRRs += IRRs[i];
      }
    }
    cout << "Read hap neighbors; phasing " << NtoPhase << " samples with >=" << minNbr
	 << " neighbors for both haps" << endl;
    cout << "Average number of neighbors to use per hap: " << totNbrs / (2.0*NtoPhase) << endl;
    meanIRRs /= NtoPhase;
  }

  // phase IRRs; impute and write output in last iter
  FileUtils::AutoGzOfstream fout; fout.openOrExit(argv[3]);
  fout << "ID\tIRRs\thap1phased\thap2phased\thap1imp\thap2imp" << endl;
  fout << std::fixed << std::setprecision(2);
  int nIters = 20;
  for (int iter = 0; iter < nIters; iter++) {
    for (int i = 0; i < N; i++) {
      double hapNbrSum[2] = {1e-9, 1e-9}; // compute mean IRRs among hap1's and hap2's neighbors
      int hapNbrNum[2] = {0, 0};
      double hapNbrMean[2];
      for (int h = 0; h < 2; h++) {
	const vector <int> &nbrs = hapNbrs[2*i+h];
	for (int k = 0; k < (int) nbrs.size(); k++)
	  if (!isnan(hapIRRs[nbrs[k]])) {
	    hapNbrSum[h] += hapIRRs[nbrs[k]];
	    hapNbrNum[h]++;
	  }
	if (hapNbrNum[h])
	  hapNbrMean[h] = hapNbrSum[h] / hapNbrNum[h];
	else
	  hapNbrMean[h] = meanIRRs/2; // mean-impute if no available neighbors
      }
      for (int h = 0; h < 2; h++) // phase sample i's IRRs: partition by neighbor means
	if (!isnan(hapIRRs[2*i+h]))
	  hapIRRs[2*i+h] = IRRs[i] * hapNbrMean[h] / (hapNbrMean[0] + hapNbrMean[1]);
      if (iter == nIters-1)
	fout << IDs[i] << "\t" << IRRs[i] << "\t" << hapIRRs[2*i] << "\t" << hapIRRs[2*i+1]
	     << "\t" << hapNbrMean[0] << "\t" << hapNbrMean[1] << endl;
    }
  }
  fout.close();

  return 0;
}
