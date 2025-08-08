// g++ -O2 -fopenmp -Wall computeIBSpbwt.cpp -o computeIBSpbwt -I/n/groups/price/poru/HSPH_SVN/src/EAGLE -I/home/pl88/boost_1_58_0/install/include -L/home/pl88/boost_1_58_0/install/lib -Wl,-Bstatic -lboost_iostreams -Wl,-Bdynamic -lz
// g++ -O2 -fopenmp -Wall -static-libgcc -static-libstdc++ computeIBSpbwt.cpp -o computeIBSpbwt -I/n/groups/price/poru/HSPH_SVN/src/EAGLE -I/home/pl88/boost_1_58_0/install/include -L/n/groups/price/poru/external_software/libstdc++/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -L/n/groups/price/poru/external_software/zlib/zlib-1.2.11 -L/home/pl88/boost_1_58_0/install/lib -Wl,-Bstatic -lboost_iostreams -lz

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "omp.h"
#include "zlib.h"

// files from https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz
#include "Types.hpp"
#include "FileUtils.cpp"
#include "StringUtils.cpp"
#include "MemoryUtils.cpp"
#include "Timer.cpp"
#include "MapInterpolater.cpp"

using namespace std;

class HapBitsT {
  uint64 *haploBitsT;
  uint64 Nhaps, M, M64;
public:
  HapBitsT(uint64 _Nhaps, uint64 _M);
  ~HapBitsT(void);
  void setBit(uint64 n, uint64 m);
  void flipBit(uint64 n, uint64 m);
  int getBit(uint64 n, uint64 m) const;
  uint64 getNhaps(void) const;
  uint64 getM(void) const;
  uint64 getHaploBitsSize(void) const;
  uint64* getHaploBitsT(void) const;
  uint64* getHaploBitsTrow(uint64 n) const;
};

void HapBitsT::setBit(uint64 n, uint64 m) { haploBitsT[n*M64 + (m>>6)] |= 1ULL<<(m&63); }

void HapBitsT::flipBit(uint64 n, uint64 m) { haploBitsT[n*M64 + (m>>6)] ^= 1ULL<<(m&63); }

int HapBitsT::getBit(uint64 n, uint64 m) const { return (haploBitsT[n*M64 + (m>>6)]>>(m&63))&1; }

uint64 HapBitsT::getNhaps(void) const { return Nhaps; }

uint64 HapBitsT::getM(void) const { return M; }

uint64 HapBitsT::getHaploBitsSize(void) const { return Nhaps * M64; }

uint64* HapBitsT::getHaploBitsT(void) const { return haploBitsT; }

uint64* HapBitsT::getHaploBitsTrow(uint64 n) const { return haploBitsT + n*M64; }

HapBitsT::HapBitsT(uint64 _Nhaps, uint64 _M) {
  Nhaps = _Nhaps;
  M = _M;
  M64 = (M+63)/64;
  haploBitsT = ALIGNED_MALLOC_UINT64S(Nhaps * M64);
  memset(haploBitsT, 0, Nhaps * M64 * sizeof(haploBitsT[0]));
}

HapBitsT::~HapBitsT(void) { ALIGNED_FREE(haploBitsT); }

void pbwt(const HapBitsT &hapBitsT, bool isFwd, int mClosest, int *lexSorts[2]) {

  int H = hapBitsT.getNhaps(), M = hapBitsT.getM();

  // initialize work arrays
  int *a1 = new int[H], *d1 = new int[H], *a = new int[H], *b = new int[H], *d = new int[H],
    *e = new int[H];
  for (int n = 0; n < H; n++) {
    a1[n] = n;
    d1[n] = M;
  }

  /***** RUN PBWT *****/
  for (int m = M-1; m >= 0; m--) {
    const int mBit = isFwd ? M-1-m : m; // process SNPs in forward or reverse order
    // compute sort order and divergence array
    int u = 0, v = 0, p = m, q = m;
    for (int i = 0; i < H; i++) {
      if (d1[i] < p) p = d1[i];
      if (d1[i] < q) q = d1[i];
      if (hapBitsT.getBit(a1[i], mBit) == 0) {
	a[u] = a1[i]; d[u] = p; u++; p = M;
      }
      else {
	b[v] = a1[i]; e[v] = q; v++; q = M;
      }
    }
    memcpy(a1, a, u * sizeof(a1[0])); memcpy(a1+u, b, v * sizeof(a1[0]));
    memcpy(d1, d, u * sizeof(d1[0])); memcpy(d1+u, e, v * sizeof(d1[0]));
    
    if (mBit == mClosest) {
      // store sort order at specified SNPs
      lexSorts[isFwd] = new int[H];
      memcpy(lexSorts[isFwd], a1, H * sizeof(int));
      break; // early exit: reached focal SNP
    }
  }
  delete[] a1;
  delete[] d1;
  delete[] a;
  delete[] b;
  delete[] d;
  delete[] e;
}


const int MAX_ERR = 2;

struct IBSmatch {
  int hap;
  double cMerrs[2][MAX_ERR];
  double len;
  bool operator < (const IBSmatch &match2) const {
    if (len != match2.len) return len > match2.len;
    else return hap < match2.hap;
  }
};

IBSmatch computeIBS(const HapBitsT &hapBitsT, int h1, int h2, int m0, const vector <double> &cMs) {
  int M = hapBitsT.getM();

  IBSmatch matchInfo; matchInfo.hap = h2;

  const double weight0err = 0.5, weightTotalLen = 0.25;

  const uint64 *haploBitsRow1 = hapBitsT.getHaploBitsTrow(h1);
  const uint64 *haploBitsRow2 = hapBitsT.getHaploBitsTrow(h2);

  for (int k = 0; k < 2; k++) {
    int dir = k==0 ? -1 : 1;
    int mEnd = dir==1 ? M : -1;
    int errs = 0;
    uint64 diffBits = 0;
    for (int m = m0; m >= 0 && m < M; m += dir) {
      // fill 64-bit xor buffer at m0 or when entering new block: m%64==0 (fwd) or m%64=63 (rev)
      if (m == m0 || (dir==1 && (m&63)==0) || (dir==-1 && (m&63)==63)) {
	diffBits = haploBitsRow1[m>>6] ^ haploBitsRow2[m>>6];
	if (diffBits == 0) { // all bits match in this block of 64 SNPs
	  if (dir == 1) m |= 63; // skip ahead to next 64-SNP block
	  else m &= ~63; // skip back to previous 64-SNP block
	  continue;
	}
      }

      //if (hapBitsT.getBit(h1, m) != hapBitsT.getBit(h2, m)) {
      if ((diffBits>>(m&63))&1) {
	matchInfo.cMerrs[k][errs] = cMs[m] - cMs[m0];
	errs++;
	if (errs == MAX_ERR)
	  break;
      }
    }
    while (errs < MAX_ERR) {
      matchInfo.cMerrs[k][errs] = cMs[mEnd-dir] - cMs[m0];
      errs++;
    }
  }  
  double lenLeft = weight0err*-matchInfo.cMerrs[0][0] + (1-weight0err)*-matchInfo.cMerrs[0][1];
  double lenRight = weight0err*matchInfo.cMerrs[1][0] + (1-weight0err)*matchInfo.cMerrs[1][1];
  double lenTot = lenLeft + lenRight;
  double lenMin = min(lenLeft, lenRight);
  // randomize for sorting
  matchInfo.len = weightTotalLen*lenTot + (1-weightTotalLen)*lenMin + (h1*h2%1000003)*1e-12;

  return matchInfo;
}

vector <string> processSample(const char *sampleFile) {
  vector <string> IDs;
  FileUtils::AutoGzIfstream fin;
  fin.openOrExit(sampleFile);
  string line; getline(fin, line); getline(fin, line);
  string ID_1, ID_2;
  while (fin >> ID_1 >> ID_2) {
    IDs.push_back(ID_1);
    getline(fin, line);
  }
  cout << "Read " << IDs.size() << " IDs from sample file" << endl << endl;
  return IDs;
}

int main(int argc, char *argv[]) {

  if (argc != 9) {
    cerr << "Usage:" << endl;
    cerr << "- arg1 = chr (integer)" << endl;
    cerr << "- arg2 = focal bp (hg19 for UKB SNP-array haplotypes)" << endl;
    cerr << "- arg3 = bgen file" << endl;
    cerr << "- arg4 = sample file" << endl;
    cerr << "- arg5 = genetic map file (in same build as bgen file)" << endl;
    cerr << "- arg6 = number of IBS neighbors per haplotype" << endl;
    cerr << "- arg7 = threads" << endl;
    cerr << "- arg8 = output file" << endl;
    exit(1);
  }
  int chr; assert(sscanf(argv[1], "%d", &chr)==1);
  int bpTarget; assert(sscanf(argv[2], "%d", &bpTarget)==1);
  const char *bgenFile = argv[3];
  const char *sampleFile = argv[4];
  const char *geneticMapFile = argv[5];
  int numNbrs; assert(sscanf(argv[6], "%d", &numNbrs)==1);
  int threads; assert(sscanf(argv[7], "%d", &threads)==1);
  const char *outFile = argv[8];

  const int pbwtBand = 10*numNbrs;

  Timer timer;
  
  cout << "Setting number of threads to " << threads << endl << endl;
  omp_set_num_threads(threads);

  vector <string> IDs = processSample(sampleFile);

  Genetics::MapInterpolater mapInterpolater(geneticMapFile);

  FILE *fin = fopen(bgenFile, "rb"); assert(fin != NULL);

  /********** read BGEN header **********/
  uint Mphased;
  {
    uint offset; fread(&offset, 4, 1, fin);
    uint L_H; fread(&L_H, 4, 1, fin);
    uint Mbgen; fread(&Mbgen, 4, 1, fin); cout << "# variants in bgen file: " << Mbgen << endl;
    assert(Mbgen != 0);
    Mphased = Mbgen;
    uint Nbgen; fread(&Nbgen, 4, 1, fin); cout << "# samples in bgen file: " << Nbgen << endl;
    if (Nbgen != IDs.size()) {
      cerr << "ERROR: Number of samples in BGEN header does not match sample file" << endl;
      exit(1);
    }
    char magic[5]; fread(magic, 1, 4, fin); magic[4] = '\0';
    fseek(fin, L_H-20, SEEK_CUR);
    uint flags; fread(&flags, 4, 1, fin);
    uint CompressedSNPBlocks = flags&3; assert(CompressedSNPBlocks==1);
    uint Layout = (flags>>2)&0xf; assert(Layout==2);
    fseek(fin, offset+4, SEEK_SET);
  }

  int H = 2*IDs.size();
  HapBitsT hapBitsT(H, Mphased);
  vector <double> cMvec;
  int mClosest = 0;

  /********** read variant blocks **********/
  {
    int bpClosest = 1<<30;

    char snpID[65536], rsID[65536], chrStr[65536];
    char *allele1, *allele0;
    uint maxLA = 65536, maxLB = 65536;
    allele1 = (char *) malloc(maxLA+1);
    allele0 = (char *) malloc(maxLB+1);

    vector <uchar> zBuf, buf; uint zBufLen, bufLen;

    for (uint mbgen = 0; mbgen < Mphased; mbgen++) {
      ushort LS; fread(&LS, 2, 1, fin);
      fread(snpID, 1, LS, fin); snpID[LS] = '\0';
      ushort LR; fread(&LR, 2, 1, fin);
      fread(rsID, 1, LR, fin); rsID[LR] = '\0';
      ushort LC; fread(&LC, 2, 1, fin);
      fread(chrStr, 1, LC, fin); chrStr[LC] = '\0';

      int bp; fread(&bp, 4, 1, fin);
      if (abs(bp - bpTarget) < abs(bpClosest - bpTarget)) {
	mClosest = cMvec.size();
	bpClosest = bp;
      }
      cMvec.push_back(100*mapInterpolater.interp(chr, bp));

      ushort K; fread(&K, 2, 1, fin); assert(K==2);
      uint LA; fread(&LA, 4, 1, fin);
      if (LA > maxLA) {
	maxLA = 2*LA;
	free(allele1);
	allele1 = (char *) malloc(maxLA+1);
      }
      fread(allele1, 1, LA, fin); allele1[LA] = '\0';
      uint LB; fread(&LB, 4, 1, fin);
      if (LB > maxLB) {
	maxLB = 2*LB;
	free(allele0);
	allele0 = (char *) malloc(maxLB+1);
      }
      fread(allele0, 1, LB, fin); allele0[LB] = '\0';

      uint C; fread(&C, 4, 1, fin);
      if (C > zBuf.size()) zBuf.resize(C-4);
      uint D; fread(&D, 4, 1, fin);
      if (D > buf.size()) buf.resize(D);
      zBufLen = C-4; bufLen = D;
      fread(&zBuf[0], 1, C-4, fin);

      uLongf destLen = bufLen;
      if (uncompress(&buf[0], &destLen, &zBuf[0], zBufLen) != Z_OK || destLen != bufLen) {
	cerr << "ERROR: uncompress() failed" << endl;
	exit(1);
      }
      uchar *bufAt = &buf[0];
      uint N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;
      assert(N == IDs.size());
      K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
      assert(K == 2U);
      uint Pmin = *bufAt; bufAt++;
      assert(Pmin == 2U);
      uint Pmax = *bufAt; bufAt++;
      assert(Pmax == 2U);
      for (uint i = 0; i < N; i++) {
	uint ploidyMiss = *bufAt; bufAt++;
	assert(ploidyMiss == 2U);
      }
      uint Phased = *bufAt; bufAt++;
      assert(Phased == 1U);
      uint bgenBits = *bufAt; bufAt++;
      assert(bgenBits == 16U);

      for (uint i = 0; i < N; i++)
	for (int h = 0; h < 2; h++) {
	  ushort hapProb = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
	  assert(hapProb == 0 || hapProb == 65535U);
	  if (hapProb)
	    hapBitsT.setBit(2*i+h, mbgen);
	}
    }

    free(allele0);
    free(allele1);
    assert(fin); assert(fgetc(fin) == EOF);
    fclose(fin);

    cout << "Found closest variant to " << bpTarget << " at " << bpClosest << endl;
  }

  cout << "\nTime for processing bgen file: " << timer.update_time() << " sec" << endl;

  // run PBWT; save lexicographic sorts at closest probe
  int *lexSorts[2];
#pragma omp parallel for
  for (int isFwd = 0; isFwd < 2; isFwd++)
    pbwt(hapBitsT, isFwd, mClosest, lexSorts);

  cout << "\nTime for PBWT: " << timer.update_time() << " sec\n" << endl;

  // open output file
  FileUtils::AutoGzOfstream fout;
  fout.openOrExit(outFile);
  fout << "ID\thap\tnbrInd\tcMlen\tcMedge\tIDnbr\thapNbr" << endl;

  /***** RUN NEIGHBOR-FINDING *****/

  // set up lookup: position of each haplotype in PBWT reverse and forward sorts

  int *revInd = new int[H];
  for (int i = 0; i < H; i++)
    revInd[lexSorts[0][i]] = i;
  int *fwdInd = new int[H];
  for (int i = 0; i < H; i++)
    fwdInd[lexSorts[1][i]] = i;

  const int hStep = 100000; // number of haplotypes to analyze before dumping a batch of output
  vector <string> outStrs(hStep);
  for (int hStart = 0; hStart < H; hStart += hStep) {
#pragma omp parallel for
    for (int h = hStart; h < min(H, hStart+hStep); h++) {
      set <IBSmatch> bestMatches;
      set <int> hNbrs0;

      const int i01[2] = {revInd[h], fwdInd[h]};
      
      for (int k = 0; k < 2; k++) {
	for (int j = max(0, i01[k]-pbwtBand); j <= min(H-1, i01[k]+pbwtBand); j++) {
	  if (j == i01[k]) continue;
	  const int hNbr = lexSorts[k][j];
	  if (k == 0) hNbrs0.insert(hNbr);
	  if (k == 1 && hNbrs0.count(hNbr)) continue; // don't retry neighbors already considered
	  bestMatches.insert(computeIBS(hapBitsT, h, hNbr, mClosest, cMvec));
	  if ((int) bestMatches.size() > numNbrs)
	    bestMatches.erase(--bestMatches.end());
	}
      }
      // store best matches
      ostringstream oss;
      //oss << IDs[h/2] << " " << h%2+1;
      oss << std::setprecision(2) << std::fixed;
      int nbrInd = 0;
      for (set <IBSmatch>::iterator it = bestMatches.begin(); it != bestMatches.end(); it++) {
	/*
	oss << " " << IDs[it->hap/2] << " " << it->hap%2+1
	    << " " << it->cMerrs[0][1] << " " << it->cMerrs[0][0]
	    << " " << it->cMerrs[1][0] << " " << it->cMerrs[1][1];
	*/
	double cMleft = -it->cMerrs[0][it->cMerrs[0][1] < it->cMerrs[0][0] - 1];
	double cMright = it->cMerrs[1][it->cMerrs[1][1] > it->cMerrs[1][0] + 1];
	nbrInd++;
	oss << IDs[h/2] << "\t" << h%2+1 << "\t" << nbrInd;
	oss << "\t" << cMleft + cMright;
	oss << "\t" << min(cMleft, cMright);
	oss << "\t" << IDs[it->hap/2];
	oss << "\t" << it->hap%2+1;
	oss << endl;
      }
      outStrs[h-hStart] = oss.str();
    }

    cout << "Finished batch " << hStart/hStep+1 << " of " << (H+hStep-1)/hStep << " in "
	 << timer.update_time() << " sec";

    for (int h = hStart; h < min(H, hStart+hStep); h++)
      fout << outStrs[h-hStart];
    
    cout << " + " << timer.update_time() << " sec writing" << endl;
  }

  fout.close();

  for (int isFwd = 0; isFwd < 2; isFwd++)
    delete[] lexSorts[isFwd];
  delete[] fwdInd;

  return 0;
}
