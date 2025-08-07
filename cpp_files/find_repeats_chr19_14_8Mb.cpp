// g++ -O3 -Wall find_repeats_chr19_14_8Mb.cpp -o find_repeats_chr19_14_8Mb

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <set>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <cassert>

using namespace std;

const char dir_strs[2][4] = {"fwd", "rev"};
const int NA = -999;

void augment_span_freqs(set < pair <string, int> > &qname_span_set, vector <int> &freqs,
			const string &qname, int c) {
  if (qname_span_set.find(make_pair(qname, c)) == qname_span_set.end()) {
    qname_span_set.insert(make_pair(qname, c)); // don't double count if mate already counted
    freqs[c]++;
  }
}

int hamming(const string &str1, int start1, const char *str2, int len2) {
  int len1 = str1.length(), dist = 0;
  for (int i = 0; i < len2; i++)
    dist += i+start1 >= 0 && i+start1 < len1 ? str1[i+start1]!=str2[i] : 1;
  return dist;
}

int hamming(const string &str1, int start1, const char *str2) {
  return hamming(str1, start1, str2, strlen(str2));
}

int main(int argc, char *argv[]) {

  if (argc != 3) {
    cerr << "Usage: arg1 = ID, arg2 = read length; stdin = reads (SAM fields 1-11)" << endl;
    return 1;
  }

  string ID = argv[1];
  int read_len; sscanf(argv[2], "%d", &read_len);

  char revComp[256];
  for (int i = 0; i < 256; i++) revComp[i] = 'N';
  revComp['A'] = 'T';
  revComp['C'] = 'G';
  revComp['G'] = 'C';
  revComp['T'] = 'A';

  string qname; int flag; string rname; int pos, mapq; string cigar, rnext; int pnext, tlen;
  string seq, qual;
  set < pair <string, int> > qname_span_set_TC;
  vector <int> freqsTC(read_len+1), freqsTC_flank(read_len+1);

  while (cin >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen
	 >> seq >> qual) {
    int len = seq.length();
    if (len != read_len || (int) qual.length() != read_len) continue;
    int is_rev = (flag>>4)&1;
    if (seq.find("AAAGAAAG") != string::npos) { // reverse complement
      reverse(qual.begin(), qual.end());
      char seq_rc[read_len+1]; seq_rc[read_len] = '\0';
      for (int i = 0; i < read_len; i++)
	seq_rc[i] = revComp[(int) seq[read_len-1-i]];
      seq = seq_rc;
      is_rev = !is_rev;
    }

    // parse cigar string to check whether read supports REF allele
    {
      int pos_cur = pos;
      istringstream iss(cigar);
      int op_len; char op;
      bool is_REF_allele = false;
      while (iss >> op_len >> op) {
	if (op == 'S' || op == 'I') // don't consume ref sequence
	  ;
	else if (op == 'D' || op == 'M') { // consume ref sequence
	  if (op=='M' && rname=="chr19" && pos_cur<=14774670 && pos_cur+op_len>14774694)
	    is_REF_allele = true; // read spans repeat in GRCh38
	  pos_cur += op_len;
	}
	else // bad cigar op; primary alignment should only have M, S, I/D
	  break;
      }

      const int c_REF = 4*3; // REF allele TTTCTTTCTTTCTTT has 3 complete repeat units
      if (is_REF_allele) {
	augment_span_freqs(qname_span_set_TC, freqsTC, qname, c_REF);
	continue;
      }
    }

    // find longest mostly-TC, mostly-4bp-repeat stretch
    int TC_start = 0, TC_end = 0;
    int left_offset = NA, right_offset = NA;

    {
      int best_diff = 0, cum_score = 0, worst_score = 0, worst_ind = -1;
      for (int i = 0; i+4 < len; i++) {
	if ((seq[i]=='T' || seq[i]=='C') && seq[i]==seq[i+4]) {
	  cum_score++;
	  if (cum_score - worst_score >= best_diff) {
	    best_diff = cum_score - worst_score;
	    TC_start = worst_ind+1;
	    TC_end = i;
	  }
	}
	else {
	  cum_score -= 2;
	  if (cum_score < worst_score) {
	    worst_score = cum_score;
	    worst_ind = i;
	  }
	}
      }
      int numT = 0, numC = 0, TC_len = TC_end-TC_start+1;
      for (int i = TC_start; i <= TC_end; i++) {
	numT += seq[i]=='T';
	numC += seq[i]=='C';
      }

      if (!(TC_len>10 && numT+numC >= 0.9*TC_len && 3*numC <= TC_len && TC_len <= 5*numC))
	continue; // require >10bp, 90% TC content, 1/3 to 1/5 C content

      if (TC_len > read_len-4 - 5) { // read supports a long allele
	for (int c = 0; c <= read_len; c++)
	  freqsTC_flank[c]++;
	continue;
      }

      // look for left flank: TTAAAT (note: can be off by 1)
      if (TC_start == 2 && hamming(seq, 0, "ATTTTC") == 0)
	left_offset = 0;
      else if (hamming(seq, TC_start-3, "AAT") == 0)
	left_offset = 0;
      else if (hamming(seq, TC_start-5, "TAAAT") <= 1)
	left_offset = 0;

      // look for right flank: CTTTTTGAGACA (note: last repeat unit won't have been matched)
      if (hamming(seq, TC_end+1, "CTTTT") <= 1 && hamming(seq, TC_end+6, "TGA") == 0)
	right_offset = 0;
    }

    if (left_offset != NA && right_offset != NA) {
      augment_span_freqs(qname_span_set_TC, freqsTC, qname,
			 (TC_end+right_offset) - (TC_start+left_offset) + 2);
    }
    else if (left_offset != NA || right_offset != NA) {
      int c_flank;
      if (left_offset != NA)
	c_flank = (TC_end+0) - (TC_start+left_offset) + 2;
      else // (right_offset != NA)
	c_flank = (TC_end+right_offset) - (TC_start+0) + 2;
      for (int c = 0; c <= c_flank; c++)
	freqsTC_flank[c]++;
    }
  }

  const int c_min = 4*17; // minimum allele length to attempt to genotype
  const int c_max = 4*29; // maximum allele length to attempt to genotype

  // find shortest allele in 17..29 with 2+ spanning reads
  // require exactly one of:
  // - 3+ reads with length >=len+7 (spanning or flank)
  // - 3+ reads with length 3..16
  int c_mid_guess = NA;
  for (int c_mid = c_min; c_mid <= c_max; c_mid += 4)
    if (freqsTC[c_mid] >= 2 // found shortest allele in 17..29 with 2+ spanning reads
	|| (c_mid >= 4*21 && freqsTC[c_mid]==1 && freqsTC[c_mid+4]>=1 && freqsTC[c_mid+4]<=2)) {
      int num_other_alleles_found = 0;
      int short_or_interrupted_allele = NA;
      for (int c = 4*3; c < c_max; c++)
	if ((c < c_min || c%4 != 0) && freqsTC[c] >= 3) {
	  num_other_alleles_found++;
	  short_or_interrupted_allele = c;
	}
      int freq_long = freqsTC_flank[c_mid+4*7]; // count flanking reads supporting len >= c_mid+7
      for (int c = c_mid+4*7; c <= read_len; c++) // include spanning reads of length >= c_mid+7
	freq_long += freqsTC[c];
      if (freq_long >= 3)
	num_other_alleles_found++;
      if (num_other_alleles_found == 1) {
	int sum_diffs = 0, num_diffs = 0, neg_diffs = 0;
	if (c_mid-4 != short_or_interrupted_allele) {
	  num_diffs += freqsTC[c_mid-4];
	  neg_diffs += freqsTC[c_mid-4];
	}
	int k_last = 0;
	for (int k = 0; k <= 7 && k <= k_last+3; k++) {
	  sum_diffs += freqsTC[c_mid+4*k] * k;
	  num_diffs += freqsTC[c_mid+4*k];
	  if (freqsTC[c_mid+4*k])
	    k_last = k;
	}
	if (num_diffs >= 3)
	  cout << ID << " SOMATIC PHENO: " << c_mid/4.0 << " +" << sum_diffs / (double) num_diffs
	       << endl;
      }
      c_mid_guess = c_mid;
      break;
    }

  // output likely germline genotype (for confident hets)
  int num_short_or_interrupted_alleles_found = 0, A1 = NA, A2 = NA;
  for (int c = 4*3; c < c_max; c++)
    if ((c < c_min || c%4 != 0) && freqsTC[c] >= 3) {
      num_short_or_interrupted_alleles_found++;
      if (A1 == NA)
	A1 = c;
      A2 = c;
    }
  const int c_long_flank = 4*35;
  int long_allele_found = false;
  if (freqsTC_flank[c_long_flank] >= 3) {
    A2 = c_long_flank;
    long_allele_found = true;
  }
  int mid_allele_found = c_mid_guess != NA;
  if (num_short_or_interrupted_alleles_found + mid_allele_found + long_allele_found == 2) {
    if (mid_allele_found) {
      if (long_allele_found)
	A1 = c_mid_guess;
      else
	A2 = c_mid_guess;
    }
    //cout << ID << " GERMLINE GENO: " << A1/4.0 << " " << A2/4.0 << endl;
  }
  /*
  // summarize spanned alleles
  cout << ID << " spanned alleles:";
  for (int c = 0; c <= read_len; c++)
    if (freqsTC[c])
      cout << " " << c/4.0 << "(" << freqsTC[c] << ")";
  cout << endl;
  */
  return 0;
}
