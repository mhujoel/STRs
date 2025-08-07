// g++ -O3 -Wall find_repeats_chr2_232_4Mb.cpp -o find_repeats_chr2_232_4Mb

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
    cerr << "Usage: arg1 = ID, arg2 = read_len; stdin = reads (SAM fields 1-11)" << endl;
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
	  cum_score -= 2;//(qual[i]=='?') ? 2 : 1;
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

      // look for left flank: GTTGGG (note: can be off by 1)
      if (TC_start == 2 && hamming(seq, 0, "GGTTTC") == 0)
	left_offset = 0;
      else if (hamming(seq, TC_start-3, "GGG") == 0 || hamming(seq, TC_start-5, "TTGGG") <= 1)
	left_offset = 0;
      else {
	int best_dist = 999;
	for (int offset = -2; offset <= -2; offset++) { // tried -10..1, but only -2 important
	  if (offset == 0) continue;
	  int dist = hamming(seq, TC_start-6 + offset, "GTTGGGTTTCTTTCTT", 6-offset);
	  if (dist <= (6-offset+1)/4 && dist <= best_dist) {
	    left_offset = offset;
	    best_dist = dist;
	  }
	}
	if (left_offset == NA && hamming(seq, TC_start-8, "AGGTTGGG") <= 2)
	  left_offset = 0;
      }

      // look for right flank: CTTTTTTTTGCTTCTTTGC (note: will match past last repeat)
      if ((TC_end == read_len-6 && hamming(seq, TC_end-4, "CTTTTTTTTG") == 0) ||
	  (TC_end == read_len-7 && hamming(seq, TC_end-4, "CTTTTTTTTGC") == 0))
	right_offset = 0;
      else if (hamming(seq, TC_end+1, "TTTT") <= 1 &&
	  (hamming(seq, TC_end+5, "GCT") == 0 || hamming(seq, TC_end+5, "GCTTC") <= 1))
	right_offset = 0;
      else {
	const char right_flank[] = "TCTTTCTTTTTTTTGCTTCTTG";
	for (int offset = 0; offset <= 2; offset++) // tried 0..10; error modes at 4 and 8
	  if (hamming(seq, TC_end+1, right_flank+10-offset) <= (strlen(right_flank+10-offset))/4)
	    right_offset = offset;
      }
    }

    if (left_offset != NA && right_offset != NA) {
      augment_span_freqs(qname_span_set_TC, freqsTC, qname,
			 (TC_end+right_offset) - (TC_start+left_offset) + 1);
    }
    else if (left_offset != NA || right_offset != NA) {
      int c_flank;
      if (left_offset != NA)
	c_flank = (TC_end+0) - (TC_start+left_offset) + 1;
      else // (right_offset != NA)
	c_flank = (TC_end+right_offset) - (TC_start+0) + 1;
      for (int c = 0; c <= c_flank; c++)
	freqsTC_flank[c]++;

    }
  }

  const int c_max = 29; // maximum allele length to attempt to genotype

  // find shortest allele with 3+ spanning reads (or first with 2+ and length>=21)
  int c_short = NA;
  for (int c = 1; c <= 4*c_max; c++)
    if (freqsTC[c] >= 3 || (freqsTC[c] >= 2 && c >= 4*21)) {
      c_short = c;
      break;
    }
  
  if (c_short != NA) { // short allele found; look for long allele
    int c_long = NA;
    // look for longer allele with different length mod 4: if found, select it
    for (int c = c_short+1; c <= 4*c_max; c++)
      if ((c-c_short)%4 != 0 && freqsTC[c] >= 2) {
	c_long = c;
	break;
      }
    if (c_long == NA) { // look for longer allele with same length mod 4
      if (freqsTC[c_short+4] >= 3) { // next allele has 3+ reads: germline or somatic?
	if ((c_short<=4*19 && 2*freqsTC[c_short+4]>=freqsTC[c_short]) ||
	    (c_short==4*20 && 3*freqsTC[c_short+4]>=2*freqsTC[c_short]) || 
	    (c_short==4*21 && 4*freqsTC[c_short+4]>=3*freqsTC[c_short]))
	  c_long = c_short+4;
      }
      else // next allele has 0,1,2 reads; look for subsequent allele with 2+ reads
	for (int c = c_short+8; c <= 4*c_max; c += 4)
	  if (freqsTC[c] >= 2) {
	    if ((c_short<=4*19 && freqsTC[c] > freqsTC[c-4]) ||
		(c_short<=4*21 && freqsTC[c] > freqsTC[c-4]+1) ||
		freqsTC[c] > freqsTC[c-4] + freqsTC[c-8]) 
	      c_long = c;
	    break;
	  }
    }

    if (c_long != NA) { // short and long allele found; output germline and somatic genotypes
      //cout << ID << " GERMLINE GENO: " << c_short/4.0 << " " << c_long/4.0 << endl;
      for (int p = 0; p < 2; p++) {
	int sum_diffs = 0, num_diffs = 0, neg_diffs = 0;
	if ((c_long-c_short)/4 >= 3 || (c_long-c_short)%4!=0) {
	  num_diffs += freqsTC[c_long-4];
	  neg_diffs += freqsTC[c_long-4];
	}
	int k_last = 0;
	for (int k = 0; k <= 7 && k <= k_last+3; k++) {
	  sum_diffs += freqsTC[c_long+4*k] * k;
	  num_diffs += freqsTC[c_long+4*k];
	  if (freqsTC[c_long+4*k])
	    k_last = k;
	}
	if (num_diffs >= 3)
	  cout << ID << " SOMATIC PHENO: " << c_long/4.0 << " +" << sum_diffs / (double) num_diffs
	       << endl;
	if ((c_long-c_short)%4 != 0)
	  swap(c_short, c_long); // also output somatic data for short allele if different mod 4
	else
	  break;
      }
    }
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
