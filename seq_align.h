/*
 PUIalign - Phylogenetically Unambiguous Indel Alignment
 
 Created by John P. McCrow - 9/5/2007
 
 Citation:
 McCrow JP, Alignment of phylogenetically unambiguous indels in Shewanella,
 Journal of Computational Biology 16 (11), 1517-1528 (2009)
*/
using namespace std;

#ifndef __SGI_STL_STRING
#include <string>
#endif

#ifndef __FSTREAM__
#include <fstream>
#endif

#ifndef __SSTREAM__
#include <sstream>
#endif

#ifndef __IOSTREAM__
#include <iostream>
#endif

#ifndef __SGI_STL_VECTOR
#include <vector>
#endif

#ifndef __SGI_STL_DEQUE
#include <deque>
#endif

#ifndef __SGI_STL_MAP
#include <map>
#endif

#include <limits.h>

#define _altype_NONE 0
#define _altype_MATCH 1
#define _altype_GAPi 2
#define _altype_GAPj 3

struct intpair {
  int i;
  int j;

  intpair() {i=j=-1;}
  intpair(int a, int b) {i=a; j=b;}
};

struct alignedgap {
  int gapseqid;
  int insseqid;
  int insposl;
  int insposr;
  int gappos;
  double bestscore;
  double gapscore;

  alignedgap() { reset(); }
  void reset(void);
  bool isgood(void);
  void copy(alignedgap *ag);
};

struct pairwise_al {
  double score;
  int seq1_id;
  int seq2_id;
  int seq1_offset;
  int seq2_offset;
  string seq1;
  string seq2;
  deque< intpair > aligned;
  alignedgap firstgap;

  pairwise_al() {score = 0.0; }
  void addaligned(int i, int j);
  bool hasgap(void);
  void print();
  int findgappos(int fseq, int a, int b);
};

struct pairwise_indel {
  int seqindex;
  int posl;
  int posr;

  pairwise_indel() {seqindex = posl = posr = -1; }
  bool isgood(void);
  void copy(pairwise_indel *pi);
};

class align_score {
 protected:
  map< string, int > M;
  double matchscore;
  double mismatchscore;
  double gapopenscore;
  double gapextscore;
 public:
  double match(char a, char b);
  void read_score_matrix(string filename);
  void set_static_match_scores(double match, double mismatch);
  void set_gap_scores(double open, double extend);
  double gapopen(void) {return gapopenscore;}
  double gapext(void) {return gapextscore;}
  double minusInfty(void) {return -LONG_MAX; }
};

class align {
 protected:
  int seq1_offset;
  int seq2_offset;
  string seq1;
  string seq2;
  align_score *scoring;
  vector<vector< double > > S;
  vector<vector< double > > E;
  vector<vector< double > > F;
  vector<vector< int > > fromi;
  vector<vector< int > > fromj;
  int maxi;
  int maxj;
  int n;
  int m;
  double maxS;
  pairwise_al al;

  double max(double a, double b) {return (b > a) ? b : a;}
  double max(double a, double b, double c) {return (b > a) ? max(b,c) : max(a,c);}
  void traceback(void);
  
 public:
  align();
  align(string s1, string s2, align_score *sc);
  void set_seqs(string s1, string s2);
  void set_scoring(align_score *s) {scoring=s;}
  void print(void);
  double get_score(void) {return maxS;}
  pairwise_al *get_alignment(void) {return &al;}
  
};

class fit_align : public align {
 protected:
  double fit(int fseq, int a, int b, int os);
  double refit(int fseq, int a, int b, int os, int diff_beg, int diff_end);

 public:
  fit_align() {}
  fit_align(string s1, string s2, align_score *sc);
  double fit();
  double fit_match(int fseq, int a, int b);
  double fit_gap(int fseq, int a, int b, int outer_size);
  double refit_match(int fseq, int a, int b, int d1, int d2);
  double refit_gap(int fseq, int a, int b, int outer_size, int d1, int d2);
};

align::align() {
}

align::align(string s1, string s2, align_score *sc) {
  set_seqs(s1,s2);
  set_scoring(sc);
}

fit_align::fit_align(string s1, string s2, align_score *sc) {
  set_seqs(s1,s2);
  set_scoring(sc);
}

