/*
 PUIalign - Phylogenetically Unambiguous Indel Alignment
 
 Created by John P. McCrow - 9/5/2007
 
 Citation:
 McCrow JP, Alignment of phylogenetically unambiguous indels in Shewanella,
 Journal of Computational Biology 16 (11), 1517-1528 (2009)
 */
using namespace std;

#include "seq_align.h"

void align_score::set_static_match_scores(double match, double mismatch) {
  matchscore = match;
  mismatchscore = mismatch;
}

void align_score::set_gap_scores(double open, double extend) {
  gapopenscore = open;
  gapextscore = extend;
}

void align_score::read_score_matrix(string filename) {
  string line;
  ifstream infile(filename.c_str());
  vector<string> colnames;
  string rowname;
  bool firstline;
  int col;

  firstline = true;
  while (infile.good() && getline(infile,line,'\n')) {  //Each line
    istringstream iss;
    iss.str(line);
    string v;
    col = 0;

    while (iss.good() && getline(iss, v, ',')) {  //Each comma separated value
      if(col == 0) {
	if(!firstline) rowname = v;
      } else {
	if(firstline) colnames.push_back(v);
	else if(rowname.length() && colnames[col-1].length()) {
	  M[rowname+colnames[col-1]] = atoi(v.c_str());  //Set M with each value
	}
      }
      col++;
    }

    firstline = false;
  }

  assert(M.size() > 1);
  
}

double align_score::match(char a, char b) {
  if(M.size() > 1) {
    //    cout<<(string (1,a))+(string (1,b))<<" = "<<M[(string (1,a))+(string (1,b))]<<endl;
    return M[(string (1,a))+(string (1,b))];
  } else {
    return (a == b) ? matchscore : mismatchscore;
  }
}

void alignedgap::reset(void) {
  gapseqid = insseqid = insposl = insposr = gappos = -1; 
}

bool alignedgap::isgood(void) {
  return (gapseqid >= 0 && 
	  insseqid >= 0 &&
	  insposl >= 0 &&
	  insposr >= 0 &&
	  gappos >= 0);
}

void alignedgap::copy(alignedgap *ag) {
  gapseqid = ag->gapseqid;
  insseqid = ag->insseqid;
  insposl = ag->insposl;
  insposr = ag->insposr;
  gappos = ag->gappos;
  bestscore = ag->bestscore;
  gapscore = ag->gapscore;
}

bool pairwise_indel::isgood(void) { 
  return seqindex >= 0 && posl >= 0 && posr >= 0; 
}

void pairwise_indel::copy(pairwise_indel *pi) {
  pi->seqindex = seqindex;
  pi->posl = posl;
  pi->posr = posr;
}

bool pairwise_al::hasgap(void) {
  int i;

  for(i=1; i<aligned.size(); i++) {
    if(aligned[i-1].i == aligned[i].i || aligned[i-1].j == aligned[i].j)
      return true;
  }
  return false;
}

void pairwise_al::print() {
  int i;
  string alignstri;
  string alignstrj;
  int n = seq1.length();

  /*
    cout<<seq1<<endl;
    cout<<seq2<<endl;

    for(i=0;i<segments.size(); i++) {
    cout<<"("<<segments[i].i<<","<<segments[i].j<<")"<<endl;
    }
  */
  
  for(i=0; i<aligned[0].i-1; i++) {
    alignstri += seq1[i];
    alignstrj += "-";
  }

  alignstri += seq1[aligned[0].i-1];
  alignstrj += seq2[aligned[0].j-1];  

  for(i=1; i<aligned.size(); i++) {
    if(aligned[i-1].i == aligned[i].i) {
      alignstri += "-";
    } else {
      alignstri += seq1[aligned[i].i-1];
    }
    if(aligned[i-1].j == aligned[i].j) {
      alignstrj += "-";
    } else {
      alignstrj += seq2[aligned[i].j-1];
    }
  }

  for(i=aligned[aligned.size()-1].i; i<n; i++) {
    alignstri += seq1[i];
    alignstrj += "-";
  }

  cout<<alignstri<<endl;
  cout<<alignstrj<<endl;
  
}

int pairwise_al::findgappos(int fseq, int a, int b) {
  int i;

  for(i=1; i<aligned.size(); i++) {
    if(fseq == 1 && seq1_offset + aligned[i].i >= a && seq1_offset + aligned[i].i <= b && 
       aligned[i].j == aligned[i-1].j) {
      return seq2_offset + aligned[i].j + 1;
    }
    if(fseq == 2 && seq2_offset + aligned[i].j >= a && seq2_offset + aligned[i].j <= b && 
       aligned[i].i == aligned[i-1].i) {
      return seq1_offset + aligned[i].i + 1;
    }
  }

  return -1;
}

void align::set_seqs(string s1, string s2) {
  seq1 = s1;
  seq2 = s2;
  n = seq1.length();
  m = seq2.length();
}

void align::print(void) {
  al.print();
}

void align::traceback(void) {
  int curi;
  int curj;
  int nexti;
  int nextj;

  al.aligned.clear();
  al.firstgap.reset();

  al.score = maxS;
  al.seq1 = seq1;
  al.seq2 = seq2;
  al.seq1_offset = seq1_offset;
  al.seq2_offset = seq2_offset;  

  curi = maxi;
  curj = maxj;

  bool followGapE = false;
  bool followGapF = false;

  while(curi > 0 && curj > 0) {
    if(followGapE) {
      if(E[curi][curj] == S[curi][curj-1] + scoring->gapopen()) {
	followGapE = false;
      }
      nexti = curi;
      nextj = curj-1;
    } else if(followGapF) {
      if(F[curi][curj] == S[curi-1][curj] + scoring->gapopen()) {
	followGapF = false;
      }
      nexti = curi-1;
      nextj = curj;
    } else if(S[curi][curj] == S[curi-1][curj-1] + scoring->match(seq1[curi-1], seq2[curj-1])) {
      nexti = curi-1;
      nextj = curj-1;
    } else if(S[curi][curj] == E[curi][curj]) {
      nexti = curi;
      nextj = curj-1;
      if(!al.firstgap.isgood()) {
	al.firstgap.insseqid = 
	al.firstgap.insposr = curj;
	al.firstgap.gappos = curi;	
      }
      if(E[curi][curj] == S[curi][curj-1] + scoring->gapopen()) {
      } else {
	followGapE = true;
      }
    } else if(S[curi][curj] == F[curi][curj]) {
      nexti = curi-1;
      nextj = curj;
      if(F[curi][curj] == S[curi-1][curj] + scoring->gapopen()) {
      } else {
	followGapF = true;
      }
    } else {
      cout<<"Corrupt Alignment"<<endl;
      return;
    }

    intpair ip(curi, curj);
    al.aligned.push_front(ip);
    
    curi = nexti;
    curj = nextj;	
  }

}

double fit_align::fit() {
  //cout<<"fit: -1,-1,-1,-1: ";
  return fit(-1, -1,-1,-1);
}

double fit_align::fit_gap(int fseq, int a, int b, int outer_size) {
  //cout<<"gap: "<<fseq<<","<<a<<","<<b<<","<<outer_size<<": ";
  return fit(fseq,a,b,outer_size);
}

double fit_align::fit_match(int fseq, int a, int b) {
  //cout<<"match: "<<fseq<<","<<a<<","<<b<<": ";
  return fit(fseq,a,b,-1);
}

double fit_align::refit_match(int fseq, int a, int b, int d1, int d2) {
  assert(S.size() > 1);
  return refit(fseq,a,b,-1,d1,d2);
}

double fit_align::refit_gap(int fseq, int a, int b, int outer_size, int d1, int d2) {
  assert(S.size() > 1);
  return refit(fseq,a,b,outer_size,d1,d2);
}

double fit_align::fit(int fseq, int a, int b, int os) {
  int i;
  int j;
  double s_match;
  
  bool force_match = (fseq == 1 || fseq == 2) && a > -1 && b > -1 && os < 0;
  bool force_gap = (fseq == 1 || fseq == 2) && a > -1 && b > -1 && os > -1;

  //cout<<force_match<<","<<force_gap<<","<<fseq<<","<<a<<","<<b<<","<<os<<endl;

  assert(force_match || force_gap || (fseq < 0 && a < 0 && b < 0 && os < 0));

  //Reset alignment matricies
  S.clear();
  E.clear();
  F.clear();
  fromi.clear();
  fromj.clear();

  //Initialize matricies
  for(i=0; i<=n; i++){
    vector<double> S_r;
    vector<double> E_r;
    vector<double> F_r;
    for(j=0; j<=m; j++){      
      S_r.push_back(0);
      E_r.push_back(0);
      F_r.push_back(0);
    }
    S.push_back(S_r);
    E.push_back(E_r);
    F.push_back(F_r);
  }
  
  //DP
  for(i=0; i<=n; i++) {
    for(j=0; j<=m; j++) {

      //Initialize row 0 and column 0
      if(i==0 || j==0) {
	E[i][j] = scoring->gapopen() - scoring->gapext();
	    
	if(i==0 && j>0) {
	  E[i][j] = E[i][j-1] + scoring->gapext();
	  S[i][j] = E[i][j];
	}
	
      } else {
	s_match = scoring->match(seq1[i-1], seq2[j-1]);

	if((force_gap //Forced Gap, boundry match
	    && ((fseq == 1 
		 && ((i >= a - os 
		      && i <= a - 1)
		     || (i >= a + b
			 && i <= a + b + os)))
		|| (fseq == 2 
		    && ((j >= a - os 
			 && j <= a - 1)
			|| (j >= a + b 
			    && j <= a + b + os)))))
	   //or Forced Match
	   || (force_match
	       && ((fseq == 1 && i >= a && i <= a + b + 1)
		   || (fseq == 2 && j >= a && j <= a + b + 1)))
	   ) {
		    
	  E[i][j] = scoring->minusInfty(); //S[i][j-1] + scoring->gapopen();
	  F[i][j] = scoring->minusInfty(); //S[i-1][j] + scoring->gapopen();
	  S[i][j] = S[i-1][j-1] + s_match;
		    
	} else if(force_gap && fseq == 1 //Forced Gap Seq 1 (Meaning Seq 2 has the deletion)
		  && i >= a && i <= a + b - 1) {
		    
	  E[i][j] = S[i][j] + scoring->gapopen();
	  F[i][j] = max(S[i-1][j] + scoring->gapopen(),
			F[i-1][j] + scoring->gapext());		    
	  S[i][j] = F[i][j];

	} else if(force_gap && fseq == 2 //Forced Gap Seq 2 (Meaning Seq 1 has the deletion)
		  && j >= a && j <= a + b - 1) {
		    
	  E[i][j] = max(S[i][j-1] + scoring->gapopen(),
			E[i][j-1] + scoring->gapext());
	  F[i][j] = S[i][j] + scoring->gapopen();
	  S[i][j] = E[i][j];

	} else { //No force
	  E[i][j] = max(S[i][j-1] + scoring->gapopen(), E[i][j-1] + scoring->gapext());
	  F[i][j] = max(S[i-1][j] + scoring->gapopen(), F[i-1][j] + scoring->gapext());
	  S[i][j] = max(S[i-1][j-1] + s_match,
			E[i][j],
			F[i][j]);

	}
      }

      //Scores should never be near minusInfty, unless E or F is set to it exactly
      assert(S[i][j] > scoring->minusInfty() + 100 &&
	     (E[i][j] > scoring->minusInfty() + 100 || E[i][j] == scoring->minusInfty()) &&
	     (E[i][j] > scoring->minusInfty() + 100 || E[i][j] == scoring->minusInfty()));
    } //j
  } //i

  //Find End of Fit Alignment
  maxS = S[n][m]-1;
  maxi = n;
  maxj = m;
  for(i=n; i>=m; i--) {
    if(S[i][m] > maxS) {
      maxS = S[i][m];
      maxi = i;
      maxj = m;
    }
  }

  traceback();

  return maxS;
}

double fit_align::refit(int fseq, int a, int b, int os, int diff_beg, int diff_end) {
  int i;
  int j;
  int i_beg;
  int i_end;
  int j_beg;
  int j_end;
  double s_match;
  
  bool force_match = (fseq == 1 || fseq == 2) && a > -1 && b > -1 && os < 0;
  bool force_gap = (fseq == 1 || fseq == 2) && a > -1 && b > -1 && os > -1;

  //cout<<force_match<<","<<force_gap<<","<<fseq<<","<<a<<","<<b<<","<<os<<","<<diff_beg<<","<<diff_end<<endl;

  assert(force_match || force_gap || (fseq < 0 && a < 0 && b < 0 && os < 0));

  //1. Skip computation of scores before update region (diff_beg -> diff_end)
  i_beg = (fseq == 1) ? diff_beg : 1;
  i_end = (fseq == 1) ? diff_end : n;
  j_beg = (fseq == 2) ? diff_beg : 1;
  j_end = (fseq == 2) ? diff_end : m;

  //2. Redo computation of update region
  for(i=i_beg; i<=i_end; i++) {
    for(j=j_beg; j<=j_end; j++) {

      //Initialize row 0 and column 0
      if(i==0 || j==0) {
	E[i][j] = scoring->gapopen() - scoring->gapext();
	    
	if(i==0 && j>0) {
	  E[i][j] = E[i][j-1] + scoring->gapext();
	  S[i][j] = E[i][j];
	}
	
      } else {
	s_match = scoring->match(seq1[i-1], seq2[j-1]);

	if((force_gap //Forced Gap, boundry match
	    && ((fseq == 1 
		 && ((i >= a - os 
		      && i <= a - 1)
		     || (i >= a + b
			 && i <= a + b + os)))
		|| (fseq == 2 
		    && ((j >= a - os 
			 && j <= a - 1)
			|| (j >= a + b 
			    && j <= a + b + os)))))
	   //or Forced Match
	   || (force_match
	       && ((fseq == 1 && i >= a && i <= a + b + 1)
		   || (fseq == 2 && j >= a && j <= a + b + 1)))
	   ) {
		    
	  E[i][j] = scoring->minusInfty(); //S[i][j-1] + scoring->gapopen();
	  F[i][j] = scoring->minusInfty(); //S[i-1][j] + scoring->gapopen();
	  S[i][j] = S[i-1][j-1] + s_match;
		    
	} else if(force_gap && fseq == 1 //Forced Gap Seq 1 (Meaning Seq 2 has the deletion)
		  && i >= a && i <= a + b - 1) {
		    
	  E[i][j] = S[i][j] + scoring->gapopen();
	  F[i][j] = max(S[i-1][j] + scoring->gapopen(),
			F[i-1][j] + scoring->gapext());		    
	  S[i][j] = F[i][j];

	} else if(force_gap && fseq == 2 //Forced Gap Seq 2 (Meaning Seq 1 has the deletion)
		  && j >= a && j <= a + b - 1) {
		    
	  E[i][j] = max(S[i][j-1] + scoring->gapopen(),
			E[i][j-1] + scoring->gapext());
	  F[i][j] = S[i][j] + scoring->gapopen();
	  S[i][j] = E[i][j];

	} else { //No force
	  E[i][j] = max(S[i][j-1] + scoring->gapopen(), E[i][j-1] + scoring->gapext());
	  F[i][j] = max(S[i-1][j] + scoring->gapopen(), F[i-1][j] + scoring->gapext());
	  S[i][j] = max(S[i-1][j-1] + s_match,
			E[i][j],
			F[i][j]);

	}
      }

      //Scores should never be near minusInfty, unless E or F is set to it exactly
      assert(S[i][j] > scoring->minusInfty() + 100 &&
	     (E[i][j] > scoring->minusInfty() + 100 || E[i][j] == scoring->minusInfty()) &&
	     (E[i][j] > scoring->minusInfty() + 100 || E[i][j] == scoring->minusInfty()));
    } //j
  } //i

  //3. Redo scores after update region that have changed; Once row and column are unchanged, stop.
  i_beg = (fseq == 1) ? diff_end+1 : 1;
  j_beg = (fseq == 2) ? diff_end+1 : 1;
  bool i_unChanged = false;
  bool j_unChanged = false;
  bool rowOrCol = true;
  int rc;
  int rc_beg;
  int rc_end;
  double curS;
  double curE;
  double curF;

  while(!(i_unChanged && j_unChanged) && i_beg <= n && j_beg <= m) {
    i_unChanged = true;
    j_unChanged = true;

    rc_beg = (rowOrCol)?i_beg:j_beg;
    rc_end = (rowOrCol)?n:m;

    //rc increments i then j (iterated by rowOrCol) to calculate one row and one column at a time, until both are unchanged
    for(rc=rc_beg; rc<=rc_end; rc++) {
      i = (rowOrCol) ? rc : i_beg;
      j = (rowOrCol) ? j_beg : rc;

      s_match = scoring->match(seq1[i-1], seq2[j-1]);
      
      curE = max(S[i][j-1] + scoring->gapopen(), E[i][j-1] + scoring->gapext());
      curF = max(S[i-1][j] + scoring->gapopen(), F[i-1][j] + scoring->gapext());
      curS = max(S[i-1][j-1] + s_match, curE, curF);
      
      if(curE != E[i][j] || curF != F[i][j] || curS != S[i][j]) {
	if(rowOrCol) i_unChanged = false;
	else j_unChanged = false;
      }

      S[i][j] = curS;
      E[i][j] = curE;
      F[i][j] = curF;

      //Scores should never be near minusInfty, unless E or F is set to it exactly
      assert(S[i][j] > scoring->minusInfty() + 100 &&
	     (E[i][j] > scoring->minusInfty() + 100 || E[i][j] == scoring->minusInfty()) &&
	     (E[i][j] > scoring->minusInfty() + 100 || E[i][j] == scoring->minusInfty()));
      
    }

    //Switch rc to next row or column that is not unchanged on the last iteration
    if(rowOrCol) { 
      j_beg++;
      if(!j_unChanged) 
	rowOrCol = false;
    } else {
      i_beg++;
      if(!i_unChanged) 
	rowOrCol = true;
    }

  } //while not both unchanged
  

  //Find End of Fit Alignment
  maxS = S[n][m]-1;
  maxi = n;
  maxj = m;
  for(i=n; i>=m; i--) {
    if(S[i][m] > maxS) {
      maxS = S[i][m];
      maxi = i;
      maxj = m;
    }
  }
 
  traceback();
 
  return maxS;
}
