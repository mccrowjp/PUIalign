/*
 PUIalign - Phylogenetically Unambiguous Indel Alignment
 
 Created by John P. McCrow - 9/5/2007
 
 Citation:
 McCrow JP, Alignment of phylogenetically unambiguous indels in Shewanella,
 Journal of Computational Biology 16 (11), 1517-1528 (2009)
 */
using namespace std;

#include <cassert>
#include "seq_align.cpp"
#include "indelfinder.h"

void groupgap::reset(void) {
  gapgroupid=insgroupid=gapsize=bestgapgroupseqid=bestinsgroupseqid=-1;
}

bool groupgap::isgood(void) {
  return (gapgroupid>=0 && 
	  insgroupid>=0 && 
	  gapsize>0 && 
	  bestgapgroupseqid>=0 && 
	  bestinsgroupseqid>=0);
}

void indelfinder::add_indel_seq(indel_seq *s) {
  seqlist.push_back(*s);
}

void indelfinder::add_indel_seq(string sp, string pi, string s) {
  indel_seq is;

  is.species = sp;
  is.prot_id = pi;
  is.prot_seq = s;

  seqlist.push_back(is);
}

void indelfinder::add_indel_pos(int pos) {
  indelposlist.push_back(pos);
}

void indelfinder::set_indel_pos_list(vector< int > *ipl) {
  indelposlist = *ipl;
}

void indelfinder::clear_alignment_vectors(void) {
  bestaligns.clear();
  matchaligns.clear();
  seqgroupid.clear();
  seqgroupsize.clear(); 
  seqgrouplist.clear();
  inbounds.clear();
}

void indelfinder::resetall(void) {
  alignradius=0;
  alignradiusext=0;
  matchradius=0;
  min_match_score=0;
  min_gap_score=0; 
  max_match_score_diff_lower=0;
  max_match_score_diff_upper=0;
  max_gap_score_diff=0; 
  seqlist.clear();

  resetindel();
}

void indelfinder::resetindel(void) {
  disjoint=false; 
  total_pairs=0; 
  total_pairwise_gapscore=0;
  total_pairwise_bestscore=0;
  lowestscore=0;
  indelposlist.clear();

  clear_alignment_vectors();
}

bool indelfinder::isready(void) {
  return (seqlist.size() > 0 &&
	  seqlist.size() == indelposlist.size() &&
	  alignradius > 0 &&
	  alignradiusext > 0 && 
	  matchradius > 0 &&
	  max_match_score_diff_lower >= 0 &&
	  max_match_score_diff_upper >= 0 &&
	  max_gap_score_diff >= 0);
}

void indelfinder::set_params(int alignradius, int alignradiusext, int matchradius, align_score *s, int min_match_score, int min_gap_score, int max_match_score_diff_lower, int max_match_score_diff_upper, int max_gap_score_diff, bool findbest, int maxparts) {
  this->alignradius = alignradius;
  this->alignradiusext = alignradiusext;
  this->matchradius = matchradius;
  this->scorematrix = s;
  this->min_match_score = min_match_score;
  this->min_gap_score = min_gap_score;
  this->max_match_score_diff_lower = max_match_score_diff_lower;
  this->max_match_score_diff_upper = max_match_score_diff_upper;
  this->max_gap_score_diff = max_gap_score_diff;
  this->findbest = findbest;
  this->maxparts = maxparts;

  assert(alignradius >= matchradius);
}

void indelfinder::find_groups(void) {
  int i,j,k,l;
  int sgi, sgj;
  pairwise_al *al;
  vector< vector< alignedgap > > algaps;
  alignedgap *algap;
  groupgap *indelsize;
  int posl, posr;
  int g1,g2;

  //Reset alignment vectors
  clear_alignment_vectors();

  for(i=0; i<seqlist.size(); i++) {
    vector< pairwise_al > ba_row;
    vector< pairwise_al > ma_row;
    vector< bool > ib_row;
    for(j=0; j<seqlist.size(); j++) {
      if(j==i) {
	pairwise_al alb;
	ba_row.push_back(alb);

	pairwise_al alm;
	ma_row.push_back(alm);

	bool ib;
	ib_row.push_back(ib);
      } else {

	//Allow alignradius to vary if it runs off the ends of the sequence
	int maxoverhangl = max(0,
			      max((int)(0 - (indelposlist[i] - alignradius - alignradiusext)),
				  (int)(0 - (indelposlist[j] - alignradius))));
	int maxoverhangr = max(0,
			       max((int)(indelposlist[i] + alignradius + alignradiusext - seqlist[i].prot_seq.length()),
				   (int)(indelposlist[j] + alignradius + alignradiusext - seqlist[j].prot_seq.length())));

	int alignradiusl = alignradius-maxoverhangl;
	int alignradiusr = alignradius-maxoverhangr;

	int s1begpos = indelposlist[i] - alignradiusl - alignradiusext;
	int s2begpos = indelposlist[j] - alignradiusl;
	int s1len = alignradiusl + alignradiusr + (alignradiusext*2);
	int s2len = alignradiusl + alignradiusr;

	if(alignradiusl >= matchradius &&
	   alignradiusr >= matchradius &&
	   s1begpos >= 0 &&
	   s2begpos >= 0 &&
	   s1begpos + s1len <= seqlist[i].prot_seq.length() &&
	   s2begpos + s2len <= seqlist[j].prot_seq.length()
	   ) {

	  ib_row.push_back(true);

	  //Set substrings for fit alignment
	  string s1 = seqlist[i].prot_seq.substr(s1begpos, s1len);
	  string s2 = seqlist[j].prot_seq.substr(s2begpos, s2len);

	  //Perform fit alignment of s2 fit to larger s1
	  fit_align fa(s1, s2, scorematrix);
	  fa.fit();
	  al = fa.get_alignment();
	  tagalignment(al, i, j, s1begpos, s2begpos);
	  ba_row.push_back(*al);

	  //al->print();

	  //If best alignment has a gap, realign with forced match before adding alignment to matchaligns
	  double matchscore;
	  if(al->hasgap()) {
	    matchscore = fa.refit_match(2,
					(indelposlist[j]-s2begpos-matchradius >= 0 ? indelposlist[j]-s2begpos-matchradius : 0),
					(matchradius*2),
					(indelposlist[j]-s2begpos-matchradius >= 0 ? indelposlist[j]-s2begpos-matchradius : 0),
					(indelposlist[j]-s2begpos-matchradius >= 0 ? indelposlist[j]-s2begpos-matchradius : 0)+(matchradius*2));

	    /*
	    //fa.get_alignment()->print();

	    double matchscore2 = fa.fit_match(2,
	    (indelposlist[j]-s2begpos-matchradius >= 0 ? indelposlist[j]-s2begpos-matchradius : 0),
	    (matchradius*2));

	    //fa.get_alignment()->print();
	    //cout<<"Refit:"<<matchscore<<" Fit:"<<matchscore2<<endl;
	  
	    assert(matchscore == matchscore2);
	    */
	  }

	  al = fa.get_alignment();
	  tagalignment(al, i, j, s1begpos, s2begpos);
	  ma_row.push_back(*al);
	} else {
	  //Pair i,j is out of bounds in this orientation
	  ib_row.push_back(false);
	  al = new(pairwise_al);
	  tagalignment(al, i, j, s1begpos, s2begpos);
	  ba_row.push_back(*al);
	  al = new(pairwise_al);
	  tagalignment(al, i, j, s1begpos, s2begpos);
	  ma_row.push_back(*al);
	}
      }      
    }
    bestaligns.push_back(ba_row);
    matchaligns.push_back(ma_row);
    inbounds.push_back(ib_row);
  }

  //Check all pairs for at least one orientation in bounds
  for(i=0; i<seqlist.size()-1; i++) {
    for(j=i+1; j<seqlist.size(); j++) {
      if(inbounds[i][j] || inbounds[j][i]) {
      }	else {
	cout<<"#Out of bounds"<<endl;
	return;
      }
    }
  }

  //Initialize all sequences in their own group
  for(i=0; i<seqlist.size(); i++) {
    seqgroupid.push_back(i);
    seqgroupsize.push_back(1);
  }

  int max_match_score_diff;
  disjoint = false;
  for(max_match_score_diff=max_match_score_diff_lower; 
      !disjoint && max_match_score_diff <= max_match_score_diff_upper; 
      max_match_score_diff++) {

    //Reset partition
    for(i=0; i<seqlist.size(); i++) {
      seqgroupid[i] = i;
      seqgroupsize[i] = 1;
    }

    //cout<<"Maxdiff: "<<max_match_score_diff<<endl;

    //Find all maximal disjoint cliques, by merging adjacent groups
    disjoint = true;
    for(i=0; i<seqlist.size()-1; i++) {
      for(j=i+1; j<seqlist.size(); j++) {
	if(seqgroupid[i] != seqgroupid[j] &&
	   ((inbounds[i][j] && bestaligns[i][j].score-matchaligns[i][j].score <= max_match_score_diff) ||
	    (inbounds[j][i] && bestaligns[j][i].score-matchaligns[j][i].score <= max_match_score_diff))) {

	  //cout<<i<<","<<j<<":"<<bestaligns[i][j].score<<"-"<<matchaligns[i][j].score<<"|"<<bestaligns[j][i].score<<"-"<<matchaligns[j][i].score<<endl;
	  //bestaligns[i][j].print();
	  //matchaligns[i][j].print();

	  //Check disjoint clique property, and merge groups
	  sgi = seqgroupid[i];
	  sgj = seqgroupid[j];
	  for(k=0; k<seqlist.size(); k++) {
	    if(seqgroupid[k] == sgi) {
	      for(l=0; l<seqlist.size(); l++) {
		if(seqgroupid[l] == sgj) {

		  //cout<<i<<","<<j<<","<<k<<","<<l<<","<<bestaligns[i][j].score-matchaligns[i][j].score<<","<<bestaligns[j][i].score-matchaligns[j][i].score<<","<<bestaligns[k][l].score-matchaligns[k][l].score<<","<<bestaligns[l][k].score-matchaligns[l][k].score<<","<<sgi<<endl;

		  if((inbounds[k][l] && bestaligns[k][l].score-matchaligns[k][l].score <= max_match_score_diff) ||
		     (inbounds[l][k] && bestaligns[l][k].score-matchaligns[l][k].score <= max_match_score_diff)) {
		  } else {
		    //cout<<i<<","<<j<<","<<k<<","<<l<<","<<seqgroupid[k]<<","<<seqgroupid[l]<<":"<<bestaligns[k][l].score-matchaligns[k][l].score<<","<<bestaligns[l][k].score-matchaligns[l][k].score<<endl;
		    //bestaligns[k][l].print();
		    //matchaligns[k][l].print();
		  
		    disjoint = false;
		    i = j = l = k = seqlist.size(); //Break out of all loops for this max_match_score_diff
		  }
		}
	      }
	    }
	  }

	  //Merge groups j->i
	  if(disjoint) {
	    for(l=0; l<seqlist.size(); l++) {
	      if(seqgroupid[l] == sgj) {
		seqgroupid[l] = sgi;
		seqgroupsize[sgi]++;
		seqgroupsize[sgj]--;
	      }
	    }
	  }

	}
      }
    }
    //cout<<max_match_score_diff<<","<<disjoint<<endl;
  }

  //If disjoint, then continue with between group gap alignments
  if(disjoint) {
    /*    
    for(i=0; i<seqlist.size(); i++) {
      cout<<seqgroupsize[i]<<",";
    }
    cout<<endl;
    for(i=0; i<seqlist.size(); i++) {
      cout<<seqgroupid[i]<<",";
    }
    cout<<endl;
    */

    //Set groupings with size > 1
    for(i=0; i<seqgroupsize.size(); i++) {
      if(seqgroupsize[i] > 1) {
	//cout<<"Group "<<i<<": ";
	seqgroup sg;
	for(j=0; j<seqgroupid.size(); j++) {
	  if(seqgroupid[j] == i) {
	    //cout<<j<<",";
	    if(sg.seqbegpos.empty()) {
	      sg.seqbegpos.push_back(0);
	    } else {
	      if(sg.seqbegpos.size() == 1) {
		sg.seqbegpos[0] = (matchaligns[sg.seqindex[0]][j].aligned[0].i + 
				   matchaligns[sg.seqindex[0]][j].seq1_offset - 1);
		//cout<<seqlist[sg.seqindex[0]].prot_seq.substr(sg.seqbegpos[0],alignradius*2)<<endl;
	      }
	      sg.seqbegpos.push_back(sg.seqbegpos[0] - 
				     matchaligns[sg.seqindex[0]][j].seq1_offset -
				     matchaligns[sg.seqindex[0]][j].aligned[0].i +
				     matchaligns[sg.seqindex[0]][j].aligned[0].j +
				     matchaligns[sg.seqindex[0]][j].seq2_offset);

	      //cout<<seqlist[j].prot_seq.substr(sg.seqbegpos[sg.seqbegpos.size()-1],alignradius*2)<<endl;
	    }

	    sg.seqindex.push_back(j);
	  }
	}
	seqgrouplist.push_back(sg);
	//cout<<endl;
      }
    }

    seqgroupsize.clear();

    //Update new seq group ids
    for(i=0; i<seqgrouplist.size(); i++) {
      for(j=0; j<seqgrouplist[i].seqindex.size(); j++) {
	seqgroupid[seqgrouplist[i].seqindex[j]]=i;
      }
    }

    if(seqgrouplist.size() > 1 && seqgrouplist.size() <= maxparts) {
      
      //Find match_diff parameter bounds
      bool mdboundlundef = true;
      bool mdbounduundef = true;
      bool lowestscoreundef=true;  
      double sc1;
      double sc2;
      double minscdiff;
      double maxbsc;
      for(i=0; i<seqgrouplist.size(); i++) {
	for(k=0; k<seqgrouplist[i].seqindex.size(); k++) {
	  for(j=0; j<seqgrouplist.size(); j++) {
	    for(l=0; l<seqgrouplist[j].seqindex.size(); l++) {	      
	      if(i!=j || k!=l) {
		sc1 = bestaligns[seqgrouplist[i].seqindex[k]][seqgrouplist[j].seqindex[l]].score
		  - matchaligns[seqgrouplist[i].seqindex[k]][seqgrouplist[j].seqindex[l]].score;
		sc2 = bestaligns[seqgrouplist[j].seqindex[l]][seqgrouplist[i].seqindex[k]].score
		  - matchaligns[seqgrouplist[j].seqindex[l]][seqgrouplist[i].seqindex[k]].score;

		if(!inbounds[seqgrouplist[i].seqindex[k]][seqgrouplist[j].seqindex[l]]) {
		  minscdiff=sc2;
		  maxbsc=bestaligns[seqgrouplist[j].seqindex[l]][seqgrouplist[i].seqindex[k]].score;
		} else if(!inbounds[seqgrouplist[j].seqindex[l]][seqgrouplist[i].seqindex[k]]) {
		  minscdiff=sc1;
		  maxbsc=bestaligns[seqgrouplist[i].seqindex[k]][seqgrouplist[j].seqindex[l]].score;
		} else {
		  minscdiff = min(sc1,sc2);
		  maxbsc = max(bestaligns[seqgrouplist[j].seqindex[l]][seqgrouplist[i].seqindex[k]].score,
			       bestaligns[seqgrouplist[i].seqindex[k]][seqgrouplist[j].seqindex[l]].score);
		}

		//cout<<i<<","<<k<<","<<j<<","<<l<<","<<seqgrouplist[i].seqindex[k]<<","<<seqgrouplist[j].seqindex[l]<<","<<sc1<<"-"<<sc2<<","<<maxmdboundl<<","<<minmdboundu<<endl;
		
		if(j==i) { //Comparing sequences in the same group
		  if(mdboundlundef || minscdiff > maxmdboundl) {		  
		    maxmdboundl = minscdiff;	
		    mdboundlundef = false;
		  }	
		} else { //Different groups
		  if(mdbounduundef || minscdiff < minmdboundu) {
		    minmdboundu = minscdiff;
		    mdbounduundef = false;
		  }
		}
		if(lowestscoreundef || maxbsc < lowestscore) {
		  lowestscore = maxbsc;
		  lowestscoreundef=false;
		}

	      }
	    }
	  }
	}
      }

      //Create group pair aligned gap matrix
      for(i=0; i<seqgrouplist.size(); i++) {
	vector< alignedgap > agv;
	for(j=0; j<seqgrouplist.size(); j++) {
	  alignedgap agi;
	  agv.push_back(agi);
	}
	algaps.push_back(agv);
      }

      algap = new alignedgap;
      indelsize = new groupgap;

      //cout<<"Groups: "<<seqgrouplist.size()<<endl;
      
      bool maxminmgsdundef = true;

      //Compare each group, find a common indel
      for(g1=0; g1<seqgrouplist.size()-1; g1++) {
	for(g2=g1+1; g2<seqgrouplist.size(); g2++) {

	  //cout<<"Comparing: "<<g1<<"-"<<g2<<endl;

	  indelsize->reset();
	  assert(!indelsize->isgood());

	  find_relative_gap_size(indelsize, g1, g2);
      
	  //cout<<g1<<","<<g2<<","<<indelsize->gapgroupid<<","<<indelsize->insgroupid<<","<<indelsize->gapsize<<","<<indelsize->bestgapgroupseqid<<","<<indelsize->bestinsgroupseqid<<endl;

	  int bestgoodpairs = 0;
	  double besttgs = 0;
	  double besttbs = 0;
	  bool foundgood = false;
	  bool foundgoodpair = false;
	  bool stopl = false;
	  bool stopr = false;
	  bool minmgsdundef = true;
	  double minmgsd;
 
	  if(indelsize->isgood()) {
	    for(i=0; i<=matchradius 
		  && !(stopl && stopr) 
		  && (findbest || !foundgoodpair); i++) { //i is the radius around indelpos to try the indel location posl to posr
	      for(j=-1; j==-1 || (j==1 && i>0); j+=2) { //Alternate between -i and +i around indelpos
		if(!((j==-1 && stopl) || (j==1 && stopr))) {

		  posl = int(indelposlist[seqgrouplist[indelsize->insgroupid].seqindex[indelsize->bestinsgroupseqid]] 
			     + (i*j)
			     - (indelsize->gapsize/2));
		  posr = posl + indelsize->gapsize - 1;
		
		  if(posl >= 0 && 
		     posr >= 0 && 
		     posr < seqlist[seqgrouplist[indelsize->insgroupid].seqindex[indelsize->bestinsgroupseqid]].prot_seq.length()) {
		  
		    int goodpairs = 0;
		    double tgs = 0;
		    double tbs = 0;
		    bool mgsdundef = true;
		    double mgsd;
		    //Try each pairwise deletion between groups
		    for(k=0; k<seqgrouplist[indelsize->insgroupid].seqindex.size(); k++) {
		      for(l=0; l<seqgrouplist[indelsize->gapgroupid].seqindex.size(); l++) {
	  
			algap->reset();
			assert(!algap->isgood());

			fit_indel(algap, indelsize->insgroupid, indelsize->gapgroupid, k, l, indelsize->bestinsgroupseqid, posl, posr, scorematrix);

			//Indel pairwise score without gap penalty must be at least min_gap_score
			//And the difference from the best pairwise score must be at most max_gap_score_diff
			if(algap->isgood() && 
			   algap->gapscore - scorematrix->gapopen() - (scorematrix->gapext() * indelsize->gapsize) >= min_gap_score) {

			  if(mgsdundef || algap->bestscore - algap->gapscore > mgsd) {
			    mgsd = algap->bestscore - algap->gapscore;
			    mgsdundef = false;
			  }

			  if(algap->bestscore - algap->gapscore <= max_gap_score_diff) {
			    goodpairs++;
			    tgs += algap->gapscore;
			    tbs += algap->bestscore;
			    
			    //cout<<algap->gapscore<<","<<tgs<<"/"<<algap->bestscore<<","<<tbs<<":"<<scorematrix->gapopen()<<","<<scorematrix->gapext()<<"*"<<indelsize->gapsize<<">="<<min_gap_score<<endl;
			  }
			  
			  //cout<<"Gap "<<g1<<","<<g2<<","<<k<<","<<l<<"="<<(algap->insposr - algap->insposl + 1);

			  if((!algaps[g1][g2].isgood()) || 
			     (tgs > algaps[g1][g2].gapscore && 
			      (algap->bestscore - algap->gapscore <= max_gap_score_diff || 
			       algaps[g1][g2].bestscore - algaps[g1][g2].gapscore > max_gap_score_diff))) {
			    algaps[g1][g2].copy(algap);
			    //cout<<"*";
			  }
			  //cout<<endl;

			}

		      }
		    }
		    
		    if(!mgsdundef && (minmgsdundef || mgsd < minmgsd)) {
		      minmgsd = mgsd;
		      minmgsdundef = false;
		    }

		    //cout<<g1<<","<<g2<<":"<<posl<<"-"<<posr<<"="<<goodpairs<<","<<tgs<<endl;
	  
		    if(!foundgood ||
		       goodpairs > bestgoodpairs || 
		       (goodpairs == bestgoodpairs && tgs > besttgs) ||
		       (goodpairs == bestgoodpairs && tgs == besttgs && tbs-tgs < besttbs-besttgs)) {
		      bestgoodpairs = goodpairs;
		      besttgs = tgs;
		      besttbs = tbs;		
		      foundgood = true;
		    }
		  		
		    if(bestgoodpairs == seqgrouplist[g1].seqindex.size() * seqgrouplist[g2].seqindex.size()) {
		      if(goodpairs < bestgoodpairs-1 && j == -1) stopl = true;
		      if(goodpairs < bestgoodpairs-1 && j == 1) stopr = true;
		    }
		    
		    if(bestgoodpairs > 0 && bestgoodpairs >= (seqgrouplist[g1].seqindex.size() * seqgrouplist[g2].seqindex.size())) foundgoodpair = true;
		    
		  }
		}
	      } //j
	    } //i

	    /*
	    if(bestgoodpairs > 0 && bestgoodpairs >= (seqgrouplist[g1].seqindex.size() * seqgrouplist[g2].seqindex.size())) {
	      seqgrouplist[g1].ispaired = true;
	      seqgrouplist[g2].ispaired = true;
	      total_pairwise_gapscore += besttgs;
	      total_pairwise_bestscore += besttbs;
	      total_pairs += (seqgrouplist[g1].seqindex.size() * seqgrouplist[g2].seqindex.size());
	    }
	    */

	    //If we display maxmingsd then we can display all partitions
	    seqgrouplist[g1].ispaired = true;
	    seqgrouplist[g2].ispaired = true;
	    total_pairwise_gapscore += besttgs;
	    total_pairwise_bestscore += besttbs;
	    total_pairs += (seqgrouplist[g1].seqindex.size() * seqgrouplist[g2].seqindex.size());
	  
	  } //indelsize isgood

	  //delete indelsize;

	  if(!minmgsdundef && (maxminmgsdundef || minmgsd < maxminmgsd)) {
	    maxminmgsd = minmgsd;
	    maxminmgsdundef = false;
	  }

	} //g2
      } //g1

      delete indelsize;
      delete algap;

      if(maxminmgsdundef) {
	cout<<"#Unable To Find Common Gap"<<endl;
      } else {

	if(total_pairs > 0) {	
	  int grpgaplen;
	  int grpgaplensign;

	  //Print indel partition match and gap diff bounds
	  cout<<maxmdboundl<<","<<minmdboundu<<","<<maxminmgsd<<","<<lowestscore<<endl;
	
	  //Print good groups
	  for(g1=0; g1<seqgrouplist.size(); g1++) {
	    if(seqgrouplist[g1].ispaired) {

	      //Print relative gap lengths between each good group
	      j=0;
	      cout<<"[";
	      for(g2=0; g2<seqgrouplist.size(); g2++) {
		if(g1==g2) {
		  if(j>0) cout<<",";
		  cout<<"0";
		  j++;
		} else if(seqgrouplist[g2].ispaired) {
		  if(j>0) cout<<",";
		  j++;
		  if((g1<g2 && algaps[g1][g2].isgood()) ||
		     (g1>g2 && algaps[g2][g1].isgood())) {
		    grpgaplen = (g1<g2 
				 ? algaps[g1][g2].insposr - algaps[g1][g2].insposl + 1
				 : algaps[g2][g1].insposr - algaps[g2][g1].insposl + 1);
		    assert(grpgaplen >= 0);
		    grpgaplensign = ((g1<g2 && seqgroupid[algaps[g1][g2].insseqid] == g1) || 
				     (g1>g2 && seqgroupid[algaps[g2][g1].insseqid] == g1)
				     ) ? 1 : -1;
		    cout<<(grpgaplen*grpgaplensign);		  
		  } else {
		    cout<<"?";
		  }
		}
	      }
	      cout<<"]:";

	      //Print member protein ids within each group
	      for(j=0; j<seqgrouplist[g1].seqindex.size(); j++) {
		if(j>0) cout<<",";
		cout<<seqlist[seqgrouplist[g1].seqindex[j]].prot_id;

		/*if(minseqtss[seqgrouplist[g1].seqindex[j]] <= maxseqtss[seqgrouplist[g1].seqindex[j]]) 
		  cout<<"("<<minseqtss[seqgrouplist[g1].seqindex[j]]<<":"<<maxseqtss[seqgrouplist[g1].seqindex[j]]<<")";
		  else cout<<"(:)";*/
	      }
	      cout<<endl;
	    }
	  }

	}
      }

    } else {
      if(seqgrouplist.size() < 1) {
	cout<<"#Only Singleton Parts"<<endl;
      } else if(seqgrouplist.size() == 1) {
	cout<<"#Only One Multisequence Part"<<endl;
      } else if(seqgrouplist.size() > maxparts) {
	cout<<"#Too Many Parts: "<<seqgrouplist.size()<<endl;
      }
    }    
  } else {
    cout<<"#No Partition"<<endl;
  }
}

void indelfinder::find_relative_gap_size(groupgap *gg, int g1, int g2) {
  int i,j,k;
  int orient;
  int gap1,gap2;
  int gx,gy,sx,sy,x,y;
  double bestscore;
  bool foundgap1;
  bool setbestscore;
  int lastsize;

  foundgap1 = false;

  for(i=0; i<seqgrouplist[g1].seqindex.size(); i++) {
    for(j=0; j<seqgrouplist[g2].seqindex.size(); j++) {
      for(orient=0; orient<=1; orient++) {

	setbestscore = false;

	if(orient == 0 && inbounds[seqgrouplist[g1].seqindex[i]][seqgrouplist[g2].seqindex[j]]) {
	  gx = g1;
	  gy = g2;
	  x = i;
	  y = j;
	  if(!foundgap1 || bestaligns[seqgrouplist[g1].seqindex[i]][seqgrouplist[g2].seqindex[j]].score > bestscore) {
	    setbestscore = true;
	  }
	} else if(orient == 1 && inbounds[seqgrouplist[g2].seqindex[j]][seqgrouplist[g1].seqindex[i]]) {
	  gx = g2;
	  gy = g1;
	  x = j;
	  y = i;
	  if(!foundgap1 || bestaligns[seqgrouplist[g2].seqindex[j]][seqgrouplist[g1].seqindex[i]].score > bestscore) {
	    setbestscore = true;	    
	  }
	}

	if(setbestscore) {
	  sx = seqgrouplist[gx].seqindex[x];
	  sy = seqgrouplist[gy].seqindex[y];
	  gap1=0;
	  gap2=0;
	  //Count relative gap length in region of indel locus +/- matchradius
	  for(k=1; k<bestaligns[sx][sy].aligned.size(); k++) {
	    if(bestaligns[sx][sy].aligned[k].j >= indelposlist[seqgrouplist[gy].seqindex[y]] - bestaligns[sx][sy].seq2_offset - matchradius &&
	       bestaligns[sx][sy].aligned[k].j <= indelposlist[seqgrouplist[gy].seqindex[y]] - bestaligns[sx][sy].seq2_offset + matchradius) {
	      
	      if(bestaligns[sx][sy].aligned[k-1].i == bestaligns[sx][sy].aligned[k].i) gap1++;
	      if(bestaligns[sx][sy].aligned[k-1].j == bestaligns[sx][sy].aligned[k].j) gap2++;
	    }
	  }

	  //Gap size is in relation to gap1
	  if(gap1 > gap2) {
	    bestscore = bestaligns[sx][sy].score;
	    gg->gapgroupid = gx;
	    gg->insgroupid = gy;
	    gg->gapsize = gap1-gap2;
	    gg->bestgapgroupseqid = x;
	    gg->bestinsgroupseqid = y;
	    
	    assert(gg->isgood());

	    foundgap1 = true;	  

	    //cout<<"Locus: "<<k<<endl;
	    //cout<<sx<<","<<sy<<":"<<orient<<","<<gap1<<","<<gap2<<","<<(gap1-gap2)<<endl;
	    //bestaligns[sx][sy].print();
	  }
	}

      }
    }
  }

}

void indelfinder::fit_indel(alignedgap *ag, int g_ins, int g_del, int gs_ins, int gs_del, int gs_ins_ref, int posl, int posr, align_score *s) {
  double gapscore;
  double bestscore;
  intpair *pos;
  int del_len;
  pairwise_al *al;

  int s1id = seqgrouplist[g_del].seqindex[gs_del];
  int s2id = seqgrouplist[g_ins].seqindex[gs_ins];

  assert(posr >= posl);
  assert(posl >= 0 && posr >= 0);
  assert(gs_ins >= 0 && gs_del >= 0 && gs_ins_ref >= 0);
  assert(g_ins >= 0 && g_del >= 0);
  assert(g_ins != g_del);
 
  //for(i=1; i<seqgrouplist[g_ins].seqindex.size(); i++) {
  //  matchaligns[seqgrouplist[g_ins].seqindex[i]][seqgrouplist[g_ins].seqindex[0]].print();
  //}

  pos = new intpair;
  convertpos(pos, g_ins, gs_ins_ref, posl, posr, gs_ins);

  del_len = pos->j - pos->i + 1;

  if(del_len > 0 && pos->i >= 0 && pos->j >= 0) {

    assert(posr - posl + 1 == del_len);

    //Do not attempt alignment if matchradius extends out of sequences
    if(indelposlist[s1id] - matchradius >= 0 && pos->i - matchradius >= 0) {
      
      int s1begpos = indelposlist[s1id] - alignradius - alignradiusext;
      int s2begpos = pos->i - alignradius;
      int s1len = (alignradius+alignradiusext)*2;
      int s2len = del_len + (alignradius*2);
      
      //Allow alignradius to vary though if it runs off the ends of the sequence
      if(s1begpos + s1len > seqlist[s1id].prot_seq.length()) s1begpos = seqlist[s1id].prot_seq.length() - s1len;
      if(s2begpos + s2len > seqlist[s2id].prot_seq.length()) s2begpos = seqlist[s2id].prot_seq.length() - s2len;
      if(s1begpos < 0) s1begpos = 0;
      if(s2begpos < 0) s2begpos = 0;
      if(s1begpos + s1len > seqlist[s1id].prot_seq.length()) s1len = seqlist[s1id].prot_seq.length() - s1begpos;
      if(s2begpos + s2len > seqlist[s2id].prot_seq.length()) s2len = seqlist[s2id].prot_seq.length() - s2begpos;

      //Do not attempt alignment if matchradius extends out of sequences
      if(s1len >= (alignradius*2) && s2len >= del_len + (alignradius*2)) {

	string s1 = seqlist[s1id].prot_seq.substr(s1begpos, s1len);
	string s2 = seqlist[s2id].prot_seq.substr(s2begpos, s2len);
	
	fit_align fa(s1,s2,s);
	
	bestscore = fa.fit();
	//fa.print();
	gapscore = fa.refit_gap(2,alignradius,del_len,matchradius,alignradius-matchradius,alignradius+del_len+matchradius);
	//assert(gapscore == fa.fit_gap(2,alignradius,del_len,matchradius));
	//fa.print();

	//cout<<s1id<<","<<s2id<<":"<<s1begpos<<","<<s1len<<","<<s2begpos<<","<<s2len<<"="<<bestscore<<"-"<<gapscore<<endl;

	ag->gapseqid = s1id;
	ag->insseqid = s2id;
	ag->insposl = pos->i;
	ag->insposr = pos->j;
	ag->gapscore = gapscore;
	ag->bestscore = bestscore;

	al = fa.get_alignment();
	tagalignment(al, s1id, s2id, s1begpos, s2begpos);	
	ag->gappos = al->findgappos(2, pos->i, pos->j);
	
	//if(!ag->isgood()) cout<<"GP:"<<ag->gappos<<endl;
	//assert(ag->isgood());
      }
    }
  }

  delete pos;

}

void indelfinder::convertpos(intpair *ip, int groupid, int groupseqid_in, int posl, int posr, int groupseqid_out) {
  int offset;
  
  offset = seqgrouplist[groupid].seqbegpos[groupseqid_out] - seqgrouplist[groupid].seqbegpos[groupseqid_in];

  ip->i = posl + offset;
  ip->j = posr + offset;

}

int indelfinder::convertpos(int groupid, int groupseqid_in, int pos, int groupseqid_out) {
  int offset;
  
  offset = seqgrouplist[groupid].seqbegpos[groupseqid_out] - seqgrouplist[groupid].seqbegpos[groupseqid_in];

  return pos + offset;
}

void indelfinder::tagalignment(pairwise_al *al, int seq1id, int seq2id, int seq1offset, int seq2offset) {
  al->seq1_id = seq1id;
  al->seq2_id = seq2id;
  al->seq1_offset = seq1offset;
  al->seq2_offset = seq2offset;
  al->firstgap.insseqid = (al->firstgap.insseqid == 1 ? seq1id : (al->firstgap.insseqid == 2 ? seq2id : al->firstgap.insseqid));
}
