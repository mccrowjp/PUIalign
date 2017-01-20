/*
 PUIalign - Phylogenetically Unambiguous Indel Alignment
 
 Created by John P. McCrow - 9/5/2007
 
 Citation:
    McCrow JP, Alignment of phylogenetically unambiguous indels in Shewanella,
    Journal of Computational Biology 16 (11), 1517-1528 (2009)
*/
#include <string>
#include <vector>
#include <map>
#include "indelfinder.cpp"
#include "progress.h"

void readParams(char *f);
int countindels(char *indelfile);

string _scoringMatrix;
double _gapopen;
double _gapext;
int _alignradius;
int _alignradiusext;
int _matchradius;
int _minmatch;
int _mingap;
int _minmatchdiffl;
int _minmatchdiffu;
int _mingapdiff;
bool _findbestgaps = false;
int _maxparts;

vector<int> split_i(string s, char d);
vector<string> split(string s, char d);

int main(int argc, char *argv[]) {
  align_score s;
  string line; 
  string gene_name;
  int indel_len;
  vector< string > indel_lp;
  vector< int > indel_pos_list;
  vector< string > seq;
  vector< indel_seq > protein_list;
  char *indelfile;
  char *paramsfile;
  indelfinder indfind;
  progressbar *pb;

  if(argc > 2) {
    indelfile = argv[1];
    paramsfile = argv[2];
  } else {
      cout<<"Phylogenetically Unambiguous Indel Alignment (PUIalign) v1.0"<<endl;
      cout<<"Created by John P. McCrow (9/5/2007)"<<endl<<endl;
      cout<<"Usage: "<<argv[0]<<" [Indel List File] [Parameters File]"<<endl<<endl;
      exit(1);
  }

  readParams(paramsfile);

  cout<<"#Indel List File: "<<indelfile<<endl;
  cout<<"#Parameters File: "<<paramsfile<<endl;
  cout<<"#Scoring Matrix: "<<_scoringMatrix<<endl;
  cout<<"#Gap Open Score: "<<_gapopen<<endl;
  cout<<"#Gap Extend Score: "<<_gapext<<endl;
  cout<<"#Pairwise Alignment Radius: "<<_alignradius<<endl;
  cout<<"#Pairwise Alignment Overhang: "<<_alignradiusext<<endl;
  cout<<"#Pairwise Match Radius: "<<_matchradius<<endl;
  cout<<"#Min Match Score: "<<_minmatch<<endl;
  cout<<"#Min Gap Score: "<<_mingap<<endl;
  cout<<"#Max Match Score Difference: "<<_minmatchdiffl<<" - "<<_minmatchdiffu<<endl;
  cout<<"#Max Gap Score Difference: "<<_mingapdiff<<endl;
  cout<<"#Find Best Gaps: "<<(_findbestgaps?"True":"False")<<endl;
  cout<<"#Max Parts: "<<_maxparts<<endl;
  
  s.read_score_matrix(_scoringMatrix);
  s.set_gap_scores(_gapopen,_gapext);

  pb = new progressbar(countindels(indelfile));

  ifstream infile(indelfile);
  int count = 0;
  int totcount = 0;
  bool restart = false;

  while (infile.good() && getline(infile,line,'\n')) {  //Each line    

    if(line.length() < 1) {

    } else if(line.substr(0,1) == ">") {
      indfind.resetall();
      gene_name = line.substr(1,line.length()-1);
      indfind.set_params(_alignradius, _alignradiusext, _matchradius, &s, _minmatch, _mingap, _minmatchdiffl, _minmatchdiffu, _mingapdiff, _findbestgaps, _maxparts);
      count = 0;

    } else if(line.substr(0,1) == "#") {
      indfind.resetindel();
      count++;
      totcount++;
      indel_lp = split(line.substr(1,line.length()-1), ':');
      indel_len = atoi(indel_lp[0].c_str());
      indel_pos_list = split_i(indel_lp[1], ',');
      indfind.set_indel_pos_list(&indel_pos_list);

      assert(indfind.isready());

      restart = true;
      //Place debug lines here for subsets of potential indels
      //if(gene_name == "SV_SO3389.aln" && count == 3) restart=true;
      if(totcount >= 3500) restart=true;
      else restart=false;

      if(restart) {
	cout<<">"<<gene_name<<","<<count<<endl;
	indfind.find_groups();
      }

      pb->update(totcount);
    } else {
      seq = split(line, '\t');
      if(seq.size() == 3 && seq[0].length() > 0 && seq[1].length() > 0) {
	indfind.add_indel_seq(seq[0], seq[1], seq[2]);
      }
    }

  }

  pb->erase();
}

vector<int> split_i(string s, char d) {
  istringstream iss;
  iss.str(s);
  vector<int> list;

  while (iss.good()) {
    string item;
    getline(iss, item, d);
    list.push_back(atoi(item.c_str()));
  }

  return list;
}

vector<string> split(string s, char d) {
  istringstream iss;
  iss.str(s);
  vector<string> list;

  while (iss.good()) {
    string item;
    getline(iss, item, d);
    list.push_back(item);
  }

  return list;
}

void readParams(char *f) {
  string line;
  ifstream infile(f);
  int col;

  while (infile.good() && getline(infile,line,'\n')) {  //Each line
    istringstream iss;
    iss.str(line);
    string v;
    string paramname;
    col = 0;

    while (iss.good() && getline(iss, v, '=')) { //Each param name and value
      if(col == 0) {
        paramname = v;
      } else {
        if(paramname.find("scoringmatrix") != string::npos)
          _scoringMatrix = v;
        else if(paramname.find("gapopen") != string::npos)
          _gapopen = atof(v.c_str());
        else if(paramname.find("gapext") != string::npos)
          _gapext = atof(v.c_str());
        else if(paramname.find("alignradiusext") != string::npos)
            _alignradiusext = atoi(v.c_str());
        else if(paramname.find("alignradius") != string::npos)
          _alignradius = atoi(v.c_str());
        else if(paramname.find("matchradius") != string::npos)
          _matchradius = atoi(v.c_str());
        else if(paramname.find("minmatchscoredifflower") != string::npos)
          _minmatchdiffl = atoi(v.c_str());
        else if(paramname.find("minmatchscorediffupper") != string::npos)
          _minmatchdiffu = atoi(v.c_str());
        else if(paramname.find("mingapscorediff") != string::npos)
          _mingapdiff = atoi(v.c_str());
        else if(paramname.find("minmatchscore") != string::npos)
          _minmatch = atoi(v.c_str());
        else if(paramname.find("mingapscore") != string::npos)
          _mingap = atoi(v.c_str());
        else if(paramname.find("findbestgaps") != string::npos)
          _findbestgaps = (v[0] == 't' || v[0] == 'T');
        else if(paramname.find("maxparts") != string::npos)
          _maxparts = atoi(v.c_str());
        else
          cout<<"#Unknown Parameter "<<paramname<<" = "<<v<<endl;
      }
      col++;
    }
  }
}

int countindels(char *indelfile) {  
  ifstream ifs(indelfile);
  int indelcount = 0;
  string line; 

  while (ifs.good() && getline(ifs,line,'\n')) {  //Each line    
    if(line.length() >= 1 && line.substr(0,1) == "#")
      indelcount++;
  }
  ifs.close();

  return indelcount;
}
