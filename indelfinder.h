/*
 PUIalign - Phylogenetically Unambiguous Indel Alignment
 
 Created by John P. McCrow - 9/5/2007
 
 Citation:
 McCrow JP, Alignment of phylogenetically unambiguous indels in Shewanella,
 Journal of Computational Biology 16 (11), 1517-1528 (2009)
 */
struct indel_seq {
  string species;
  string prot_id;
  string prot_seq;
};

struct seqgroup {
  bool ispaired;
  vector< int > seqindex;
  vector< int > seqbegpos;

  seqgroup() { ispaired = false; }
};

struct groupgap {
  int gapgroupid;
  int insgroupid;
  int gapsize;
  int bestgapgroupseqid;
  int bestinsgroupseqid;

  groupgap() { reset();}
  void reset(void);
  bool isgood(void);
};

class indelfinder {
protected:
  int alignradius;
  int alignradiusext;
  int matchradius;
  int min_match_score;
  int min_gap_score;
  int max_match_score_diff_lower;
  int max_match_score_diff_upper;
  int max_gap_score_diff;
  bool findbest;
  int maxparts;
  align_score *scorematrix;
  bool disjoint;
  double total_pairwise_gapscore;
  double total_pairwise_bestscore;
  int total_pairs;
  double maxmdboundl;
  double minmdboundu;
  double maxminmgsd;
  double lowestscore;
  
  vector< indel_seq > seqlist;
  vector< int > indelposlist;
  vector< vector< pairwise_al > > bestaligns;
  vector< vector< pairwise_al > > matchaligns;
  vector< int > seqgroupid;
  vector< int > seqgroupsize; 
  vector< seqgroup > seqgrouplist;
  vector< vector< bool > > inbounds;

  void find_relative_gap_size(groupgap *gg, int g1, int g2);
  void fit_indel(alignedgap *ag, int g_ins, int g_del, int gs_ins, int gs_gap, int gs_ins_ref, int posl, int posr, align_score *s);
  void convertpos(intpair *ip, int groupid, int groupseqid_in, int posl, int posr, int groupseqid_out);
  int convertpos(int groupid, int groupseqid_in, int pos, int groupseqid_out);
  void tagalignment(pairwise_al *al, int seq1id, int seq2id, int seq1offset, int seq2offset);

public:
  indelfinder() { resetall(); }
  void add_indel_seq(indel_seq *s);
  void add_indel_seq(string sp, string pi, string s);
  void add_indel_pos(int pos);
  void set_indel_pos_list(vector< int > *ipl);
  void clear_alignment_vectors(void);
  void resetall(void);
  void resetindel(void);
  bool isdisjoint(void) {return disjoint; }
  void get_groups(vector< seqgroup > *sg) {sg = &seqgrouplist; }
  void set_params(int alignradius, int alignradiusext, int matchradius, align_score *s, int min_match_score, int min_gap_score, int max_match_score_diff_lower, int max_match_score_diff_upper, int max_gap_score_diff, bool findbest, int maxparts);
  void find_groups(void);
  bool isready(void);
  double get_total_pairwise_gapscore(void) {return total_pairwise_gapscore;}
  double get_total_pairwise_bestscore(void) {return total_pairwise_bestscore;}
  int get_total_pairs(void) {return total_pairs;}
  double get_mdlbound(void) {return maxmdboundl;}
  double get_mdubound(void) {return minmdboundu;}
  double get_gdlbound(void) {return maxminmgsd;}
};
