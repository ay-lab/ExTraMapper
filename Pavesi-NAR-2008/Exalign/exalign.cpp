//I added the below two includes and it worked -- Ferhat
#include <cstring>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <list>
#include <math.h>
#include <map>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>

//#define DEBUG1
//#define DEBUG2
//#define DEBUG4
//#define DEBUG5
//#define DEBUG6
//#define DEBUG7
//#define DEBUG8
//#define DEBUG9
//#define DEBUG10
//#define CFREQ
//#define HACK1
#define CIDER

using namespace std;

typedef float m_point;

const char SEPARATORI[1] = {','}; //Aggiungere qui ulteriori separatori di campo del file da parsare, se necessario.

const short int SEPSIZE = 1; //E ricordarsi di aggiornare questo numero al nuovo numero di separatori!!!

const float BIN_STEP = 0.00125;

int EXLINE = 10, FILTER = 0, S_SKIP = 1, E_SKIP = 1, BEST = 5, TOL = 0, N20, SIM = 12, V_BEST = 5, MAX_CYCLE = 3, MAX_EX, MAX_MERGE = 250, RES_COUNT = 0;

int NWS = 0, PWS = 0;

float OGAP = -2, CGAP = -1, EXACT_MATCH_BONUS = 0;

bool WEB = false, M3 = true, local = false, semilocal = false, freq = false, query = true, table = false, gensel = false, querysel = false,tabapp = false, outapp = false, score_out = false, at_least_one = false, last_one = false, dynamic_gap = true, custom = false, blasthitonly = false, useall = false;

const int INF = 2000000000, MAX_DIFF = 21;

vector<int> TOT_EXONS, TOT_GENES;

const string Hs = "Homo_sapiens.rf", Mm = "Mus_musculus.rf", Rn = "Rattus_norvegicus.rf", Dr = "Danio_rerio.rf", Fr = "Fugu_rubripes.rf", 
	     Dm = "Drosophila_melanogaster.rf",Gg = "Gallus_gallus.rf", Xt = "Xenopus_tropicalis.rf", Ci = "Ciona_intestinalis.rf",
	     Ce = "Caenorhabditis_elegans.rf", Pt = "Pan_troglodites.rf", Bt = "Bos_taurus.rf", Md = "Monodelphis_domestica.rf",
             Tn = "Tetraodon_nigroviridis.rf", Am = "Apis_mellifera.rf", Sp = "Strongylocentrotus_purpuratus.rf";

const string HsA = "Homo_sapiens.all", MmA = "Mus_musculus.all", RnA = "Rattus_norvegicus.all", DrA = "Danio_rerio.all", FrA = "Fugu_rubripes.all",                               DmA = "Drosophila_melanogaster.all",GgA = "Gallus_gallus.all", XtA = "Xenopus_tropicalis.all", CiA = "Ciona_intestinalis.all",                                     CeA = "Caenorhabditis_elegans.all", PtA = "Pan_troglodites.all", BtA = "Bos_taurus.all", MdA = "Monodelphis_domestica.all",                                        TnA = "Tetraodon_nigroviridis.all", AmA = "Apis_mellifera.all", SpA = "Strongylocentrotus_purpuratus.all";

#ifndef CIDER
const string HTML_DIR = "/home/federico/websrv/htdocs/ex_html/";
#else
const string HTML_DIR = "/Library/WebServer/Documents/exalign/ex_html/";
#endif

string queryfile, genesfile, outfile, freqfile, tabfile, gsel, qsel, QID, wtmpfile, tmps, actual_query_n, actual_query_p;
vector<string> genesfiles,freqfiles; //new

ostringstream wfile;

int x,y;
enum{X,Y,MX};

struct gene 
{
	bool merged;
	bool tri_merged;
	vector<string> merged_with;
	string name;
	string name2;
	string protein_acc;
	string chr;
	string family;
	char strand;
	int cdsS, cdsE, cdsD[4];
	int exN;
	vector<int> exL;
	vector<int> ex_or_L;
	vector<short int> exF;
	vector<int> merg_pos;
//	vector<bool> m_flag;
	m_point or_score;
	int filenum;
	string fileorg;
};

struct match 
{
	int pos[2];
	int exN[2];
};

struct Frame_Freq
{
	int filenum;
	vector<m_point> Framefreq;
};

vector<gene> parser(const char*,int,int,int);
string linecleaner(string);
vector<int> exL_parser(string,string, int, char, int, int, vector<short int>*, int*);
//vector<short int> exF_parser(string, int, char); 
vector<match> sizecomp(vector<gene>::iterator, vector<gene>::iterator, vector<gene>::iterator, vector<gene>::iterator);
//vector<gene> fus_generator(vector<gene>);
vector<int> match_merge(int, vector<match>::iterator,vector<match>::iterator,vector<match>::iterator,vector<match>::iterator,vector<match>::iterator,vector<match>::iterator);
void main_aligner(vector<gene>::iterator, vector<gene>::iterator, gene, vector<int>::iterator, vector<int>::iterator,
                  vector< map<int,double> >::iterator, vector< map<int,double> >::iterator, vector<Frame_Freq>*); 
m_point aligner(gene*,gene*,bool,double*,long double*,double,m_point*,bool);
void m_init(m_point***,m_point****); 
void m_filler(m_point***,m_point****,vector<vector <int> >,double*,double,vector<vector <short int> >, m_point*);
vector<int> remFL(vector<int>);
vector<short int> remFL(vector<short int>);
//m_point c_check(int,int, double*, double);
m_point c_check(int,int,double*,short int,short int,m_point*,int,double max_diff);
m_point freq_c_check(int,int,double*,short int,short int,m_point*,int,double max_diff);
m_point max(vector<m_point>,int*);
void m_pmax(m_point***,int,int,int*);
vector<string> m_traceback(m_point***, m_point****,vector<int>,vector<int>,int*,vector<int>,vector<int>,long double*,int,vector<short int>*,vector<short int>*, m_point*);
void m_crawler(m_point***, m_point*,int*);
void s_out(m_point**, vector<string>, int*, gene, gene, m_point);
vector<string> s_scan(string);
vector<int> best_pos(vector<m_point>,int);
void clineparser(int, char**);
void displayhelp();
void argcheck();
m_point lmax(vector<m_point>,int*);
m_point m_max(m_point***,int*);
int gap_count(vector<string>);
//vector<gene> l_merged(vector<gene>);
void freqgen(vector<gene>::iterator,vector<gene>::iterator,int);
map<int,double> getfreq(m_point*, bool);
vector<gene> fus_match(vector<gene>, vector<match>, int);
//vector<gene> gene_fusion(gene,vector<gene>,double*,int,bool);
vector<gene> gene_fusion(gene,vector<gene>,int,bool);
vector<gene> gene_tri_fusion(gene,vector<gene>,int,bool);
vector<gene> gene_merge(vector<gene>,vector<gene>);
vector<int> match_merge(vector<int>,int,int);
double* map_to_array(map<int,double>,double*);
//double* map_to_array(map<int,double>,double*,int*);
int get_min_freq(vector<double>);
int rex(int);
bool qmode(string);
int qfind(vector<gene>, string);
m_point v_sum(vector<m_point>,int);
double dev_std(double, int, vector<m_point>);
float percento(short int, short int);
void tab_init();
const string tof(bool);
int v_pos_max(vector<m_point>);
bool merge_check(vector<int>, int);
//long double get_pvalue(int,m_point,double);
//long double get_new_pvalue(int,int,m_point,int);
//int get_gaps(string, m_point*);
//m_point get_coeff(int,int);
vector<gene> v_gene_shrink(vector<gene>);
bool gene_unique(gene,vector<gene>);
double factorial(double);
bool iso_check(string);
bool merge_gap_check(vector<string>);
vector<vector <gene> > sep_genes(vector<gene>);
vector<short int> exF_calculus(vector<int>*,int,int,int,int);
void frame_freq(vector<gene>*, m_point*);
map<int,double> write_new_freqfile(map<int,double>, vector<int> *, int, int, m_point*);
vector<string> string_to_vector(string);
bool fam_check(vector<string>);
vector<int> dup_wash(vector<int>, vector<gene>, vector<int>);
string launch_blastp(string, string, string, string, bool);
string launch_blastn(string, string, string, string, bool);
string seq_find(string,string,bool);
int frame_min_freq(m_point*);
unsigned short int blast_hit_counter(string);
vector<string> post_process_merger(vector<string>,gene*,gene*);
string brack_clean(string,int*);
void name_change(gene*,string,int);
bool merged_with_check(string,vector<string>*);
#ifdef DEBUG10
void gene_view(gene);
#endif

int main(int argc, char *argv[])
{
	vector<vector <gene> > qfg;
	vector<gene> genes, Qgenes, qgenes;
	vector<match> matches, fus_matches1, fus_matches2;
	vector<vector <int> > mqmatches;	
	vector<map <int,double> > Fr;
	vector<Frame_Freq> V_Framefreq;
	genesfiles.push_back(Hs); //new
	clineparser(argc, argv);	
	ofstream out;
	int fc = 0;

	if(WEB || QID.size() != 0)
	{
		tmps = "ex_work_dir/ex.out.";	
		tmps += QID;
	//	nice(2);
	}

	cerr << "\nParsing genes file...";
	for(int gf = 0; gf < genesfiles.size(); gf++)
	{
		vector<gene> temp_genes;
		temp_genes = parser(genesfiles[gf].c_str(), gf,0,INF);
		genes = gene_merge(genes,temp_genes);
	}

	cerr << "done" << endl;

	for(int gf = 0; gf < genesfiles.size(); gf++)
	{
		bool oldfreq = freq;
		m_point framefreq[5];
	//	vector<m_point> Framefreq[3];
		genesfile = genesfiles[gf];
		freqfile = freqfiles[gf];

		if(!freq)
			Fr.push_back(getfreq(framefreq,false));

		Frame_Freq one_more_framefreq;

		one_more_framefreq.filenum = gf;
	
		if(freq)
		{
			freqgen(genes.begin(), genes.end(),gf);

			getfreq(framefreq,true);

			Fr.push_back(getfreq(framefreq,false));

			if(!query)
				exit(EXIT_SUCCESS);
		}

		for(int i = 0; i < 5; i++)
			one_more_framefreq.Framefreq.push_back(framefreq[i]);

		V_Framefreq.push_back(one_more_framefreq);

		freq = oldfreq;
		cerr << endl;
	}


	#ifdef DEBUG7
	cerr << "Fr.size() = " << Fr.size() << endl; 
	for(int m1 = 0; m1 < Fr.size(); m1++)
		cout << m1 << '\t' << Fr[m1] << endl;
	#endif

	int qf = qfind(genes,queryfile);
	
	if(qf == -1)
	{
		cerr << "Parsing query file...";
		Qgenes = parser(queryfile.c_str(),0,5,200);
		cerr << "done" << endl;
	}
	else
		Qgenes.push_back(genes[qf]);

//	if(gensel && querysel)
	for(int r = 0; r < Fr.size(); r++)
	{
		m_point framefreq[5];
		string old_genes_file = genesfiles[r];

		if(!custom)
			genesfile = queryfile;

		else if(custom && Qgenes[0].fileorg != "Other")
			genesfile = Qgenes[0].fileorg;
			
	//	freqfile = freqfiles[1];
		freqfile = queryfile;
		freqfile += ".freq";

		map<int,double>  qfr;

		qfr = getfreq(framefreq, false);

		for(int i = 0; i < 5; i++)
		{
			V_Framefreq[r].Framefreq[i] += framefreq[i];
			V_Framefreq[r].Framefreq[i] /= 2;
		}

		for(int i = 0; i < qfr.size(); i++)
		{
			Fr[r][i] += qfr[i];
			Fr[r][i] /= 2;
		}

		genesfile = old_genes_file;
	}

	if(gensel)
	{
		qf = qfind(genes,gsel);
		
		if(qf == -1)
		{
			cerr << "Can't find " << gsel << endl;
			exit(EXIT_FAILURE);
		}

		gene gS = genes[qf];
		genes.clear();
		genes.push_back(gS);
	}

	if(querysel && Qgenes.size() > 1)
	{
		qf = qfind(Qgenes,qsel);

		if(qf == -1)
		{
			cerr << "Can't find " << qsel << " in " << queryfile << endl;
			exit(EXIT_FAILURE);
		}
		
		gene qS = Qgenes[qf];
		Qgenes.clear();
		Qgenes.push_back(qS);
	}

	//OUTINIT
	if(outfile.size() > 0 && !outapp)
	{
		out.open(outfile.c_str(),ios::out);
		out /*<< "Exons length database file: " << genesfile */ << " (" << genes.size()
                    << " genes)" << endl << "Query file: " << queryfile << " (" << qgenes.size() << ")"
		    << endl << "Alignment type : ";
		if(local)
			out << "local";
		else
			out << "global";
		out << "\nOpen gap penalty = " << OGAP << "\nContinuing gap penalty = " << CGAP
                    << "\nPre filter : " << FILTER << "\nStart exons skipped: " << S_SKIP << "\nEnd exons skipped: " << E_SKIP
		    << "\nTolerance : " << TOL << endl; 
	}
	//END

	if(table)
		tab_init();

	if(outfile.size() > 0)
	{
		out.close();
		out.open(outfile.c_str(),ios::app);
	}

//	fus_genes = fus_generator(genes);

	if(Qgenes.size() < 1)
	{
		cerr << endl << "I need a query to work..." << endl;
		exit(EXIT_FAILURE);
	} 

	else if(Qgenes.size() == 1)
		cerr << "\nWorking on " << genes.size() << " genes and " << Qgenes.size() << " query." << endl;
	else
		cerr << "\nWorking on " << genes.size() << " genes and " << Qgenes.size() << " queries." << endl;

	for(int Q = 0; Q < Qgenes.size(); Q++)
	{
		qfg.clear();
		qgenes.clear();
		matches.clear();
		fus_matches1.clear();
		fus_matches2.clear();
		mqmatches.clear();
		qgenes.push_back(Qgenes[Q]);

/*		if(qmerge)
		{
			fus_qgenes = fus_generator(qgenes);
			fus_matches1 = sizecomp(genes.begin(), genes.end(), fus_qgenes.begin(), fus_qgenes.end());
		}*/

		matches = sizecomp(genes.begin(), genes.end(), qgenes.begin(), qgenes.end());
		
//		if(gmerge)
//			fus_matches2 = sizecomp(fus_genes.begin(), fus_genes.end(), qgenes.begin(), qgenes.end());	


		for(int i = 0; i < qgenes.size(); i++)
			mqmatches.push_back(match_merge(i,matches.begin(),matches.end(),fus_matches1.begin(),
						fus_matches1.end(),fus_matches2.begin(),fus_matches2.end())); 

		matches.clear();
		fus_matches1.clear();

		fus_matches2.clear();

		#ifdef DEBUG4
		for(int db2 = 0; db2 < mqmatches.size(); db2++)
		{
			cout << "QUERY N " << db2+1 << " size = " << mqmatches[db2].size() << endl;
			for(int db3 = 0; db3 < mqmatches[db2].size(); db3++)
				cout << mqmatches[db2][db3] << endl;
			cout << endl << "****************" << endl;
		}
		#endif

		int nmq = 0;

		for(int i = 0; i < qgenes.size(); i++)
		{
			bool merged = true;

			if(!qgenes[i].merged && i != 0)
			{	
				nmq++;
				merged = false;
			}

			if(mqmatches[nmq].size() >= 20)
               			N20 = (mqmatches[nmq].size()/20);
       			else
               			N20 = 1;

			if(!WEB)
			cout << "\nNumber of genes to be aligned with query (" << qgenes[i].name << "): " 
             			<< mqmatches[nmq].size() << endl;	
			if(outfile.size() > 0)
				out << "\nNumber of genes to be aligned with query (" << qgenes[i].name << "): " 
		    			<< mqmatches[nmq].size() << endl;

			if(mqmatches[nmq].size() > 0)
//				if(!gmerge)
				main_aligner(genes.begin(),genes.end(),qgenes[i],mqmatches[nmq].begin(),mqmatches[nmq].end(),Fr.begin(),Fr.end(),&V_Framefreq); 
//					else
//						main_aligner(qfg[nmq].begin(),qfg[nmq].end(),qgenes[i],mqmatches[nmq].begin(),mqmatches[nmq].end(),Fr.begin(),Fr.end());
			else
				cerr << "\nSkipping...\n\n";

			if(!merged && nmq > 1)
				mqmatches[nmq - 1].clear();
		}

		if(outfile.size() > 0)
			out.close();
	}

	if(table)
		cerr << "\nTable file: " << tabfile << endl << endl;

	return EXIT_SUCCESS;
}

vector<gene> parser(const char* ifile, int filenum, int minex, int maxex)
{
	ifstream in(ifile);
	string line;
	vector<gene> genes;
	int count = 0; 

	if(!in)
	{
		if(!WEB)
			cout << "No such file: " << ifile << endl;
		exit(EXIT_FAILURE);
	}
	
	while(getline(in,line))
	{
		if(line[0] == '#' || line.size() == 0)
			continue;

		gene g1;
		g1.merged = false;
		g1.tri_merged = false;
		g1.filenum = filenum;

		string exonS,exonE;

		istringstream st1(line);

		st1 >> g1.name;

		#ifdef DEBUG1
		cout << "Parsing gene: " << g1.name << endl;
		#endif

		if(g1.name.find("CUSTOM_") == string::npos)
			st1 >> g1.chr >> g1.strand  >> g1.cdsS >> g1.cdsE >> g1.exN >> exonS >> exonE >> g1.name2 >> g1.protein_acc;

		else
		{
			g1.chr = "0";
			g1.strand = '+';
			custom = true;

			st1 >> g1.cdsS >> g1.cdsE >> g1.exN >> exonS >> g1.fileorg;

			if(querysel && gensel)
			{
				if(g1.fileorg != "Other")
				{
					freqfiles[1] = g1.fileorg;
					freqfiles[1] += ".freq";
				}
				else
					freqfiles[1] = freqfiles[0];
			}

			g1.protein_acc = g1.name;

		}

		if(g1.protein_acc.find("SINFRUP") != string::npos)
			g1.protein_acc = g1.name;

		g1.family = g1.name;
		g1.exL = exL_parser(exonS,exonE, g1.exN, g1.strand, g1.cdsS, g1.cdsE, &g1.exF, g1.cdsD);
		
		if(WEB)
			g1.ex_or_L = g1.exL;

		#ifdef HACK1
		if(g1.exL.size() >= minex && g1.exL.size() <= maxex)
			genes.push_back(g1);
		#else
		genes.push_back(g1);
               	#endif
	}

	return genes;
}

string linecleaner(string line)
{
	for(int i = 0; i < line.size(); i++)
	{
		for(short int sep = 0; sep < SEPSIZE; sep++)
			if(line[i] == SEPARATORI[sep])
				line[i] = ' ';
	}
	
	return line;
}

vector<int> exL_parser(string exonS,string exonE, int exN, char strand, int cdsS, int cdsE, vector<short int> *EXF, int *cdsD)
{
	vector<int> exL(exN);

	exonS = linecleaner(exonS);
	exonE = linecleaner(exonE);
	int fcE = -1, fcP = -1, lcE = -1, lcP = -1;

	if(!qmode(exonE))
	{
		istringstream st1(exonS), st2(exonE);

		if(strand == '+')
			for(int i = 0; i < exN; i++)
			{
				int exS,exE;
				st1 >> exS;
				st2 >> exE;
				exL[i] = exE - exS;
				#ifdef DEBUG2
				cout << "exS = " << exS << " exE = " << exE << " i = " << i << " exL[" << i << "] = " << exL[i] << endl;
				#endif
				if(exE >= cdsS && fcE < 0)
				{
					fcE = i;
					fcP = exL[i] - (exE - cdsS) + 1;
				}
				if(exE >= cdsE && lcE < 0)
				{
					lcE = i;
					lcP = exL[i] - (exE - cdsE) + 1;
				}
					
			}
		else if (strand == '-')
			for(int i = exN - 1; i >=0; i--)
			{
				int exS,exE;
                        	st1 >> exS;
				st2 >> exE;
                        	exL[i] = exE - exS;
				#ifdef DEBUG2
				cout << "exS = " << exS << " exE = " << exE << " i = " << i << " exL[" << i << "] = " << exL[i] << endl;
				#endif

				if(exE >= cdsE && fcE < 0)
                                {
                                        fcE = i;
                                        fcP = exL[i] - (exL[i] - (exE - cdsE)) + 1;
                                }
				if(exE >= cdsS && lcE < 0)
                                {
                                        lcE = i;
                                        lcP = exL[i] - (exL[i] - (exE - cdsS)) + 1;
                                }
			}

		else
			{
				cerr << "Something went wrong in the exons parsing process..." << endl;
				exit(EXIT_FAILURE);
			}	

		cdsD[0] = fcE;
		cdsD[1] = fcP;
		cdsD[2] = lcE;
		cdsD[3] = lcP;

		*EXF = exF_calculus(&exL, fcE, fcP, lcE, lcP); 

//		cerr << endl;

/*		for(int i = 0; i < EXF->size(); i++)
			cerr << EXF->at(i) << '\t';*/
	}

	else
	{
		istringstream st1(exonS);

		for(int i = 0; i < exN; i++)
		{
			int ex;
			st1 >> ex;
			if(ex <= 0)
			{
				cerr << "Something went wrong in the exons parsing process..." << endl;
				exit(EXIT_FAILURE);
			}
			exL[i] = ex;
		}

		if(cdsS < 0 || cdsE <= 0 || cdsE <= cdsS)
		{
			lcE = 1000;
			lcP = 1000;
			cdsD[0] = fcE;
                        cdsD[1] = fcP;
                        cdsD[2] = lcE;
                        cdsD[3] = lcP;
			*EXF = exF_calculus(&exL, fcE, fcP, lcE, lcP);
		}
		else	
		{
			int sum = 0, sum_1 = 0;
			bool sflag = false, eflag = false;
			for(int c = 0; c < exL.size(); c++)
			{
				sum += exL[c];

				if(cdsS <= sum && !sflag)
				{
					fcE = c;
					fcP = cdsS - sum_1;
					sflag = true;
				}

				if(cdsE <= sum && !eflag)
				{
					lcE = c;
					lcP = cdsE - sum_1;  
					eflag = true;
				}

				sum_1 += exL[c];
			}

			cdsD[0] = fcE;
                	cdsD[1] = fcP;
                	cdsD[2] = lcE;
                	cdsD[3] = lcP;

			*EXF = exF_calculus(&exL, fcE, fcP, lcE, lcP);
		}
	}
 
	return exL;
}

vector<short int> exF_calculus(vector<int> *exL, int fcE, int fcP, int lcE, int lcP)
{
	vector<short int> exF;
	int actL = 0, i = -1, cdsL = 0;
	short int nextfr;

	if(fcE == -1 && fcP == -1 && lcE == 1000 && lcP == 1000)
	{
		for(int v = 0; v < exL->size(); v++)
			exF.push_back(-1);

		return exF;
	}

	for(vector<int>::iterator exLi = exL->begin(); exLi < exL->end(); exLi++)
	{
		i++;	
//		cerr << endl << cdsL << '\t' << i << '\t' << *exLi;
		if(i < fcE)
			continue;
		else if (i == fcE && i < lcE)
		{
			cdsL = *exLi - fcP;
			continue;
		}
		else if(i == fcE && i == lcE)
		{
			cdsL = lcP - fcP;
			continue;	
		}
		else if(i > fcE && i < lcE)
		{
			cdsL += *exLi;
			continue;
		}
		else if(i > fcE && i == lcE)
		{
			cdsL += lcP;
			continue;
		}
		else if(i > fcE && i > lcE)
			break;

	}

//	cerr << "\ncdsL = " << cdsL << " fcE = " << fcE << " fcP = " << fcP << " lcE = " << lcE << " lcP = " << lcP;

/*	if(cdsL%3 != 0)
		return exF;*/

	i = -1;
		
	for(vector<int>::iterator exLi = exL->begin(); exLi < exL->end(); exLi++)
	{
		i++;

		if(i < fcE)
		{
			exF.push_back(-1);
			continue;
		}
		else if(i == fcE && cdsL > (*exLi-fcP))
		{
			actL = *exLi - fcP + 1;
			exF.push_back(actL%3);
			nextfr = actL%3; 
//			cerr << endl << fcP << '\t' << actL << '\t' << nextfr;
			if(nextfr == 1)
			{
				nextfr = 2;
				continue;
			}
			else if(nextfr == 2)
			{
				nextfr = 1;
				continue;
			}

			continue;
		}
		else if(i == fcE && cdsL <= (*exLi-fcP))
                {
                        exF.push_back(cdsL%3);
                        actL = cdsL;
                        nextfr = -1;
			continue;
                } 
		else if(i > fcE && cdsL > (actL + *exLi))
		{
			exF.push_back(nextfr);
			actL += *exLi;
			nextfr = (*exLi - nextfr)%3;
			if(nextfr == 1)
                        {
                                nextfr = 2;
                                continue;
                        }
                        else if(nextfr == 2)
                        {
                                nextfr = 1;
                                continue;
                        }
			
			continue;
		} 
		else if(i > fcE && cdsL <= (actL + *exLi))
		{
			exF.push_back(nextfr);
			actL = cdsL;
			nextfr = -1;
			continue;
		}
	}

	return exF;
}


/*vector<short int> exF_parser(string frames, int exN, char strand)
{
	vector<short int> exF(exN);

	frames = linecleaner(frames);
	istringstream st1(frames);
	#ifdef DEBUG2
	cout << "frames.size() = " << frames.size() << " frames = " << frames << endl;
	#endif

	if(strand == '+')
		for(int i = 0; i < exN; i++)
		 {
			short int frame;
			st1 >> frame;
			exF[i] = frame;
		}
	else
		for(int i = exN-1; i >= 0; i--)
		{
			short int frame;
			st1 >> frame;
			exF[i] = frame;
		}

	return exF;
}*/

vector<match> sizecomp(vector<gene>::iterator gS, vector<gene>::iterator gE, vector<gene>::iterator gqS, vector<gene>::iterator gqE)
{
	vector<gene> genes(gS,gE);
	vector<gene> qgenes(gqS,gqE);
	vector<match> matches;
	bool flag;

	for(int i1 = 0; i1 < genes.size(); i1++)
		for(int i2 = 0; i2 < qgenes.size(); i2++)
		{
			flag = false;
			int count = 0;
			for(int i3 = 0; i3 < genes[i1].exL.size(); i3++)
			{
				for(int i4 = 0; i4 < qgenes[i2].exL.size(); i4++)
				{
					if(genes[i1].exL[i3] == qgenes[i2].exL[i4] || !FILTER) 
					{
						count++;
						if(count == FILTER || !FILTER)
						{
							match m1;
							m1.pos[0] = i1;
							m1.pos[1] = i2;
							m1.exN[0] = i3;
							m1.exN[1] = i4;
							matches.push_back(m1);
							flag = true;
							break;
						}
					}
				}
				if(flag)
					break;
			}
		}
						
	return matches;
}

/*vector<gene> fus_generator(vector<gene> genes)
{
	for(int i1 = 0; i1  < genes.size(); i1++)
	{
		vector<int> fusion;
		for(int i2 = 0; i2 < genes[i1].exL.size() - 1; i2++)
			fusion.push_back(genes[i1].exL[i2] + genes[i1].exL[i2+1]);

		genes[i1].merged = true;
		genes[i1].exL = fusion;
		genes[i1].exN = fusion.size();
	}

	return genes;
}*/	

vector<int> match_merge(int qg, vector<match>::iterator mS, vector<match>::iterator mE, vector<match>::iterator f1mS, vector<match>::iterator f1mE, vector<match>::iterator f2mS, vector<match>::iterator f2mE)
{
	vector<match> matches(mS,mE);
	vector<match> f1matches(f1mS,f1mE);
	vector<match> f2matches(f2mS,f2mE);
	list<int> Lmatches;
	vector<int> merged;
		
	for(int i = 0; i < matches.size(); i++)
		if(matches[i].pos[1] == qg)
			Lmatches.push_back(matches[i].pos[0]);

	for(int i = 0; i < f1matches.size(); i++)
                if(f1matches[i].pos[1] == qg)
                        Lmatches.push_back(f1matches[i].pos[0]);


	for(int i = 0; i < f2matches.size(); i++)
                if(f2matches[i].pos[1] == qg)
                        Lmatches.push_back(f2matches[i].pos[0]);

	Lmatches.sort();
	Lmatches.unique();

	for(list<int>::iterator i = Lmatches.begin(); i != Lmatches.end(); i++)
		merged.push_back(*i);

	return merged;
}

vector<vector <gene> > gene_separator(vector<gene> genes)
{
	vector<vector <gene> > sep_genes;

	for(int i = 0; i < genesfiles.size(); i++)
	{
		vector<gene> g;
		sep_genes.push_back(g);
	}

	for(int i = 0; i < genes.size(); i++)
		sep_genes[genes[i].filenum].push_back(genes[i]);

	return sep_genes;
}	
	
void main_aligner(vector<gene>::iterator gS, vector<gene>::iterator gE, gene qgene, vector<int>::iterator mS, vector<int>::iterator mE,
                  vector< map<int,double> >::iterator frI, vector< map<int,double> >::iterator frE, vector<Frame_Freq> *V_framefreq) 
{
	vector<gene> all_genes(gS,gE); 
	vector<vector <gene> > sep_genes;

	sep_genes = gene_separator(all_genes);
	int align_count_check = 0, tmpss, tmpes;

	for(int z = 0; z < sep_genes.size(); z++)
	{
		vector<gene> genes, qmerged, bestg, bmerged;
//		vector<int> matches(mS,mE); 
		vector<int> mpos,bpos,bbpos;
		vector<int> mmpos[2];
		vector<m_point> score, Score;
		vector < map<int,double> > Fr(frI,frE);
		double *FrV, AvgScore, DevStd, min_freq;
	//	long double Pvalue;
		long double vvec[3] = {0,0,0};
		int more_al, moreN20, more_count = 0, oss = S_SKIP, oes = E_SKIP, merge_count_check = 0;
		m_point Framefreq[5];
		enum{Q,G};
		
		genes = sep_genes[z];

		if(genes.size() > 20)
			N20 = (int)all_genes.size()/20;
		else
			N20 = 1;
		
		if((int)qgene.exL.size() - S_SKIP - E_SKIP < 1)
		{
			S_SKIP = 0;
			E_SKIP = 0;
		}

		FrV = map_to_array(Fr[0],&min_freq);

		for(int q = 0; q < 5; q++)
			Framefreq[q] = V_framefreq->at(0).Framefreq[q]; 

		cout.precision(4);
	
//		if(!qmerge)
		qmerged.push_back(qgene);

		ofstream out;
		ofstream wout;

		if(outfile.size() > 0)
			out.open(outfile.c_str(), ios::app);

		cerr << "\nAligning";

		if(WEB && z == 0)
		{
			wtmpfile = "ex_work_dir/.";
			wtmpfile += QID;
			wtmpfile += ".tmp";
			wout.open(wtmpfile.c_str());
			wout << "#1\t" << all_genes.size() << endl;
			wout.close();
		}

		for(int i = 0; i < genes.size(); i++)
		{
			if(i%N20 == 0 && i > 0 || (i == all_genes.size()-1 && i/N20 < 20))
			{
				if(WEB && align_count_check < 20)
				{
	//				wout.open(wtmpfile.c_str(),ios::app);
	//				wout << "<img src=\"../gif/adv.gif\" hspace=0 >\n";
					align_count_check++;
	//				wout.close();
				}
				
				else		
					cerr << '.';
			}

			if(i > 0)
				if(genes[i].filenum != genes[i-1].filenum)
				{
					delete[] FrV;
					FrV = map_to_array(Fr[genes[i].filenum],&min_freq);

					for(int q = 0; q < 5; q++)
						Framefreq[q] = V_framefreq->at(genes[i].filenum).Framefreq[q];
				}			

			genes[i].or_score = aligner(&qgene, &genes[i],false,FrV,vvec,min_freq,Framefreq,false);
			score.push_back(genes[i].or_score);
			mpos.push_back(i);
		}
	
		cerr << "done" << endl;

		if(score_out)
		{
			ofstream sout("scores.txt",ios::app);
			sout << '-' << endl;
			for(int sc = 0; sc < score.size(); sc++)
				sout << score[sc] << endl;
			sout.close();
		}

/*		if(WEB && align_count_check < 20 && z == (int)genesfiles.size() - 1)
		{
			wout.open(wtmpfile.c_str(),ios::app);

			for(int ac = 0; ac < (20 - align_count_check); ac++)
				wout << "<img src=\"../gif/adv.gif\" hspace=0 >\n";

			wout.close();
		}*/


//		AvgScore = (double)v_sum(score, score.size())/(double)genes.size();	
//		DevStd = dev_std(AvgScore,genes.size(),score);

		if(!blasthitonly) //|| !WEB)
		{
			bpos = best_pos(score, BEST);

			for(int ap = 0; ap < bpos.size(); ap++)
				bestg.push_back(genes[bpos[ap]]);
		}
		
		else
		{
			bpos = best_pos(score, 2*BEST);

			for(int z = 0; z < bpos.size(); z++)
			{
				string tabp, tabn;
				unsigned short int phit,nhit;

				tabp = launch_blastp(qgene.protein_acc,queryfile,genes[bpos[z]].protein_acc,genesfiles[genes[bpos[z]].filenum], false);	
				tabn = launch_blastn(qgene.family,queryfile,genes[bpos[z]].family,genesfiles[genes[bpos[z]].filenum], false);

				phit = blast_hit_counter(tabp);			
				nhit = blast_hit_counter(tabn); 

				if(phit || nhit)
					bestg.push_back(genes[bpos[z]]);

				if(z == bpos.size() - 1)
				{
					ostringstream rm1, rm2;
					rm1 << "rm " << tabp;
					rm2 << "rm " << tabn;

					system(rm1.str().c_str());
					system(rm2.str().c_str());	
				}
			}

			if(bestg.size() == 0)
				bestg.push_back(genes[bpos[0]]);
		}

//		cerr << "\nbestg.size() = " << bestg.size();

		if(MAX_CYCLE > 0)
		{
			cerr << "\nMerging exons step...";

			vector<gene> qqmerged = qmerged;

		//	qlmerged.push_back(qgene);
			for(int d = 0; d < bestg.size(); d++)
			{
				int sc = 0, sc2 = 0, max_cycle = MAX_CYCLE;
                		bool fflag = false;	

				if(bestg[d].exN <= S_SKIP + E_SKIP)
                        	{
                                	tmpss = S_SKIP;
                              	        tmpes = E_SKIP;

                           	        S_SKIP = 0;
                                        E_SKIP = 0;
                       		 }

				while(max_cycle)
               			{
					vector<gene> bufg,bufq;
					bufg.push_back(bestg[d]);
					bufq = qqmerged;

					bufq = gene_fusion(bestg[d],bufq,sc2,fflag);

					if(max_cycle == MAX_CYCLE)
                                        {
                                                vector<gene> tbufq;
				
						tbufq = gene_merge(qqmerged,bufq);
                                                tbufq = gene_tri_fusion(bestg[d],tbufq,sc2,fflag);

                                                for(int z = 0; z < tbufq.size(); z++)
                                                        bufq.push_back(tbufq[z]);
                                        }

			
                        		bufq = v_gene_shrink(bufq); //xxx
					
					vector<gene>::iterator bqi= bufq.begin();

					for(int zq = 0; zq < bufq.size(); zq++)
					{
						m_point tsc;

						tsc = aligner(&bufq[zq],&bestg[d],false,FrV,vvec,min_freq,Framefreq,true);

						if(tsc == (m_point)-INF)
						{
							bufq.erase(bqi);
							zq--;
						}
						else
							bqi++;
					}

                        		if(sc == bufq.size() || bufq.size() - sc > MAX_MERGE)
                                		break;

                        		else if(sc == 0 && !fflag)
                                		fflag = true;

                        		else if(sc >= 0 && fflag)
                        		{
                                		sc2 = sc;
                                		sc = bufq.size();
                                	}

					qqmerged = gene_merge(qqmerged,bufq);

                        		max_cycle--;

					qqmerged = v_gene_shrink(qqmerged);

		/*			cerr << endl << "qqmerged:" << endl;
					for(int i = 0; i < qqmerged.size(); i++)
        	                		cerr << endl << qqmerged[i].name;
					cerr << "\n***" << endl;*/
                		}


				qmerged = gene_merge(qqmerged,qmerged);

				qmerged = v_gene_shrink(qmerged);
		
				if(bestg[d].exN <= S_SKIP + E_SKIP)
                                {
                                        S_SKIP = tmpss;
                                        E_SKIP = tmpes; 
                                }

			}
		} 


//		for(int i = 0; i < qmerged.size(); i++)
//			cerr << endl << qmerged[i].name;

		vector<gene> bufbest = bestg;
		int sc = 0, sc2 = 0, max_cycle = MAX_CYCLE;
		bool fflag = false;

       		while(max_cycle)
       		{

       			bestg = gene_fusion(qgene,bestg,sc2,fflag);

			if(max_cycle == MAX_CYCLE)
			{
				vector<gene> tbestg;

				tbestg = gene_merge(bufbest,bestg);

				tbestg = gene_tri_fusion(qgene,tbestg,sc2,fflag);

				for(int z = 0; z < tbestg.size(); z++)
					bestg.push_back(tbestg[z]);
			}

			bestg = v_gene_shrink(bestg); //xxx

			vector<gene>::iterator bgi = bestg.begin();

			for(int zq = 0; zq < bestg.size(); zq++)
			{
				m_point tsc;

				tsc = aligner(&qgene,&bestg[zq],false,FrV,vvec,min_freq,Framefreq,true);	


				if(tsc == -INF)
				{
	//				cerr << "\n\ngerasing: " << bestg[zq].name << endl << endl;
//					cerr << "\n\ngerasing: " << bgi->name << endl << endl;

					bestg.erase(bgi);
					zq--;
				}
				else
					bgi++;

//				cerr << "\nbestg.size()= " << bestg.size() << endl;
//				cerr << "\nmax_cycle = " << max_cycle << endl;
			}		

     			if(sc == bestg.size() || (bestg.size() - sc > MAX_MERGE))
     			 	break;

       			else if(sc == 0 && !fflag)
       				fflag = true;

       			else if(sc >= 0 && fflag)
       			{
       				sc2 = sc;
       				sc = bestg.size();
      			}

			max_cycle--;
       		}

		bestg = gene_merge(bufbest,bestg);


//		cerr << "\nbestg.size() = " << bestg.size();

		bufbest.clear();

		cerr << "done" << endl;

	#ifdef DEBUG10
	for(int dbx = 0; dbx < bestg.size(); dbx++)
                gene_view(bestg[dbx]);

	cerr << "\nbestg.size() = " << bestg.size() << endl;
	#endif

		bestg = v_gene_shrink(bestg);


//		cerr << "\nbestg.size() = " << bestg.size();

	#ifdef DEBUG10
	for(int dbx = 0; dbx < bestg.size(); dbx++)
                gene_view(bestg[dbx]);

	cerr << "\nbestg.size() = " << bestg.size() << endl;
	#endif

		more_al = qmerged.size() * bestg.size();
		int adjust = 0;

		if(more_al > 20)
		{
			moreN20 = more_al / 20;
		}
		else
			moreN20 = 1;

		delete[] FrV;
		FrV = map_to_array(Fr[bestg[0].filenum],&min_freq);
		for(int q = 0; q < 5; q++)
			Framefreq[q] = V_framefreq->at(bestg[0].filenum).Framefreq[q];

		if(WEB && genesfiles.size() <= 1)
       	 	{
       	 		wout.open(wtmpfile.c_str(),ios::app);
       	 		wout << "#2\t" << more_al << endl;
       	 		wout.close();
       	 	}
		else if(WEB && genesfiles.size() > 1 && z == (int)genesfiles.size() - 1)
		{
			wout.open(wtmpfile.c_str(),ios::app);
                        wout << "#2\tMULTIPLE" << endl;
                        wout.close();
		}

		for(int bp = 0; bp < bestg.size(); bp++)
		{
			gene Qgene;
			vector<m_point> sscore;
			int psscore;

			if(bp > 0)
				if(bestg[bp].filenum != bestg[bp-1].filenum)
				{
					delete[] FrV;
					FrV = map_to_array(Fr[bestg[bp].filenum],&min_freq);
					for(int q = 0; q < 5; q++)
						Framefreq[q] = V_framefreq->at(bestg[bp].filenum).Framefreq[q];
				}

			if(bp == 0 && MAX_CYCLE > 0)
			{
				if(!WEB)	
					cerr << "\nPerforming " << more_al << " more alignments";
			}

			for(int t0 = 0; t0 < qmerged.size(); t0++)
			{
				if((!qmerged[t0].merged && !qmerged[t0].tri_merged) && (!bestg[bp].merged && !bestg[bp].tri_merged) ||
				    (merged_with_check(qmerged[t0].name,&bestg[bp].merged_with)) || merged_with_check(bestg[bp].name,&qmerged[t0].merged_with))
						sscore.push_back(aligner(&qmerged[t0], &bestg[bp], false, FrV,vvec,min_freq,Framefreq,false)); //NEW
				else
		
						sscore.push_back(-INF);


//				cerr << endl << "Aligning " << qmerged[t0].name << " with " << bestg[bp].name << " score: " << sscore[sscore.size()-1];
				more_count++;

				if(more_count % moreN20 == 0 && MAX_CYCLE > 0 && genesfiles.size() <= 1)
				{
					if(WEB && merge_count_check < 20 && genesfiles.size() <= 1)
					{
					/*	wout.open(wtmpfile.c_str(),ios::app);
						wout << "<img src=\"../gif/adv2.gif\" hspace=0 >\n";
						merge_count_check++;
						wout.close();*/
					}
						
					else
                               			cerr << '.';
				}
			}


			psscore = v_pos_max(sscore);
			Qgene = qmerged[psscore];
			mmpos[Q].push_back(psscore);
			mmpos[G].push_back(bp);

			if(bp == bestg.size() - 1 && MAX_CYCLE > 0)
				cerr << "done" << endl;

//			else
//				Qgene = qgene;
		
//			if(qmerge && gmerge && !WEB)
  //              	        cout << "\nScore = " << score[bpos[bp]];

//			if(!gensel && (qmerge || gmerge) && !WEB)
//				cout << " Avg.Score = " << AvgScore << " Dev.Std = " << DevStd << endl;

//			else if((qmerge || gmerge) && !WEB)
//				cout << endl;

			if(outfile.size() > 0)
				out << "\nScore = " << score[bpos[bp]] << endl;

//			if(qmerge || gmerge)
//				aligner(Qgene,genes[bpos[bp]],true,FrV,vvec,min_freq,Framefreq);
//			else
			Score.push_back(aligner(&Qgene,&bestg[bp],false,FrV,vvec,min_freq,Framefreq,false));
		}

	
/*		if(WEB)
		{
			if(merge_count_check < 20 && genesfiles.size() <= 1)
			{
				wout.open(wtmpfile.c_str(),ios::app);
		
				for(int cc = 0; cc < (20 - merge_count_check); cc++)
               	                	wout << "<img src=\"../gif/adv2.gif\" hspace=0 >\n";
		
				wout.close();
			}
		} */

//		if(!qmerge && !gmerge)
//		{
		bbpos = best_pos(Score, (int)Score.size());

//		cerr << endl << "bbpos.size() = " << bbpos.size() << endl;
//		cerr << endl << "bestg.size() = " << bestg.size() << endl;


		bbpos = dup_wash(bbpos, bestg, mmpos[G]);

//		cerr << endl;//DEBUG
//		for(int i = 0; i < bestg.size(); i++)//DEBUG
//			cerr << bestg[i].name << endl;//DEBUG

//		cerr << endl << "bbpos.size() = " << bbpos.size() << endl;

//		cerr << endl;//DEBUG
  //              for(int i = 0; i < bbpos.size(); i++)//DEBUG
    //                    cerr << bestg[bbpos[i]].name << endl;//DEBUG
		
		int old_vbest = V_BEST, act_vbest;

		if(V_BEST > bbpos.size())
			V_BEST = bbpos.size();

		act_vbest = V_BEST;

		delete[] FrV;
		FrV = map_to_array(Fr[bestg[mmpos[G][bbpos[0]]].filenum],&min_freq);
		for(int q = 0; q < 5; q++)
			Framefreq[q] = V_framefreq->at(bestg[mmpos[G][bbpos[0]]].filenum).Framefreq[q];

		for(int jk = 0; jk < act_vbest; jk++)
		{
		//	double df = 1, tmp;
			double tmp;
			ostringstream tout;
			tout.precision(4);

			if(jk == (int)bbpos.size() - 1)
				last_one = true;

			if(jk > 0)
				if(bestg[mmpos[G][bbpos[jk]]].filenum != bestg[mmpos[G][bbpos[jk-1]]].filenum)
				{
					delete[] FrV;
					FrV = map_to_array(Fr[bestg[mmpos[G][bbpos[jk]]].filenum],&min_freq);
					for(int q = 0; q < 5; q++)
						Framefreq[q] = V_framefreq->at(bestg[mmpos[G][bbpos[jk]]].filenum).Framefreq[q];
				}
				
			if(!WEB)
				tout << "\nScore = " << Score[bbpos[jk]] << " ( from " << bestg[mmpos[G][bbpos[jk]]].or_score;

			if(!gensel && !WEB) 
				tout << " Avg.Score = " << AvgScore << " Dev.Std = " << DevStd << " )" << endl;
			else if(!WEB)
				tout << " )" << endl;

			vvec[1] = AvgScore;
			vvec[2] = DevStd;

			if(bestg[mmpos[G][bbpos[jk]]].exN <= S_SKIP + E_SKIP)
			{
				tmpss = S_SKIP;
				tmpes = E_SKIP;

				S_SKIP = 0;
				E_SKIP = 0;
			}

//			cerr << endl << "---" << bestg[mmpos[G][bbpos[jk]]].name;

			tmp = aligner(&qmerged[mmpos[Q][bbpos[jk]]],&bestg[mmpos[G][bbpos[jk]]],true,FrV,vvec,min_freq,Framefreq,false);

			if(bestg[mmpos[G][bbpos[jk]]].exN <= S_SKIP + E_SKIP)
                        {
                                S_SKIP = tmpss;
                                E_SKIP = tmpes;
                        }

			if(tmp == -INF && jk < (int)bbpos.size() - 1)
			{
				if(act_vbest < (int)bbpos.size() - 1)
					act_vbest++;	

				continue;
			}
			else if(tmp == -INF && jk == (int)bbpos.size() - 1)
				continue;
			else
			{
				tout << "Evalue = " << vvec[0] << endl;	

				if(!WEB)
					cout << tout.str() << "**************" << endl << endl;
			}
		}

		V_BEST = old_vbest;
//	}
		
		mpos.clear();
		bpos.clear();
		score.clear();
		mmpos[Q].clear();
		mmpos[G].clear();
		bbpos.clear();
		Score.clear();
		bestg.clear();

		delete[] FrV;

		if((int)qgene.exL.size() - oss - oes < 1)
		{
			S_SKIP = oss;
			E_SKIP = oes;
		}
	}

	ofstream wfileout(tmps.c_str(),ios::out);
	wfileout << wfile.str();
	wfileout.close();
	
	return;
}

m_point aligner(gene *q, gene *g, bool flag, double *FrV, long double *vvec, double min_freq, m_point *Framefreq, bool bflag)
{
	vector<vector <int> > s(2);
	vector<vector <short int> > f(2);

	if(dynamic_gap)
	{
		int ml, ng;
		if(q->exL.size() < g->exL.size())
			ml = q->exL.size();
		else
			ml = g->exL.size();

		if(ml > 0 && ml < 15)
			ng = -2;
		else if(ml >= 15 && ml < 30)
			ng = -3;
		else if(ml >= 30)
			ng = -4;

		OGAP = CGAP = ng;
	}
	
	s[0] = remFL(q->exL);
	s[1] = remFL(g->exL);	
	f[0] = remFL(q->exF);
	f[1] = remFL(g->exF);

	m_point **matrix[3], ***patrix[3], score;

        int mx, mpos[3];
	
	x = (int)s[0].size() + 1;
        y = (int)s[1].size() + 1;

	for(short int m = 0; m < 3; m++)
        {
                matrix[m] = new m_point*[x];
                patrix[m] = new m_point**[x];

                for(int a = 0; a < x; a++)
                {
                        matrix[m][a] = new m_point[y];
                        patrix[m][a] = new m_point*[y];
                }
        }

	m_init(matrix,patrix);
	m_filler(matrix,patrix,s,FrV,min_freq,f,Framefreq);

	m_pmax(matrix, x-1, y-1, &mx);
	
	if(!local)
		score = matrix[mx][x-1][y-1];
	else
		score = m_max(matrix,mpos);

	#ifdef DEBUG5
	cout << "Query = " << q.name << " Gene = " << g.name << " Score = " << matrix[mx][x-1][y-1] << endl;
	#endif	

	if(!flag && !bflag)	
	{
		for(short int m = 0; m < 3; m++)
		{
			for(int a = 0; a < x; a++)
			{
				delete[] matrix[m][a];	
				delete[] patrix[m][a];
			}

			delete[] matrix[m];
			delete[] patrix[m];
		}

		s[0].clear();
		s[1].clear();
		s.clear();
		return score;
	}

	else

	{
		vector<string> rs(3);
		
		if(!local)
		{
			mpos[MX] = 0;
			mpos[X] = x - 1;
			mpos[Y] = y - 1;
		}
		else
			m_max(matrix,mpos);

		rs = m_traceback(matrix,patrix,s[0],s[1],mpos,q->merg_pos,g->merg_pos,vvec,g->filenum,&q->exF,&g->exF, Framefreq);

		if(bflag)
		{
			vector<string> q_a, g_a, a_a;

//			q->m_flag.clear();
//			g->m_flag.clear();

			q_a = string_to_vector(rs[0]);
			a_a = string_to_vector(rs[1]);
			g_a = string_to_vector(rs[2]);

		/*	for(int i = 0; i < a_a.size(); i++)
			{
				if(a_a[i] == "G")
				{
					if(q_a[i] != "-")
						q->m_flag.push_back(false);
					if(g_a[i] != "-")
						g->m_flag.push_back(false);
			
					continue;
				}
				else if(a_a[i] == "|")
				{
					q->m_flag.push_back(true);
					g->m_flag.push_back(true);				
					continue;
				}
				else
				{
					q->m_flag.push_back(false);
                                        g->m_flag.push_back(false);
				}
			}*/

//			cerr << endl << q->name << '\t' << g->name << endl;

			for(int i = 0; i < a_a.size(); i++)
			{
			//	cerr << endl << q_a[i] << '\t' << g_a[i] << '\t' << a_a[i];
				if((q_a[i].find("M") != string::npos || g_a[i].find("M") != string::npos) && a_a[i] != "|")
				{
				//	cerr << "\nqai = " << q_a[i] << endl << "gai = " << g_a[i] << endl << "aai = " << a_a[i] << endl;
					score = -INF;
				//	cerr << '\t' << "score = " << score;
				}
			}

			for(short int m = 0; m < 3; m++)     //NEW
                	{
                        	for(int a = 0; a < x; a++)
                        	{
                                	delete[] matrix[m][a];
                                	delete[] patrix[m][a];
                       	 	}

                        	delete[] matrix[m];
                        	delete[] patrix[m]; //END
                	}
	
			return score;  
		}

		if(!merge_gap_check(rs) && (at_least_one || !last_one))
		{
		//	cerr << "\nlast_one = " << last_one;
		//	cerr << "\nat_least_one = " << at_least_one;
			for(short int m = 0; m < 3; m++)
            	        {
                        	for(int a = 0; a < x; a++)
                        	{
                                	delete[] matrix[m][a];
                                	delete[] patrix[m][a];
                        	}

                        	delete[] matrix[m];
                        	delete[] patrix[m];
                	}

			return -INF;
		}

		#ifdef DEBUG6
		cout << "Query = " << q.name << " gene = " << g.name << endl 
		     << "rs[0] = " << rs[0] << endl << "rs[1] = " << rs[1] << endl << "rs[2] = " << rs[2] << endl;

		for(int db1 = 0; db1 < q.exL.size(); db1++)
			cout << q.exL[db1] << '\t';
		cout << endl;
 		#endif

		ofstream out;

		if(outfile.size() > 0)
			out.open(outfile.c_str(),ios::app);

		if(gensel && querysel)
                                rs = post_process_merger(rs,q,g);
		
		if(!WEB)
			cout << "\nQuery = " << q->name << " Gene = " << g->name << endl;
		else
		{
		//	ofstream wfile(tmps.c_str(),ios::app);	
			wfile << "->" << q->name << "#" << q->name2 << "\n->" << g->name << "#" << g->name2 << endl;
		//	wfile.close();
		}

		if(outfile.size() > 0)
			out << "\nQuery = " << q->name << " Gene = " << g->name << endl;


		s_out(matrix[0],rs,mpos,*q,*g,vvec[0]);	

		at_least_one = true;

		for(short int m = 0; m < 3; m++)
		{
                        for(int a = 0; a < x; a++)
                        {
                                delete[] matrix[m][a];
                                delete[] patrix[m][a];
                        }

			delete[] matrix[m];
			delete[] patrix[m];
		}

		rs.clear();
//		rs[0].clear();
//		rs[1].clear();
//		rs[2].clear();
		return 0;
	}
}

void m_init(m_point ***matrix, m_point ****patrix)
{
	for(short int m = 0; m < 3; m++)
        {
                matrix[m][0][0] = 0;
                patrix[m][0][0] = &matrix[m][0][0];
        }

        for(int i = 1; i < x; i++)
        {
                matrix[0][i][0] = -INF;
                patrix[0][i][0] = &matrix[0][i-1][0];
		if(!semilocal && !local)
                	matrix[1][i][0] = OGAP+((i-1)*CGAP);
		else
			matrix[1][i][0] = 0; 
                patrix[1][i][0] = &matrix[1][i-1][0];
                matrix[2][i][0] = -INF;
                patrix[2][i][0] = &matrix[2][i-1][0];
        }

        for(int i = 1; i < y; i++)
        {
                matrix[0][0][i] = -INF;
                patrix[0][0][i] = &matrix[0][0][i-1];
                matrix[1][0][i] = -INF;
                patrix[1][0][i] = &matrix[1][0][i-1];
		if(!semilocal && !local)
                	matrix[2][0][i] = OGAP+((i-1)*CGAP);
		else
			matrix[2][0][i] = 0; 
                patrix[2][0][i] = &matrix[2][0][i-1];
        }

        return;
}

void m_filler(m_point ***matrix, m_point ****patrix, vector<vector <int> > s, double *FrV, double min_freq, 
	      vector<vector <short int> > f, m_point *Framefreq)
{
	int frame_mf = frame_min_freq(Framefreq);

        for(int iy = 1; iy < y; iy++)
        {
                for(int ix = 1; ix < x; ix++)
                {
                    //    list<m_point> L;
			vector<m_point> L;
                        int mx;
			m_point check_res;

                        if(local)
                                L.push_back(0);

			check_res = c_check(s[0][ix-1],s[1][iy-1],FrV,f[0][ix-1],f[1][iy-1],Framefreq,frame_mf,min_freq);

                        L.push_back(matrix[0][ix-1][iy-1] + check_res);
                        L.push_back(matrix[1][ix-1][iy-1] + check_res);
                        L.push_back(matrix[2][ix-1][iy-1] + check_res);

                        if(!local)
                        {
                                matrix[0][ix][iy] = max(L,&mx);
                                patrix[0][ix][iy] = &matrix[mx][ix-1][iy-1];
                        }

                        else
                        {
                                matrix[0][ix][iy] = lmax(L,&mx);
                                if(matrix[0][ix][iy] > 0)
                                        patrix[0][ix][iy] = &matrix[mx][ix-1][iy-1];
                                else
                                        patrix[0][ix][iy] = &matrix[mx][ix][iy];
                        }

                        L.clear();

			if((iy < y - 1 && ix < x - 1) || (!semilocal && !local))
			{
                        	L.push_back(matrix[0][ix-1][iy] + OGAP);
                       		L.push_back(matrix[1][ix-1][iy] + CGAP);
                        	L.push_back(matrix[2][ix-1][iy] + OGAP);
			}
			else
			{
				L.push_back(matrix[0][ix-1][iy]);
                                L.push_back(matrix[1][ix-1][iy]);
                                L.push_back(matrix[2][ix-1][iy]);
			}

                        matrix[1][ix][iy] = max(L,&mx);
                        patrix[1][ix][iy] = &matrix[mx][ix-1][iy];
                        L.clear();

			if((iy < y - 1 && ix < x - 1) || (!semilocal && !local))
			{
				L.push_back(matrix[0][ix][iy-1] + OGAP);
                        	L.push_back(matrix[1][ix][iy-1] + OGAP);
                        	L.push_back(matrix[2][ix][iy-1] + CGAP);
			}
			else
			{
				L.push_back(matrix[0][ix][iy-1]);
                                L.push_back(matrix[1][ix][iy-1]);
                                L.push_back(matrix[2][ix][iy-1]);
			}

                        matrix[2][ix][iy] = max(L,&mx);
                        patrix[2][ix][iy] = &matrix[mx][ix][iy-1];

                        L.clear();
                }
        }

        return;
}

vector<int> remFL(vector<int> exL)
{
	if((int)exL.size() < 1 + S_SKIP + E_SKIP)
		return exL;
	else
	{
		for(int is = 0; is < S_SKIP; is++)
			exL.erase(exL.begin());

		for(int ie = 0; ie < E_SKIP; ie++)
		{
			vector<int>::iterator iE = exL.end();
			iE--;
			exL.erase(iE);
		}
		return exL;
	}
}
vector<short int> remFL(vector<short int> exL)
{
        if((int)exL.size() < 1 + S_SKIP + E_SKIP)
                return exL;
        else
        {
                for(int is = 0; is < S_SKIP; is++)
                        exL.erase(exL.begin());

                for(int ie = 0; ie < E_SKIP; ie++)
                {
                        vector<short int>::iterator iE = exL.end();
                        iE--;
                        exL.erase(iE);
                }
                return exL;
        }
}
/*m_point c_check(int a, int b, double *Fr, double min_freq)
{
	m_point n = 0;

	if(a > b)
        {
        	int buf = b;
        	b = a;
        	a = buf;
        }

	if(a >= b - TOL && a <= b + TOL)
	{

		if(Fr[a] != 0)
			return (m_point)(-log(Fr[a]) + EXACT_MATCH_BONUS);
		else
			return (m_point)(-log(min_freq) + EXACT_MATCH_BONUS);
	}
		
	else if(abs(a-b) % 3 == 0 && M3) 
		for(int i = a; i <= b; i += 3)
		{
			if(i < MAX_EX)
                		n += Fr[i];
		}

	else
		for(int i = a; i <= b; i++)
		{
			if(i < MAX_EX)
				n += Fr[i];
		}	

	if(n != 0)
		n = (m_point)-log(n);
	else
		n = (m_point)-log(min_freq);

	return n;
}*/

/*m_point c_check(int a, int b, double *Fr, double min_freq, short int af, short int bf, m_point *Framefreq)
{
	m_point n = 0, multipler = 1;

	if(af != -1 && bf != -1 && af == bf)
		multipler = Framefreq[af];	

	if(a > b)
        {
        	int buf = b;
        	b = a;
        	a = buf;
        }

	if(a >= b - TOL && a <= b + TOL)
	{

		if(Fr[a] != 0)
			return (m_point)(-log(Fr[a]*multipler) + EXACT_MATCH_BONUS);
		else
			return (m_point)(-log(min_freq*multipler) + EXACT_MATCH_BONUS);
	}
		
	else
		for(int i = a; i <= b; i++)
		{
			if(i < MAX_EX)
				n += Fr[i] * multipler;
		}	

	if(n != 0)
		n = (m_point)-log(n);
	else
		n = (m_point)-log(min_freq*multipler);

	return n;
}*/

m_point c_check(int a, int b, double *Fr,short int af, short int bf, m_point *Framefreq, int frame_min, double max_diff)
{
        m_point multipler = 1;

        if(af != -1 && bf != -1 && af == bf)
                multipler = pow(Framefreq[af],2);

        //cerr << Fr[0] << "\t" << Framefreq[frame_min] << "\t" << log(Fr[0]*Framefreq[frame_min])/2 << endl;

        if(abs(a-b) < max_diff)
	{
//		cerr << "\na = " << a << " b = " << b << " return = " << (-log(Fr[abs(a-b)]*multipler) - Framefreq[3])/Framefreq[4];
//		cerr << "\nmax_diff = " << max_diff << " multipler = " << multipler << " Framefreq[3] = " << Framefreq[3] << " Framfreq[af] = " << Framefreq[af] << " Fr[abs(a-b) = " << Fr[abs(a-b)] << endl;
        //      return (-log(Fr[abs(a-b)]*multipler) - 1.67)/1.16;
                return (-log(Fr[abs(a-b)]*multipler) - Framefreq[3])/Framefreq[4];
//              return (-log(Fr[abs(a-b)]*multipler));
                //return -log(Fr[abs(a-b)]*multipler) + log(Fr[0]*Framefreq[frame_min])/2;
	}
        else
//              return -1.67/1.16;
                return -Framefreq[3]/Framefreq[4];
//              return 0;
                //return log(Fr[0])/2;
                //return log(Fr[0]*Framefreq[frame_min])/2;

}
m_point freq_c_check(int a, int b, double *Fr,short int af, short int bf, m_point *Framefreq, int frame_min, double max_diff)
{
        m_point multipler = 1;

        if(af != -1 && bf != -1 && af == bf)
                multipler = pow(Framefreq[af],2);

        if(abs(a-b) < max_diff)
                return (-log(Fr[abs(a-b)]*multipler));
        else
                return 0;

}

/*m_point max(list <m_point> L, int *mx)
{
        list<m_point>::iterator i = L.begin();
        m_point v = *i, c = 0;
	int m = 0;

        for(i = L.begin(); i !=  L.end(); i++)
        {
                if(*i > v)
                {
                        v = *i;
                        m = (int)c;
                }
                c++;
        }

        L.sort();
        i = L.end();
        i--;
        *mx = m;
	c = *i;
	L.clear();
        return c;
}*/

m_point max(vector <m_point> L, int *mx)
{
	m_point max = L[0];
	int pos = 0;
	
	for(int i = 1; i < L.size(); i++)
	{
		if(L[i] > max)
		{
			max = L[i];
			pos = i;
		}
	}

	*mx = pos;

	return max;
}

void m_pmax(m_point ***matrix,int px, int py, int *mx)
{
      //  list<m_point> L;
	vector<m_point> L;
        m_point m;

        L.push_back(matrix[0][px][py]);
        L.push_back(matrix[1][px][py]);
        L.push_back(matrix[2][px][py]);

        m =  max(L,mx);

        return;
}

vector<string> m_traceback(m_point ***matrix, m_point ****patrix, vector<int>  seq1, vector<int> seq2, int *mpos, vector<int> q_merg, vector<int> g_merg,
			   long double *vvec, int gfilenum, vector<short int> *qF, vector<short int> *gF, m_point *Framefreq)
{
	vector<string> s(4);
	ostringstream str3;

        s[0].clear();
        s[1].clear();
        s[2].clear();
	s[3].clear();

	int mp[3], mx, opos[3], cn = 0,sx,sy,gps;
	short int exact = 0, q_exact = 0, n_exact = 0, gaps = 0;
	m_point *px,score;

	for(short int i = 0; i < 3; i ++)
		mp[i] = opos[i] = mpos[i];

//	x = opos[X];
//	sy = opos[Y];
	
	m_pmax(matrix, opos[X], opos[Y], &mx);
        opos[MX] = mx;

	px = patrix[opos[MX]][opos[X]][opos[Y]];
	
	score = matrix[opos[MX]][opos[X]][opos[Y]];
	str3 << score << '\t'; 

	while(cn < 2 && (!local || matrix[opos[MX]][opos[X]][opos[Y]] != 0))
        {
		ostringstream str1,str2;
		string sbuf;
                m_crawler(matrix,px,mp);

/*			cerr << endl << opos[MX] << '\t' << opos[X] + S_SKIP + 1<< '\t' << opos[Y] + S_SKIP + 1<< endl;
			cerr << endl << mp[MX] << '\t' << mp[X] + S_SKIP + 1<< '\t' << mp[Y] + S_SKIP + 1<< endl;
			cerr << "\nOldscore = " << matrix[opos[MX]][opos[X]][opos[Y]] << endl;
			cerr << "\nscore = " << matrix[mp[MX]][mp[X]][mp[Y]] << endl << endl << endl;
*/
//		if(local && matrix[mp[MX]][mp[X]][mp[Y]] <= 0)
//		{
//			break;
//		}

                if(mp[X] != opos[X] && mp[Y] != opos[Y])
                {
				str1 << seq1[mp[X]];

				if(merge_check(q_merg,mp[X]))
					str1 << 'M';
	
				if(local)// && sx > S_SKIP)
					str1 << '[' << (mp[X]+1+S_SKIP) << ']' << ' ';
				else
					str1 << ' ';
				sbuf = str1.str();
				s[0].insert(s[0].begin(),sbuf.begin(),sbuf.end());

				str2 << seq2[mp[Y]];

				if(merge_check(g_merg,mp[Y]))
					str2 << 'M';

				if(local) //&& sy > S_SKIP)
					str2 << '[' << (mp[Y]+1+S_SKIP) << ']' << ' ';
				else
					str2 << ' ';

				sbuf = str2.str();
				s[2].insert(s[2].begin(),sbuf.begin(),sbuf.end());

				s[1].insert(s[1].begin(),' ');	

				if(seq1[mp[X]] == seq2[mp[Y]])
				{
					s[1].insert(s[1].begin(),'|');
					exact++;
				}
				else
				{
					if((seq1[mp[X]] >= seq2[mp[Y]] - SIM && seq1[mp[X]] < seq2[mp[Y]]) ||
                                            (seq1[mp[X]] <= seq2[mp[Y]] + SIM && seq1[mp[X]] > seq2[mp[Y]]))
					{	
						#ifdef DEBUG9	
						cerr << "seq1[mp[X]] = " << seq1[mp[X]] << endl;
						cerr << "seq2[mp[Y]] = " << seq2[mp[Y]] << endl;
						#endif
						q_exact++;
		//				if(abs(seq1[mp[X]] - seq2[mp[Y]]) % 3 != 0)
						if(qF->at(mp[X]+S_SKIP) != gF->at(mp[Y]+S_SKIP))
							s[1].insert(s[1].begin(),':');
						else
						{
							s[1].insert(s[1].begin(),'.');
							s[1].insert(s[1].begin(),':');
						}
					}
					else
					{
						n_exact++;
	//					if(abs(seq1[mp[X]] - seq2[mp[Y]]) % 3 != 0)
						if(qF->at(mp[X]+S_SKIP) != gF->at(mp[Y]+S_SKIP))
							s[1].insert(s[1].begin(),'*');
						else
                                                {
                                                        s[1].insert(s[1].begin(),'.');
                                                        s[1].insert(s[1].begin(),'*');
                                                }

					}
				}
				
                        }

                else if(mp[X] == opos[X] && mp[Y] != opos[Y])
                {
				str1 << "- ";
				sbuf = str1.str();
				s[0].insert(s[0].begin(),sbuf.begin(),sbuf.end());

				s[1].insert(s[1].begin(),' ');
                                s[1].insert(s[1].begin(),'G');

				str2 << seq2[mp[Y]];

				if(merge_check(g_merg,mp[Y]))
					str2 << 'M';
				
				if(local)// && sy > S_SKIP)
                                        str2 << '[' << (mp[Y]+1+S_SKIP) << ']' << ' ';
                                else
                                        str2 << ' ';

				sbuf = str2.str();
				s[2].insert(s[2].begin(),sbuf.begin(),sbuf.end());
				gaps++;
                }

                else if(mp[X] != opos[X] && mp[Y] == opos[Y])
                {
				str1 << seq1[mp[X]];

				if(merge_check(q_merg,mp[X]))
					str1 << 'M';
				
				if(local)// && sx > S_SKIP)
                                        str1 << '[' << (mp[X]+1+S_SKIP) << ']' << ' ';
                                else
                                        str1 << ' ';

				sbuf = str1.str();
				s[0].insert(s[0].begin(),sbuf.begin(),sbuf.end());

                                s[1].insert(s[1].begin(),' ');
				s[1].insert(s[1].begin(),'G');

				str2 << "- ";
				sbuf = str2.str();
				s[2].insert(s[2].begin(),sbuf.begin(),sbuf.end());
				gaps++;
                }

                px = patrix[mp[MX]][mp[X]][mp[Y]];

                if(px == patrix[0][0][0] || px == patrix[1][0][0] || px == patrix[2][0][0])
                        cn++;

                for(short int t = 0; t < 3; t++)
                        opos[t] = mp[t];
        }

	if(local || semilocal)
		vvec[0] = exp(-score) * TOT_EXONS[gfilenum] * seq1.size(); 
	else
		vvec[0] = exp(-score) * TOT_GENES[gfilenum];

/*	vvec[0] = 1 - exp(-vvec[0]);

	vvec[0] = ((score / (exact + q_exact + n_exact)) - Framefreq[3])/Framefreq[4]/sqrt(exact + q_exact + n_exact);

	vvec[0] = gsl_cdf_ugaussian_Q(vvec[0]);*/

	str3 << exact << '\t' << q_exact << '\t' << n_exact << '\t' << gaps << '\t' << vvec[0] << '\t' << vvec[1] << '\t' << vvec[2] 
	     << '\t' << iso_check(s[1]);

	s[3] = str3.str();

        return s;
}

void m_crawler(m_point ***matrix, m_point *px, int *c)
{
        for(int m = 0; m < 3; m++)
        {
                for(int ix = c[X];ix >= c[X] - 1 && ix >= 0; ix--)
                {
                        for(int iy = c[Y]; iy >= c[Y] -1 && iy >= 0; iy--)
                        {
                                if(&matrix[m][ix][iy] == px)
                                {
                                        c[X] = ix;
                                        c[Y] = iy;
                                        c[MX] = m;
                                        return;
                                }
                        }
                }
        }
}

void s_out(m_point **matrix, vector<string> s, int *mpos, gene q, gene g,m_point p_value)
{
	bool oflag = false;
	ofstream out;
	if(outfile.size() > 0)
	{
                out.open(outfile.c_str(), ios::app);
		out << endl << endl;
		oflag = true;
	}

	if(!WEB)
        	cout << endl << endl;
	else	
	{
		RES_COUNT++;
		string tabn, tabp;

		tabp = launch_blastp(q.protein_acc,queryfile,g.protein_acc,genesfiles[g.filenum], true);

		tabn = launch_blastn(q.family,queryfile,g.family,genesfiles[g.filenum], true);

		for(int d = 0; d < s.size(); d++)
		{
			if(d == 0)
			{
				for(int p = 0; p < S_SKIP; p++)
					wfile << q.exL[p] << ' ';
			}
			else if(d == 1)
			{
				for(int p = 0; p < S_SKIP; p++)
					wfile << "S ";
			}
			else if(d == 2)
			{
				for(int p = 0; p < S_SKIP; p++)
					wfile << g.exL[p] << ' ';	
			}

			wfile << s[d];

			if(d == 0)
                        {
                                for(int p = q.exL.size() - E_SKIP; p < q.exL.size(); p++)
                                        wfile << ' ' << q.exL[p] << ' ';

                        }
			else if(d == 1)
                        {
                                for(int p = 0; p < E_SKIP; p++)
                                        wfile << ' ' << "S ";
                        }
			else if(d == 2)
                        {
                                for(int p = g.exL.size() - E_SKIP; p < g.exL.size(); p++)
                                        wfile << ' ' << g.exL[p] << ' ';
                        }

			wfile << endl;
		}

		for(int i = 0; i < q.exF.size(); i++)
			wfile << q.exF[i] << '\t';

		wfile << endl;

		for(int i = 0; i < g.exF.size(); i++)
                        wfile << g.exF[i] << '\t';

		wfile << endl;
	
		wfile << q.cdsD[0] << '\t' << q.cdsD[1] << '\t' << q.cdsD[2] << '\t' << q.cdsD[3];

		wfile << endl;

		wfile << g.cdsD[0] << '\t' << g.cdsD[1] << '\t' << g.cdsD[2] << '\t' << g.cdsD[3];

		wfile << endl;
	
		for(int i = 0; i < q.exL.size(); i++)
			wfile << q.exL[i] << '\t';
		
		wfile << endl;

		for(int i = 0; i < g.exL.size(); i++)
			wfile << g.exL[i] << '\t';

		wfile << endl << tabn << endl << tabp << endl;

		for(int i = 0; i < q.ex_or_L.size(); i++)
			wfile << q.ex_or_L.at(i) << '\t';

		wfile << endl;

		for(int i = 0; i < g.ex_or_L.size(); i++)
                        wfile << g.ex_or_L.at(i) << '\t';

		wfile << endl;

		wfile << "<-" << endl;
//		wfile.close();
	}

	bool flag = true;
	short int c = 0,exact,q_exact,n_exact,gaps,tm = 0;
	m_point score;
	vector<vector <string> > vvs;
	vector<int>::iterator qf = q.exL.begin(), ql = q.exL.end(), gf = g.exL.begin(), gl = g.exL.end();

	for(short int ie = 0; ie < E_SKIP; ie++)
	{
		ql--;
		gl--;
	}

	for(short int i = 0; i < 3; i++)
		vvs.push_back(s_scan(s[i]));

	istringstream is(s[3]);

	is >> score >> exact >> q_exact >> n_exact >> gaps;
	tm = exact + q_exact + n_exact; 

	while(flag)
	{
		if (c != 0)
		{	
			if(!WEB)
				cout << endl;
			if(oflag)
				out << endl;
		}

		for(short int m = 0; m < 3; m++)	
		{
			for(short int is = 0; is < S_SKIP && c == 0 && !local && gap_count(vvs[m]) > S_SKIP + E_SKIP; is++)
			{
	 			if(m == 0)
				{
					if(!WEB)
                				cout << '(' << *qf << ')' << '\t';
					if(oflag)
						out << '(' << *qf << ')' << '\t';
					qf++;
				}

                		else if (m == 2)
				{
					if(!WEB)
                        			cout << '(' << *gf << ')' << '\t';
					if(oflag)
						out << '(' << *gf << ')' << '\t';
					gf++;
				}
			
				else
				{
					if(!WEB)
						cout << '\t';
					if(oflag)
						out << '\t';
				}
			}

			if(gap_count(vvs[m]) < S_SKIP + E_SKIP)
			{
				for(int g = 0; g < S_SKIP; g++)
				{
					if(!WEB)
						cout << '\t';
					if(oflag)
						out << '\t';
				}
			}

			for(short int i = c * EXLINE; ; i++)
			{	
				if(i == vvs[0].size())
				{	
					flag = false;
					break;
				}
				if(i - EXLINE == EXLINE *c)
					break;

				if(vvs[m][i] != "G")
				{
					if(i == c * EXLINE  && c != 0 && !local)
						for(short int sk = 0; sk < S_SKIP; sk++)
						{
							if(!WEB)
								cout << '\t';
							if(oflag)
								out << '\t';
						}
					
					if(!WEB)
						cout << vvs[m][i] << '\t';
					if(oflag)	
						out << vvs[m][i] << '\t';
				}
				else
				{
					if(i == c * EXLINE  && c != 0 && !local && m == 1)
					{
						if(!WEB)
							cout << ' ' << '\t';
                                        	if(oflag)
                                                	out << ' ' << '\t';
					}

					if(!WEB)
						cout << ' ' << '\t';
					if(oflag)
						out << ' ' << '\t';
				}
			}
			
			for(short int ie = 0; ie < E_SKIP && !flag && !local && gap_count(vvs[m]) > S_SKIP + E_SKIP; ie++)	
			{
				if(m == 0)
				{
					if(!WEB)
                        			cout << '(' << *ql << ')' << '\t';
					if(oflag)
						out << '(' << *ql << ')' << '\t';
					if(ie == E_SKIP-1)
					{
						if(!WEB)
							cout  << '('  << q.strand << ')';
						if(oflag)
							out  << '('  << q.strand << ')';
					}
					ql++;
				}

                		else if(m == 1)
				{
					if(!WEB)
                        			cout << '\t';
					if(oflag)
						out << '\t';
				}

                		else if(m == 2)
				{
					if(!WEB)
                        			cout << '(' << *gl << ')' << '\t';
					if(oflag)
						out << '(' << *gl << ')' << '\t';
					if(ie == E_SKIP-1)
					{
						if(!WEB)
							cout << '(' << g.strand << ')';
						if(oflag)
							out << '(' << g.strand << ')';
					}
					gl++;
				}
			}

		if(!WEB)
			cout << endl;
		if(oflag)
			out << endl;
		}
		c++;
	}
	
	if(!WEB)	
	{
		cout.precision(4);
		cout 	<< endl << "Exact matches:\t\t" << exact << " on " << tm << " (" << percento(exact,tm) << "%)"
 	     		<< endl << "Nearly exact matches:\t" << q_exact << " on " << tm << " (" << percento(q_exact,tm) << "%)"
	     		<< endl << "Other matches:\t\t" << n_exact << " on " << tm << " (" << percento(n_exact,tm) << "%)"
	     		<< endl; //<< "**************" << endl;
	}

	if(oflag)
	{
		out << endl << "Exact matches:\t\t" << exact << " on " << tm << " (" << percento(exact,tm) << "%)"
        	    << endl << "Almost exact matches:\t" << q_exact << " on " << tm << " (" << percento(q_exact,tm) << "%)"
		    << endl << "Other matches:\t\t" << n_exact << " on " << tm << " (" << percento(n_exact,tm) << "%)"
		    << endl; //<< "**************" << endl;
		out.close();
	}

	if(table)
	{
		ofstream tout(tabfile.c_str(),ios::app);
		tout.precision(4);

		tout << q.name << '\t' 
		     << q.exL.size()  << '\t' 
                   //  << g.name.substr(0,g.name.find("[")) << '\t' 
		     << g.name << '\t'
                     << g.exL.size() << '\t'
                     << exact << '\t' 
		     << percento(exact,tm) << '\t'
		     << q_exact << '\t' 
                     << percento(q_exact,tm) << '\t'
		     << n_exact << '\t' 
                     << percento(n_exact,tm) << '\t'
		     << tm << '\t' << gaps 
                     << '\t' << score 
		 //    << '\t' << get_pvalue(score, 1 + abs((int)g.exL.size() - (int)q.exL.size()))
		     << '\t' << p_value
                     << endl; 

		tout.close();
	}
 

	vvs.clear();
        return;

}

vector<string> s_scan(string s)
{
	string buf;
	vector<string> vs;	
	istringstream str(s);

	do
	{
		str >> buf;
		vs.push_back(buf);
	}
	while(str);

	vs.erase(vs.end());

	return vs;
}
	
vector<int> best_pos(vector<m_point> score, int best)
{
	vector<int> bpos;

	if(best > score.size())
		best = score.size();

	for(short int b = 0; b < best; b++)
	{
		int p = 0;
		m_point bscore = -INF;

		for(int i = 0; i < score.size(); i++)
		{
			if(score[i] > bscore)
			{
				bscore = score[i];
				p = i;
			}
		}
		
		bpos.push_back(p);
		score[p] = -INF;
	}

	return bpos;
}
	
void clineparser(int argc, char **argv)
{
	if(argc == 1)
		displayhelp();

	for(int i = 1; i < argc; i++)
	{
		string buf = argv[i];

		if(buf == "-q")
		{
			if(i < argc)
				queryfile = argv[++i];

			cerr << "\nqueryfile = " << queryfile << endl;

			continue;
		}
		else if(buf == "-O")
		{
			genesfiles.clear();

			do
			{
				i++;
				buf = argv[i];
				if(i < argc && !useall)
				{
					if(buf == "HS")
						genesfiles.push_back(Hs);
					else if(buf == "MM")
						genesfiles.push_back(Mm);
					else if(buf == "RN")
						genesfiles.push_back(Rn);
					else if(buf == "DR")
						genesfiles.push_back(Dr);
					else if(buf == "DM")
						genesfiles.push_back(Dm);
					else if(buf == "XT")
						genesfiles.push_back(Xt);
					else if(buf == "GG")
						genesfiles.push_back(Gg);
					else if(buf == "FR")
						genesfiles.push_back(Fr);
					else if(buf == "CI")
						genesfiles.push_back(Ci);
					else if(buf == "CE")
						genesfiles.push_back(Ce);
					else if(buf == "PT")
                                                genesfiles.push_back(Pt);
					else if(buf == "SP")
						genesfiles.push_back(Sp);
					else if(buf == "BT")
						genesfiles.push_back(Bt);
					else if(buf == "MD")
						genesfiles.push_back(Md);
					else if(buf == "TN")
						genesfiles.push_back(Tn);
					else if(buf == "AM")
						genesfiles.push_back(Am);
					else if(buf[0] != '-')
                                                genesfiles.push_back(buf);
				}
				else if(i < argc && useall)
				{
					if(buf == "HS")
                                                genesfiles.push_back(HsA);
                                        else if(buf == "MM")
                                                genesfiles.push_back(MmA);
                                        else if(buf == "RN")
                                                genesfiles.push_back(RnA);
                                        else if(buf == "DR")
                                                genesfiles.push_back(DrA);
                                        else if(buf == "DM")
                                                genesfiles.push_back(DmA);
                                        else if(buf == "XT")
                                                genesfiles.push_back(XtA);
                                        else if(buf == "GG")
                                                genesfiles.push_back(GgA);
                                        else if(buf == "FR")
                                                genesfiles.push_back(FrA);
                                        else if(buf == "CI")
                                                genesfiles.push_back(CiA);
                                        else if(buf == "CE")
                                                genesfiles.push_back(CeA);
                                        else if(buf == "PT")
                                                genesfiles.push_back(PtA);
                                        else if(buf == "SP")
                                                genesfiles.push_back(SpA);
                                        else if(buf == "BT")
                                                genesfiles.push_back(BtA);
                                        else if(buf == "MD")
                                                genesfiles.push_back(MdA);
                                        else if(buf == "TN")
                                                genesfiles.push_back(TnA);
                                        else if(buf == "AM")
                                                genesfiles.push_back(AmA);
					else if(buf[0] != '-')
						genesfiles.push_back(buf);
				}
			} while(buf[0] != '-' && i < argc - 1);

			if(buf[0] == '-')
				i--;
			
			continue;
		}

		else if(buf == "-l")
		{
			local = true;
			semilocal = false;
			continue;
		}

		else if(buf == "-gl")
		{
			semilocal = true;
			local = false;
			continue;
		}

		else if(buf == "-og" && i < argc - 1)
		{
			i++;
			buf = argv[i];
			istringstream str(buf);
			str >> OGAP;
			dynamic_gap = false;
			continue;
		}

		
		else if(buf == "-cg" && i < argc - 1)
		{
			i++;
			buf = argv[i];
			istringstream str(buf);
			str >> CGAP;
			dynamic_gap = false;
			continue;
		}

		
/*		else if(buf == "-si" && i < argc - 1)
		{
			i++;
			buf = argv[i];
			istringstream str(buf);
			str >> SIM;
			continue;
		}*/

		else if(buf == "-el" && i < argc - 1)
		{
			i++;
			buf = argv[i];
			istringstream str(buf);
			str >> EXLINE;
			continue;
		}

		else if(buf == "-b" && i < argc - 1)
		{
			i++;
			buf = argv[i];
			istringstream str(buf);
			str >> V_BEST;
			continue;
		}

		else if(buf == "-B" && i < argc - 1)
		{
			i++;
			buf = argv[i];
			istringstream str(buf);
			str >> BEST;
			continue;
		}
		
		else if(buf == "-ss" && i < argc - 1)
		{
			i++;
			buf = argv[i];
			istringstream str(buf);
			str >> S_SKIP;
			continue;
		}

		else if(buf == "-es" && i < argc - 1)
		{
			i++;
			buf = argv[i];
			istringstream str(buf);
			str >> E_SKIP;
			continue;
		}

/*		else if(buf == "-mc" && i < argc - 1)
		{
			i++;
                        buf = argv[i];
                        istringstream str(buf);
                        str >> MAX_CYCLE;
                        continue;
		}*/

/*		else if(buf == "-nws" && i < argc - 1)
                {
                        i++;
                        buf = argv[i];
                        istringstream str(buf);
                        str >> NWS;
			if(NWS < 0)
				NWS = 0;
                        continue;
                }

		else if(buf == "-pws" && i < argc - 1)
                {
                        i++;
                        buf = argv[i];
                        istringstream str(buf);
                        str >> PWS;
			if(PWS < 0)
				NWS = 0;
                        continue;
                }

		else if(buf == "-mb" && i < argc - 1)
		{
			i++;
			buf = argv[i];
			istringstream str(buf);
			str >> EXACT_MATCH_BONUS;
			continue;
		}*/

		else if(buf == "-gs" && i < argc - 1)
		{
			 if(i < argc)
                                gsel = argv[++i];
			gensel = true;
                        continue;
		}

		
		else if(buf == "-qs" && i < argc - 1)
		{
			 if(i < argc)
                                qsel = argv[++i];
			querysel = true;
			cerr << "\nquerygene = " << qsel << endl;
                        continue;
		}

/*		else if(buf == "-qm")
		{
			qmerge = true;
			continue;
		}*/

		else if(buf == "-ta")
		{
			tabapp = true;
			continue;
		}

/*		else if(buf == "-gm")
                {
                        gmerge = true;
                        continue;
                }*/

		else if(buf == "-t")
                {
                        table = true;
                        continue;
                }

/*		else if(buf == "-noM3")
		{
			M3 = false;
			continue;
		}

		else if(buf == "-boh")
                {
                        blasthitonly = true;
                        continue;
                }

		else if(buf == "-all")
                {
                        useall = true;
                        continue;
                }*/

		else if(buf == "-freq")
                {
                       	freq = true;
                        continue;
                }
/*		else if(buf == "-web")
		{
			WEB = true;
			continue;
		}
		else if(buf == "-score_out")
                {
                        score_out = true;
                        continue;
                }
		else if(buf == "-qid" && i < argc - 1)
		{
			 if(i < argc)
                                QID = argv[++i];

                        continue;
		}*/

		else if(buf == "-h")
			displayhelp();
		
		else
		{
			cerr << "\nBad argument: " << buf << endl << endl;
			exit(EXIT_FAILURE);
		}
	}

	argcheck();
	return;
}

void displayhelp()
{
        cerr << "\n\nEXALIGN quick help. Please read documentation for further help.\n\nUse:\n\texalign -q query_file -O organism_file [OPTIONS]\n\n"
             << "OPTIONS:\n\n";

        cerr << "\t-q ...\tThe file containing the query gene structure(s).\n"
             << "\t-O\tThe organism file containing the gene structure database to search against.\n"
             << "\t-l\tToggle local alignment instead of global alignment.\n"
             << "\t-gl\tToggle glocal alignment instead of global alignment.\n"
             << "\t-t\tA file containing a summary table of the results is written.\n"
             << "\t-qs ...\tOnly the given gene structure of the query will be aligned (it must be contained in the query file).\n"
             << "\t-gs ...\tAlign the query only with the given gene structure (it must be contained in the organism file).\n"
//           << "\t-mc num\tSet the max number of merging cycles (effects on merging process speed and deepness). (default: " << MAX_CYCLE << ")\n"
             << "\t-og num\tSet the OPENING GAP penalty to num. (default: \"dynamic gap penalty\")" << endl
             << "\t-cg num\tSet the CONTINUING GAP penalty to num. (default: \"dynamic gap penalty \")" << endl
             << "\t-el num\tSet the number of exons per line displayed in the output.\n\t\t"
             << "Try to lower this value if your output is messed up. (default: " << EXLINE << ')' << endl
             << "\t-b num\tChoose the maximum number of alignments you want in your output.\n\t\t"
             << "Alignments are sorted by descending score. (default: " << V_BEST << ')' << endl
             << "\t-B num\tHow many of the best scoring alignments will proceed to the merging exons step. (default: " << BEST << ")\n\t\t"
             << "WARNING: increasing this will increase computational time.\n"
             << "\t-ss num\tSet the number of start exons that are not aligned\n\t\t"
             << "(default: " << S_SKIP << ')' << endl
             << "\t-es num\tSet the number of end exons that are not aligned\n\t\t"
             <<"(default: " << E_SKIP << ')' << endl
             << "\t-freq\tThe program generates the new frequency file for the specified Organism file.\n"
             << "\t-h\tDisplay this help." << endl << endl;

             exit(EXIT_SUCCESS);
}
			
void argcheck()
{
	if(queryfile.size() == 0 && !freq)
	{
		cerr << "\n\nYou must specify a query file\n\n";
		exit(EXIT_FAILURE);
	}
	else if(queryfile.size() == 0 && freq)
		query = false;

	if(BEST < 1)
	{
		cerr << "\n\n-B value must be greater or equal than 1.\n\n";
		exit(EXIT_FAILURE);
	}
	if(V_BEST < 1)
	{
		cerr << "\n\n-b value must be greater than 1.\n\n";
		exit(EXIT_FAILURE);
	}
	if(WEB)
	{
		BEST = 2 * V_BEST;
	}
	if(CGAP < OGAP || CGAP > 0)
	{
		cerr << "\n\nCONTINUING GAP must be smaller or equal than OPENING GAP and they must be smaller than 0.\n\n";
		exit(EXIT_FAILURE);
	}
	if(CGAP == 0 || OGAP == 0)
	{
		dynamic_gap = true;
	}
	if(FILTER < 0)
	{
		cerr  << "\n\nPRE-FILTER value must be greater than 0.\n\n";
		exit(EXIT_FAILURE);
	}
	if(EXLINE < 5)
	{
		cerr << "\n\n-el value must be greater than 5.\n\n";
		exit(EXIT_FAILURE);
	}
	if(SIM < 1)
	{
		cerr << "\n\n-si value must be > 0.\n\n";
		exit(EXIT_FAILURE);
	}
	if(MAX_CYCLE < 0)
	{
		cerr << "\n\n-mc value must be >= 0.\n\n";
		exit(EXIT_FAILURE);
	}

	for(int gf = 0; gf < genesfiles.size(); gf++)
	{
		freqfiles.push_back(genesfiles[gf]);

		if(freqfiles[freqfiles.size() - 1].find(".kg") == string::npos)
        		freqfiles[freqfiles.size() - 1] += ".freq";
		else
		{
			freqfiles[freqfiles.size() - 1] = freqfiles[freqfiles.size() - 1].substr(0,freqfiles[freqfiles.size() - 1].find(".kg"));
			freqfiles[freqfiles.size() - 1] += ".rf.freq";
		}
	}

	if(querysel && gensel)
	{
		freqfiles.push_back(queryfile);

		if(freqfiles[freqfiles.size() - 1].find(".kg") == string::npos)
			freqfiles[freqfiles.size() - 1] += ".freq";
		else
                {
                        freqfiles[freqfiles.size() - 1] = freqfiles[freqfiles.size() - 1].substr(0,freqfiles[freqfiles.size() - 1].find(".kg"));
                        freqfiles[freqfiles.size() - 1] += ".rf.freq";
                }

		V_BEST = 1;
		BEST = 1;
	}
}

/*m_point lmax(list <m_point> L, int *mx)
{
        list<m_point> L2 = L;
        list<m_point>::iterator i = L.begin(), i2;
        m_point v = *i;
	int m = 0, c = 0;

        L2.sort();
        i2 = L2.end();
        i2--;

        if(*i2 == 0)
        {
                *mx = 0;
                return 0;
        }
        else
                L.erase(L.begin());

        for(i = L.begin(); i !=  L.end(); i++)
        {
                if(*i > v)
                {
                        v = *i;
                        m = c;
                }
                c++;
        }

        L.sort();
        i = L.end();
        i--;
        *mx = m;
        return *i;
}*/
m_point lmax(vector <m_point> L, int *mx)
{
	m_point max = L[0];
	int pos = 0;	

	for(int i = 1; i < L.size(); i++)
	{
		if(L[i] > max)
		{
			max = L[i];
			pos = i;
		}
	}

	if(max == 0)
	{
		*mx = 0;
		return 0;
	}
	else
	{
		*mx = pos - 1;
		return max;
	}
}
	
                                            
m_point m_max(m_point ***matrix, int *mpos)
{
        m_point value = -INF;
        for(int m = 0; m < 3; m++)
        {
                for(int ix = 1; ix < x; ix++)
                {
                        for(int iy = 1; iy < y; iy++)
                        {
                                if(matrix[m][ix][iy] > value)
                                {
                                        mpos[MX] = m;
                                        mpos[X] = ix;
                                        mpos[Y] = iy;
                                        value = matrix[m][ix][iy];
                                }
                        }
                }
        }

        return value;
}

int gap_count(vector<string> s)
{
	int counter = 0;
	for(int i = 0; i < s.size(); i++)
	{
		if(s[i] != "-")
			counter++;
	}
	return counter;
}

/*vector<gene> l_merged(vector<gene> genes)
{
	vector<gene> merged;

	for(int i1 = 0; i1 < genes.size(); i1++)
	{
		merged.push_back(genes[i1]);

		for(int i2 = S_SKIP; i2 < genes[i1].exL.size() - 1 - E_SKIP; i2++)
		{
			ostringstream str1;
			str1 << ' ' << i2 + 1 << '(' << genes[i1].exL[i2] << ')' << " + " << i2 + 2 << '(' << genes[i1].exL[i2+1] << ')';
			gene g;
			g.merged = true;
			g.merg_pos.push_back(i2 - S_SKIP);
			g.chr = genes[i1].chr;
			g.strand = genes[i1].strand;
			g.name = genes[i1].name;
			g.filenum = genes[i1].filenum;
			g.name += str1.str(); 

			for(int i3 = 0; i3 < genes[i1].exL.size(); i3++)
			{
				if(i3 == i2)
					g.exL.push_back(genes[i1].exL[i2]+genes[i1].exL[i2+1]);

				else if(i3 == i2 + 1)
					continue;
				else 
					g.exL.push_back(genes[i1].exL[i3]);
			}

			g.exN = g.exL.size();

			if(g.exN > 0)
				merged.push_back(g);
		}			
	}
		
	return merged;
}*/

void frame_freq(vector<gene> *genes, m_point *framefreq)
{
	int i = 0, count = 0;

	for(vector<gene>::iterator vgi = genes->begin(); vgi < genes->end(); vgi++)
	{
		for(int a = 0; a < genes->at(i).exF.size(); a++)
		{	
			if(genes->at(i).exF[a] >= 0)
			{
				count++;
				if(genes->at(i).exF[a] == 0)
					framefreq[0]++;
				if(genes->at(i).exF[a] == 1)	
					framefreq[1]++;
				if(genes->at(i).exF[a] == 2)
					framefreq[2]++;
			}
		}

		i++;
	}

	framefreq[0] /= (m_point)count;
	framefreq[1] /= (m_point)count;
	framefreq[2] /= (m_point)count;
//	cerr << "\ncount = "<< count << endl;
}
	
void freqgen(vector<gene>::iterator iS, vector<gene>::iterator iE, int filenum)
{
	vector<gene> genes(iS,iE);
	list<int> exL;
	int nex, pvsize = 0, count = 0, gcount = 0, en20; 
        m_point framefreq[5] = {0,0,0,0,0};
	list<int>::iterator i;
	
	cerr << "Generating frequencies file:" << endl;

	frame_freq(&genes, framefreq);

	ofstream out(freqfile.c_str());

	out.precision(20);

	for(int i1 = 0; i1 < genes.size(); i1++)
	{
		if(genes[i1].exN > S_SKIP + E_SKIP && genes[i1].filenum == filenum)
			for(int i2 = S_SKIP; i2 < int(genes[i1].exL.size() - E_SKIP); i2++)
				exL.push_back(genes[i1].exL[i2]);
		else if(genes[i1].exN <= S_SKIP + E_SKIP && genes[i1].filenum == filenum)
			for(int i2 = 0; i2 < genes[i1].exN; i2++)
				exL.push_back(genes[i1].exL[i2]);
	}

	cerr << "Total exons: " << exL.size() << " ss = " << S_SKIP << " es = " << E_SKIP << endl;
	out  << exL.size() << '\t' << genes.size() << " ss = " << S_SKIP << " es = " << E_SKIP << '\n' << framefreq[0] << '\t' << framefreq[1] 
            << '\t' << framefreq[2] << endl; 
	nex = exL.size();
	en20 = nex/20;

	if (en20 == 0)
		en20++;

	cerr << "Sorting...";
	exL.sort();
	cerr << "done" << endl;

	exL.push_back(0);

	#ifdef CFREQ
	int prevf = pvsize;
	#endif
	
	cerr << "Counting";
	for(i = exL.begin(); i != exL.end(); i++)
	{

		count++;
		if(gcount%en20 == 0)
			cerr << '.';
		gcount++;

		if(*i  == pvsize)
			continue;
		else 
		{
			#ifdef CFREQ
			if(pvsize - prevf > 1)
			{
				int a = pvsize - prevf;
				int b = 1;
				while(a > 1)
				{
					out << prevf + b << '\t' << 0 << '\t' << 0 << endl;
					a--;
					b++;
				}
			}
			prevf = pvsize;
			#endif

			if(pvsize)
				out << pvsize << '\t' << double((double)count/(double)nex) << '\t' << count << endl;
	
			pvsize = *i;
			count = 0;
		}
	}

	cerr << endl;
	out.close();
	cerr << "Frequencies written in " << freqfile << endl;

	return;
}

map<int,double> getfreq(m_point *framefreq, bool write_new)
{
	ifstream in(freqfile.c_str());

	map<int,double> Fr;
	string line, trash;
	int tot_exons,tot_genes;
	vector<int> lenghts;
	
	if(!in)
	{
//		freq = true;
//		return Fr;
		cerr << "\nFrequency file not found: " << freqfile << endl << endl;
		exit(EXIT_FAILURE);
		
	}

	getline(in,line);

	istringstream str0(line),str2;
	str0 >> tot_exons >> tot_genes;	

	TOT_EXONS.push_back(tot_exons);
	TOT_GENES.push_back(tot_genes);

	getline(in,line);

	str2.str(line);

	str2 >> framefreq[0] >> framefreq[1] >> framefreq[2] >> framefreq[3] >> framefreq[4];

	while(getline(in,line))
	{
		istringstream str1(line);
		int L;
		double F;

		str1 >> L >> F;

		Fr[L] = F;
		lenghts.push_back(L);
	}

	in.close();

	if(write_new)
		Fr = write_new_freqfile(Fr, &lenghts, tot_exons, tot_genes, framefreq);	

	return Fr;
}

map<int,double> write_new_freqfile(map<int,double> Fr, vector<int> *L, int tot_exons, int tot_genes, m_point *framefreq)
{
	map<int, map<int,double> > GSS;	
	map<int, double> N_FR_map;
	vector<double> diffcount, diff_cum; 
	double *FrV, *N_FrV, min_freq;
	double average = 0, var = 0;
	int m_size = (L->at(L->size() - 1)), fmf;

	framefreq[3] = 0;

	fmf = frame_min_freq(framefreq); 

	FrV = map_to_array(Fr,&min_freq);

	for(int i = 0; i < m_size; i++)
		diffcount.push_back(0);

	for(int i = 0; i < Fr.size(); i++)
		for(int h = i; h < Fr.size(); h++)
		{
			if(L->at(h) - L->at(i))
				diffcount[L->at(h) - L->at(i)] += (FrV[L->at(i)] * FrV[L->at(h)] * tot_exons * tot_exons);
			else
				diffcount[0] += ((FrV[L->at(i)] * tot_exons * ((FrV[L->at(i)] * tot_exons) - 1))/2);
		}

//	for(int i = 0; i < diffcount.size(); i++)
//		cout << i << '\t' << diffcount[i] << endl;

	diff_cum.push_back(diffcount[0]);

	for(int i = 1; i < diffcount.size(); i++)
		diff_cum.push_back(diff_cum[i - 1] + diffcount[i]);

//	for(int i = 0; i < diff_cum.size(); i++)
//		cout << i << '\t' << diff_cum[i] << endl;

	N_FrV = new double[diff_cum.size()];

	for(int i = 0; i < diff_cum.size(); i++)
		N_FrV[i] = (double)diff_cum[i]/(double)diff_cum[diff_cum.size() - 1];

	double frame_square[3];

        frame_square[0] = pow(framefreq[0],2);
        frame_square[1] = pow(framefreq[1],2);
        frame_square[2] = pow(framefreq[2],2);

	for(int i = 0; i < L->size(); i++)
        {
                for(int h = 0; h < L->size(); h++)
                {
                        GSS[i][h] = freq_c_check(L->at(i),L->at(h), N_FrV, -1, -1, framefreq, fmf, diffcount.size());

			average += FrV[i] * FrV[h] * frame_square[0] * (GSS[i][h] - log(frame_square[0]));
                        average += FrV[i] * FrV[h] * frame_square[1] * (GSS[i][h] - log(frame_square[1]));
                        average += FrV[i] * FrV[h] * frame_square[2] * (GSS[i][h] - log(frame_square[2]));
                        average += FrV[i] * FrV[h] * (1 - frame_square[0] - frame_square[1] - frame_square[2]) * GSS[i][h];
                }
        }	
	
	for(int i = 0; i < L->size(); i++)
                for(int h = 0; h < L->size(); h++)
                        var += FrV[i] * FrV[h] * pow((GSS[i][h] - average),2);

	framefreq[3] = average;
//	framefreq[3] = -log(N_FrV[0]*framefreq[fmf])/2;
	framefreq[4] = sqrt(var);

	ofstream out(freqfile.c_str(),ios::out);

	out.precision(10);

	out << tot_exons << '\t' << tot_genes << endl;
	out << framefreq[0] << '\t' << framefreq[1] << '\t' << framefreq[2] << '\t' << framefreq[3] << '\t' << framefreq[4] << endl;
	
	for(int i = 0; i < diff_cum.size(); i++)
	{
		out << i << '\t' << N_FrV[i] << endl;
		N_FR_map[i] = N_FrV[i];
	}

	out.close();

	cerr << endl << "Average score = " << average << endl;
	cerr << endl << "Dev. Std. = " << framefreq[4] << endl;

	delete[] N_FrV;

	return N_FR_map;
}
							

vector<gene> fus_match(vector<gene> genes, vector<match> matches, int q)
{
	vector<gene> fgenes;

	for(int i = 0; i < matches.size(); i++)
		if(matches[i].pos[1] == q)
			fgenes.push_back(genes[matches[i].pos[0]]);

	return fgenes;
}

vector<gene> gene_fusion(gene q, vector<gene> genes, int sc, bool fflag)
{
	vector<gene> fgenes;

        if(fflag)
                fgenes = genes;

        for(int i1 = sc; i1 < genes.size(); i1++)
        {
                vector<int> s_ex;

                for(int i2 = S_SKIP; i2 < (int)genes[i1].exL.size() - E_SKIP - 1; i2++)
                {
                        for(int i3 = S_SKIP; i3 < (int)q.exL.size() - E_SKIP; i3++)
                        {
                                if(genes[i1].exL[i2]+genes[i1].exL[i2+1] == q.exL[i3]) //&& (!genes[i1].m_flag[i2] || !genes[i1].m_flag[i2+1])) 
                                        s_ex.push_back(i2);
                        }
                }
	
	    	if((int)genes[i1].exL.size() > S_SKIP + E_SKIP + 1)
		{
			for(int g0 = 0; g0 < s_ex.size(); g0++)
			{
				int ex = s_ex[g0];

				gene g;
                        	ostringstream str1;

                       /* 	str1 << ' ' << '[' << rex(ex) << '(' << genes[i1].exL[ex] << ')'
                             	     << '+' << rex(ex) + 1 << '(' << genes[i1].exL[ex+1] << ')' << ']'
                                     << ' ' << '=' << ' ' << genes[i1].exL[ex]+genes[i1].exL[ex+1] ;*/

				str1 << ' ' << rex(ex) << '(' << genes[i1].exL[ex] << ')'
                                     << '+' << rex(ex) + 1 << '(' << genes[i1].exL[ex+1] << ')'; 

                        	g = genes[i1];

                        	g.exL.clear();
				g.exF.clear();

                       // 	g.name += str1.str();

                        	g.merged = true;

				g.filenum = genes[i1].filenum;
				g.family = genes[i1].family;
				g.protein_acc = genes[i1].protein_acc;
				g.merg_pos = genes[i1].merg_pos;
				g.merged_with.push_back(q.name);

				if(WEB)
					g.ex_or_L = genes[i1].ex_or_L; //new

                        	for(int h0 = 0; h0 < S_SKIP; h0++)
				{
                                	g.exL.push_back(genes[i1].exL[h0]);
					g.exF.push_back(genes[i1].exF[h0]);
				}

                        	for(int i4 = S_SKIP; i4 < (int)genes[i1].exL.size() - E_SKIP; i4++)
                        	{
                                	if(i4 != ex)
					{
                                        	g.exL.push_back(genes[i1].exL[i4]);
						g.exF.push_back(genes[i1].exF[i4]);
					}
                                	else
                                	{
						g.exF.push_back(genes[i1].exF[i4]);
                                        	g.exL.push_back(genes[i1].exL[i4] + genes[i1].exL[++i4]);

						if(g.cdsD[0] > i4)
							g.cdsD[0]--;
						if(g.cdsD[2] > i4)
							g.cdsD[2]--;	

						if(g.merg_pos.size() == 0)
						{
                                        		g.merg_pos.push_back(g.exL.size() - 1 - S_SKIP);
							g.name += str1.str();
						}
							
						else
						{
							for(int bb = 0; bb < g.merg_pos.size(); bb++)
								if(g.merg_pos[bb] >= g.exL.size() - 1 - S_SKIP)
								{
									istringstream ibuf(g.name);
									ostringstream buf1;	
									string sbuf;
									
									ibuf >> sbuf;

									buf1 << sbuf << ' ' << str1.str() << '-';

									ibuf >> sbuf;

									buf1 << sbuf;

									g.name = buf1.str();

									for(int gg = 0; gg < g.merg_pos.size(); gg++)
									{
										if((g.merg_pos[gg] >= g.exL.size() - 1 - S_SKIP))
											g.merg_pos[gg]--;
									}
	
									break;
								}
								else         //new
                        					{
                                					g.name += "-";
                                					g.name += str1.str().substr(1);
									break;
                        					}
									
							g.merg_pos.push_back(g.exL.size() - 1 - S_SKIP);
						}
                                	}
                        	}

                        	for(int h1 = 0; h1 < E_SKIP; h1++)
				{
                                	g.exL.push_back(genes[i1].exL[genes[i1].exL.size() - 1 - h1]);
					g.exF.push_back(genes[i1].exF[genes[i1].exF.size() - 1 - h1]);
				}

                                g.exN = g.exL.size();

				fgenes.push_back(g);
			}
                }
        }

	return fgenes;
}	

vector<gene> gene_tri_fusion(gene q, vector<gene> genes, int sc, bool fflag)
{
	vector<gene> fgenes;

//	cerr << "sc = " << sc << endl ;
//	cerr << "genes.size = " << genes.size() << endl; 

        if(fflag)
                fgenes = genes;

//	cerr << "genes.size = " << genes.size() << endl;

        for(int i1 = sc; i1 < genes.size(); i1++)
        {
                vector<int> s_ex;

                for(int i2 = S_SKIP; i2 < (int)genes[i1].exL.size() - E_SKIP - 2; i2++)
                {
                        for(int i3 = S_SKIP; i3 < (int)q.exL.size() - E_SKIP; i3++)
                        {
//				cerr << (genes[i1].exL[i2]+genes[i1].exL[i2+1]+genes[i1].exL[i2+2]) << '\t' << q.exL[i3] << endl;
                                if(genes[i1].exL[i2]+genes[i1].exL[i2+1]+genes[i1].exL[i2+2] == q.exL[i3]) //&& (!genes[i1].m_flag[i2] || !genes[i1].m_flag[i2+1]))
				{
                                        s_ex.push_back(i2);
				}
                        }
                }
	
	    	if((int)genes[i1].exL.size() > S_SKIP + E_SKIP + 2)
		{
			for(int g0 = 0; g0 < s_ex.size(); g0++)
			{
				int ex = s_ex[g0];

				gene g;
                        	ostringstream str1;

                       /* 	str1 << ' ' << '[' << rex(ex) << '(' << genes[i1].exL[ex] << ')'
                             	     << '+' << rex(ex) + 1 << '(' << genes[i1].exL[ex+1] << ')' << ']'
                                     << ' ' << '=' << ' ' << genes[i1].exL[ex]+genes[i1].exL[ex+1] ;*/

				str1 << ' ' << rex(ex) << '(' << genes[i1].exL[ex] << ')'
                                     << '+' << rex(ex) + 1 << '(' << genes[i1].exL[ex+1] << ')'
                                     << '+' << rex(ex) + 2 << '(' << genes[i1].exL[ex+2] << ')'; 

                        	g = genes[i1];

                        	g.exL.clear();
				g.exF.clear();

                       // 	g.name += str1.str();

				g.tri_merged = true;

				g.filenum = genes[i1].filenum;
				g.family = genes[i1].family;
				g.protein_acc = genes[i1].protein_acc;
				g.merg_pos = genes[i1].merg_pos;
				g.merged_with.push_back(q.name);

				if(WEB)
					g.ex_or_L = genes[i1].ex_or_L; //new

                        	for(int h0 = 0; h0 < S_SKIP; h0++)
				{
                                	g.exL.push_back(genes[i1].exL[h0]);
					g.exF.push_back(genes[i1].exF[h0]);
				}

                        	for(int i4 = S_SKIP; i4 < (int)genes[i1].exL.size() - E_SKIP; i4++)
                        	{
                                	if(i4 != ex)
					{
                                        	g.exL.push_back(genes[i1].exL[i4]);
						g.exF.push_back(genes[i1].exF[i4]);
					}
                                	else
                                	{
						g.exF.push_back(genes[i1].exF[i4]);
                                        	g.exL.push_back(genes[i1].exL[i4] + genes[i1].exL[++i4] + genes[i1].exL[++i4]);

						if(g.cdsD[0] > i4)
						{
							g.cdsD[0]--;
							g.cdsD[0]--;
						}
						if(g.cdsD[2] > i4)
						{
							g.cdsD[2]--;	
							g.cdsD[2]--;
						}

						if(g.merg_pos.size() == 0)
						{
                                       			g.merg_pos.push_back(g.exL.size() - 1 - S_SKIP);
							g.name += str1.str();
						}
						else
						{
                                                        for(int bb = 0; bb < g.merg_pos.size(); bb++)
                                                                if(g.merg_pos[bb] >= g.exL.size() - 1 - S_SKIP)
                                                                {
                                                                        istringstream ibuf(g.name);
                                                                        ostringstream buf1;
                                                                        string sbuf;

                                                                        ibuf >> sbuf;

                                                                        buf1 << sbuf << ' ' << str1.str() << '-';

                                                                        ibuf >> sbuf;

                                                                        buf1 << sbuf;

                                                                        g.name = buf1.str();

                                                                        for(int gg = 0; gg < g.merg_pos.size(); gg++)
                                                                        {
                                                                                if((g.merg_pos[gg] >= g.exL.size() - 1 - S_SKIP))
										{
                                                                                        g.merg_pos[gg]--;
											g.merg_pos[gg]--;
										}
                                                                        }

                                                                        break;
                                                                }
								else         //new
                                                                {
                                                                        g.name += "-";
                                                                        g.name += str1.str().substr(1);
									break;
                                                                }

                                                        g.merg_pos.push_back(g.exL.size() - 1 - S_SKIP);
                                                }
                                	}
                        	}

                        	for(int h1 = 0; h1 < E_SKIP; h1++)
				{
                                	g.exL.push_back(genes[i1].exL[genes[i1].exL.size() - 1 - h1]);
					g.exF.push_back(genes[i1].exF[genes[i1].exF.size() - 1 - h1]);
				}

                                g.exN = g.exL.size();

				fgenes.push_back(g);

//				cerr << g.name << '\t';

				for(int ww = 0; ww < g.merg_pos.size(); ww++)
					cerr << g.merg_pos[ww] << '\t';
				cerr << endl;
			}
                }
        }

//	cerr << "\nfgenes.size = " << fgenes.size() << endl;
	return fgenes;
}	
/*vector<gene> gene_fusion(gene q, vector<gene> genes, double *FrV, int sc, bool fflag)
{
	vector<gene> fgenes;

	if(fflag)
		fgenes = genes;
	
	for(int i1 = sc; i1 < genes.size(); i1++)	
	{
		vector<int> s_ex;
		vector<double> sf_ex;

		for(int i2 = S_SKIP; i2 < (int)genes[i1].exL.size() - E_SKIP - 1; i2++)
		{
			for(int i3 = S_SKIP; i3 < (int)q.exL.size() - E_SKIP; i3++)
			{
				if(genes[i1].exL[i2]+genes[i1].exL[i2+1] == q.exL[i3])
				{
					s_ex.push_back(i2);
					sf_ex.push_back(FrV[genes[i1].exL[i2] + genes[i1].exL[i2+1]]);
				}
			}
		}

		if(s_ex.size() > 0 && (int)genes[i1].exL.size() > S_SKIP + E_SKIP + 1)
		{
			int ex = s_ex[get_min_freq(sf_ex)];
			gene g;
                        ostringstream str1;

                        str1 << ' ' << '[' << rex(ex) << '(' << genes[i1].exL[ex] << ')'
                             << '+' << rex(ex) + 1 << '(' << genes[i1].exL[ex+1] << ')' << ']'
                             << ' ' << '=' << ' ' << genes[i1].exL[ex]+genes[i1].exL[ex+1] ;

			g = genes[i1];

                        g.exL.clear();

                        g.name += str1.str();

                        g.merged = true;
			
			for(int h0 = 0; h0 < S_SKIP; h0++)
				g.exL.push_back(genes[i1].exL[h0]);

                        for(int i4 = S_SKIP; i4 < (int)genes[i1].exL.size() - E_SKIP; i4++)
			{
				if(i4 != ex)
					g.exL.push_back(genes[i1].exL[i4]);
				else
				{
					g.exL.push_back(genes[i1].exL[i4] + genes[i1].exL[++i4]);	
					g.merg_pos.push_back(g.exL.size() - 1 - S_SKIP);
				}
			}

			for(int h1 = 0; h1 < E_SKIP; h1++)
				g.exL.push_back(genes[i1].exL[genes[i1].exL.size() - 1  - h1]);

			g.exN = g.exL.size();

			fgenes.push_back(g);
		}
	}
		
	return fgenes;
}*/
					
vector<gene> gene_merge(vector<gene> genes, vector<gene> fusion)
{
	for(int i = 0; i < (int)fusion.size(); i++)
		genes.push_back(fusion[i]);

	return genes;
}

vector<int> match_merge(vector<int> matches, int gn, int gfn)
{
	for(int i = gn; i < gfn; i++)
		matches.push_back(i);

	return matches;
}

double* map_to_array(map<int,double> Fr, double *min_freq)
{
	vector<double> fv;
	double *FrV;	

	for(int fs = 0; fs < Fr.size(); fs++)
                fv.push_back(Fr[fs]);

        FrV = new double[fv.size()];

	MAX_EX = fv.size();

        for(int ff = 0; ff < Fr.size(); ff++)
                FrV[ff] = Fr[ff];

//	*min_freq = fv[get_min_freq(fv)];
	*min_freq = Fr.size();

	Fr.clear();
	fv.clear();
	
	return FrV;
}

/*double* map_to_array(map<int,double> Fr, double *min_freq, int *size)
{
        vector<double> fv;
        double *FrV;

        for(int fs = 0; fs < Fr.size(); fs++)
                fv.push_back(Fr[fs]);

        FrV = new double[fv.size()];

        MAX_EX = fv.size();

        for(int ff = 0; ff < Fr.size(); ff++)
                FrV[ff] = Fr[ff];

        *min_freq = fv[get_min_freq(fv)];

	*size = fv.size();

        Fr.clear();
        fv.clear();

        return FrV;
}*/


int get_min_freq(vector<double> fr)
{
	int pos;
	double min = INF;	
	
	for(int i = 0; i < fr.size(); i++)
	{
		if(fr[i] < min && fr[i] > 0)
		{
			min = fr[i];
			pos = i;
		}
	}

	return pos;
} 
	
int rex(int ex)
{
	return ex+1-S_SKIP;
}

bool qmode(string s)
{
	if(s.size() == 0)
		return true;

	istringstream str1(s);

	while(str1)
	{
		int n;
		str1 >> n;

		if(n <= 0)
			return true;

		if(n > 2)
			return false;
	}

	return true;
}

int qfind(vector<gene> genes, string s)
{
	for(int i = 0; i < genes.size(); i++)
	{
		if(genes[i].name == s)
			return i;
	}

	return -1;
}

m_point v_sum(vector<m_point> V, int n)
{
	m_point sum = 0;

	for(int i = 0; i < n; i++)
		sum += V[i];

	return sum;
}

double dev_std(double Avg, int n, vector<m_point> score)
{
	double Var = 0;

	for(int i = 0; i < score.size(); i++)
		Var += pow((double)score[i]-Avg,2);	

	Var *= (double)1/(double)n;

	return sqrt(Var);
}	

float percento(short int fat, short int div)
{
	return (float)(((float)fat/(float)div)*100);
}

const string tof(bool v)
{
	if(v)
		return "true";
	else
		return "false";
}

void tab_init()
{
	tabfile = queryfile;
	tabfile += ".tab";

	if(tabapp)
		return;

	ofstream out(tabfile.c_str());

	out << "Query = " << queryfile << "\nExons file = " << genesfile << "\nFrequency file = " << freqfile << endl
            << "-mc = " << MAX_CYCLE << "\tlocal = " << tof(local) << endl
	    << "-ss = " << S_SKIP << "\t-es = " << E_SKIP << endl 
	    << "-B = " << BEST << "\t -b = " << V_BEST
	    << "Opening gap = " << OGAP << "\tContinuing gap = " << CGAP << "\nSimilarity = " << SIM << endl << endl;

	out << "QID\tQEXN\tGID\tGEXN\tEX_M\tEX_M%\tQEX_M\tQEX_M%\tO_M\tO_M%\tTOT_M\tGAP\tSCORE\tEVALUE\n";

	out.close();
	return;
}

int v_pos_max(vector<m_point> v)
{
	int pos = 0;
	for(int i = 0; i < v.size(); i++)
		if(v[i] > v[pos])
			pos = i;
	return pos;
}

bool merge_check(vector<int> v, int m)
{
	for(int i = 0; i < v.size(); i++)
		if(v[i] == m)
			return true;

	return false;
}

/*long double get_pvalue(int filenum, m_point score, double df)
{

	df = 1; //comment this for gap correction; 
	long double p = exp(-(score)), pv;

	pv = get_coeff(1,filenum) * p;

//	cout << "\nscore/df = " << score/df << "\nprepv = " << pv;

//	pv *= pow(((double)1 - p),(TOT_EXONS[filenum]*df)-1);
//	pv = p;

//	pv *= df;

//	cout << "\npv = " << pv << "\np = " << p << "\ndf = " << df << "\nTOT_EXONS[filenum] = " << TOT_EXONS[filenum] << "\nScore = " << score << endl;
//	cout << "pow =" << pow(((double)1 - p),(TOT_EXONS[filenum]*df)-1) << endl;

	return pv;
}*/

/*m_point get_coeff(int df,int filenum)
{
	return (m_point)TOT_EXONS[filenum]*df;
}*/

#ifdef DEBUG10
void gene_view(gene g)
{
	cerr << "\nGene = " << g.name << endl;

	for(int i = 0; i < g.exL.size(); i++)
		cerr << g.exL[i] << "-";

	cerr << endl;

	return;
}
#endif

bool gene_unique(gene g, vector<gene> v)
{
	for(int i = 0; i < v.size(); i++)
		if(g.name == v[i].name && g.filenum == v[i].filenum && g.cdsS == v[i].cdsS && g.cdsE == v[i].cdsE && v[i].merged_with == g.merged_with)
			return false;

	return true;
}

vector<gene> v_gene_shrink(vector<gene> v)		
{
	vector<gene> s;

	for(int i = 0; i < v.size(); i++)
		if(gene_unique(v[i],s))
			s.push_back(v[i]);

	return s;
}


double factorial(double n)
{
	double f = 1;

	if(n == 0)
		return 1;

	while(n > 0)
	{
		f *= n;
		n--;
	}
	
	return f;
} 

bool iso_check(string s)
{
	istringstream str(s);
	string buf;
	bool flag = false;

	while(str)
	{
		str >> buf;

		if(buf != "S" && buf != "|" && buf != "G")
			return false;
		if(buf == "G")
			flag = true;
	}

	return flag;
}

/*int get_gaps(string s, m_point *score)
{
	istringstream str(s);
	string buf;
	bool flag = false;
	int count = 0;
	m_point lscore = *score;

	while(str)
	{
		str >> buf;

		if(!str)
			break;

		if(buf == "G" && flag)
                {
                        lscore -= CGAP;
                }

		else if(buf == "G" && !flag)
		{
			flag = true;
			lscore -= OGAP;	
			count++;
		}
		else
			flag = false;
	}

	*score = lscore;

	return count;
}*/

/*long double get_new_pvalue(int ex, int gap, m_point score, int tot_ex)
{
	long double p = exp(-score);

	double num = factorial((double)ex);	

	double den = factorial((double)gap);

	double tmp = ex - gap;

	if(tmp < 0)
		tmp = 0;

	tmp = factorial(tmp); 

	den *= tmp;

	num /= den;

	p *= num;

	p *= tot_ex;

	return p;
}*/

bool merge_gap_check(vector<string> vs)
{
	istringstream str0(vs[0]), str1(vs[1]), str2(vs[2]);
	
	while(str0)
	{
		string buf0,buf1,buf2;

		str0 >> buf0;
		str1 >> buf1;
		str2 >> buf2;

		if((buf0.find("M") != string::npos || buf2.find("M") != string::npos) && buf1 == "G")	
			return false;
	}

	return true;
}		

vector<string> string_to_vector(string s)
{
	istringstream is(s);
	vector<string> vs;

	for(int i = 0; i < S_SKIP; i++)
		vs.push_back(".");

	while(is)
	{
		string buf;

		is >> buf;

		if(buf.size())
			vs.push_back(buf);
	}

	for(int i = 0; i < E_SKIP; i++)
		vs.push_back(".");
	
	return vs;
}
	
bool fam_check(vector<string> v)
{
	if(v.size() <= 1)
		return true;

	for(int i = 0; i < v.size() - 1; i++)		
	{
		if(v[i] == v[v.size() - 1])
		{
//			cerr << endl << v[i] << '\t' << v[v.size() -1];
			return false;
		}
	}

	return true;
}
				
vector<int> dup_wash(vector<int> bbpos, vector<gene> bestg, vector<int> mmpos)
{
	vector<string> fams;
	vector<int> washed;

	for(int i = 0; i < bbpos.size(); i++)
	{
		fams.push_back(bestg[mmpos[bbpos[i]]].family);

		if(fam_check(fams))						
			washed.push_back(bbpos[i]);
	}

/*	for(int i = 0; i < washed.size(); i++)
	{
		cout << endl << "washed = " << washed[i] << '\t' << " gene = " << bestg.at(washed[i]).name;
	}*/

	return washed;
}

string launch_blastp(string qp, string qfileorg, string gp, string gfileorg, bool flag)
{
	ostringstream command, command2, tabpfile;

	string seq1, seq2;

	if(actual_query_p != qp)
		seq1 = seq_find(qp, qfileorg, false);
	else
	{
		ostringstream outfilename;
		outfilename << "ex_work_dir/." << QID << "." << qp << ".fasta";	
		seq1 = outfilename.str();
		actual_query_p = qp;
	}

	seq2 = seq_find(gp, gfileorg, false);

	if(flag)
	{
		command << "bl2seq -i " << seq1 << " -j " << seq2 << " -o " << HTML_DIR << "ex.out." << QID << "_" << RES_COUNT 
			<< ".pblast" << " -e 0.001 -p blastp" << " -W " << PWS;

		system(command.str().c_str());

		command2 << "bl2seq -i " << seq1 << " -j " << seq2 << " -o " << "ex_work_dir/." << "ex.out." << QID << "_" << RES_COUNT 
                         << ".pblast" << " -e 0.001 -p blastp -D 1" << " -W " << PWS;

		tabpfile << "ex_work_dir/." << "ex.out." << QID << "_" << RES_COUNT << ".pblast";
	}

	else
	{
		command2 << "bl2seq -i " << seq1 << " -j " << seq2 << " -o " << "ex_work_dir/." << "ex.out." << QID << "_" << "tmp" 
                         << ".pblast" << " -e 0.001 -p blastp -D 1" << " -W " << PWS;

		tabpfile << "ex_work_dir/." << "ex.out." << QID << "_" << "tmp" << ".pblast";
	}

	system(command2.str().c_str());

	
	return tabpfile.str();
}		

string launch_blastn(string qn, string qfileorg, string gn, string gfileorg, bool flag)
{
        ostringstream command, command2, tabnfile;

	string seq1, seq2;

	if(actual_query_n != qn)
        	seq1 = seq_find(qn, qfileorg, true);
	else
        {
                ostringstream outfilename;
                outfilename << "ex_work_dir/." << QID << "." << qn << ".fasta";
                seq1 = outfilename.str();
                actual_query_n = qn;
        }

        seq2 = seq_find(gn, gfileorg, true);


	if(flag)
	{
		command << "bl2seq -i " << seq1 << " -j " << seq2 << " -o " << HTML_DIR << "ex.out." << QID << "_" << RES_COUNT 
                        << ".nblast" << " -e 0.001 -p blastn" << " -W " << NWS;

        	system(command.str().c_str());

		command2 << "bl2seq -i " << seq1 << " -j " << seq2 << " -o " << "ex_work_dir/." << "ex.out." << QID << "_" << RES_COUNT 
                         << ".nblast" << " -e 0.001 -p blastn -D 1" << " -W " << NWS;

		tabnfile << "ex_work_dir/." << "ex.out." << QID << "_" << RES_COUNT << ".nblast";
	}

	else
	{
		command2 << "bl2seq -i " << seq1 << " -j " << seq2 << " -o " << "ex_work_dir/." << "ex.out." << QID << "_" << "tmp"
                         << ".nblast" << " -e 0.001 -p blastn -D 1" << " -W " << NWS;
		
		tabnfile << "ex_work_dir/." << "ex.out." << QID << "_" << "tmp" << ".nblast";

	}

	system(command2.str().c_str());

        return tabnfile.str();
}
	
string seq_find(string seqname, string seqfile,bool nucleotide)
{
	string line;

	if(seqname.find("CUSTOM") != string::npos)
		seqfile = "ex_work_dir/zzz_custom";

	if(nucleotide)
		seqfile += ".nuc";
	else
		seqfile += ".pro";

	ifstream in(seqfile.c_str());

	ostringstream outfilename;

	outfilename << "ex_work_dir/." << QID << "." << seqname << ".fasta";	

	ofstream out(outfilename.str().c_str());
	
        while(getline(in,line))
        {
                if(line[0] != '>')
                        continue;

                else
                    //    if(line.find(seqname) != string::npos)
			if(line.substr(1,line.find(".")-1) == seqname || (line.find(seqname) != string::npos && line.find(".") == string::npos))
                        {
                                out << line << endl;

                                while(getline(in,line))
                                {
                                        if(line[0] != '>')
                                                out << line;
                                        else
                                        {
                                                in.close();
						out.close();
                                                return outfilename.str();
                                        }
                                }

				in.close();
                                out.close();
                                return outfilename.str();
                        }
                        else
                                continue;
        }

        in.close();

        return "BADNESS";
}
	
int frame_min_freq(m_point *framefreq)
{
	int k = 0; 
	m_point min = framefreq[0];

	for(int i = 0; i < 3; i++)
	{
		if(framefreq[i] <= min)
		{
			min = framefreq[i];
			k = i;	
		}
	}

	return k;
}

unsigned short int blast_hit_counter(string filename)
{
        string line;
        ifstream in(filename.c_str());
        unsigned short int count = 0;

        while(getline(in,line))
        {
                if(line[0] != '#')
                        count++;
        }

        in.close();

        return count;
}

vector<string> post_process_merger(vector<string> al, gene *q, gene *g)
{
	int stq = S_SKIP + 1, stg = S_SKIP + 1;

//	for(int i = 0; i < al.size(); i++)
//		cout << al[i] << endl;

	if(local)
	{
		al[0] = brack_clean(al[0],&stq);
		al[2] = brack_clean(al[2],&stg);
	}

//	cout << endl;

	vector<vector<string> > val;

	for(int i = 0; i < 3; i++)
	{
		istringstream vss(al[i]);
		vector<string> vtmp;

		while(vss)
		{
			string tmp;
			vss >> tmp;

			if(!tmp.empty())
				vtmp.push_back(tmp);	
		}

		val.push_back(vtmp);
	}	

//	cerr << endl << "valsize = " << val.size() << endl;

	vector<string>::iterator si0 = val[0].begin(), si1 = val[1].begin(), si2 = val[2].begin();

	for(int i = 0; i < (int)val[1].size() - 3; i++)
	{
		if(val[1][i] == "|" && val[1][i+3] == "|" && val[0][i+1].find("M") == string::npos && val[0][i+2].find("M") == string::npos
							  && val[2][i+1].find("M") == string::npos && val[2][i+2].find("M") == string::npos)
		{
			if((val[1][i+1] == "G" 
                          && 
			  (val[1][i+2].find(":") != string::npos || val[1][i+2].find("*") != string::npos))
			|| 
			  (val[1][i+2] == "G" 
                          && 
                          (val[1][i+1].find(":") != string::npos || val[1][i+1].find("*") != string::npos)))	
			{
				int diff, sum,gap = 0;
				ostringstream str1;

				if(val[0][i+1] == "-")
				{
					sum = atoi(val[2][i+1].c_str())+atoi(val[2][i+2].c_str());
					diff = abs(atoi(val[0][i+2].c_str()) - sum);
					
					if(diff <= MAX_DIFF && diff%3 == 0)
					{
						val[0].erase(si0+1);
						*(si1+2) = ":.";
						val[1].erase(si1+1);

						for(int r = 0; r < i; r++)
							if(val[2][r] == "-")
								gap++;
						
						str1 << ' ' << (i+1-gap+S_SKIP+(stg-S_SKIP-1)) << '(' << g->exL[i+1-gap+S_SKIP+(stg-S_SKIP-1)] << ')'
                                     		     << '+' << (i+2-gap+S_SKIP+(stg-S_SKIP-1)) << '(' << g->exL[i+2-gap+S_SKIP+(stg-S_SKIP-1)] << ')';

						name_change(g,str1.str(),i+1);

						ostringstream os;
						os << sum << 'M';

						*(si2 + 2) = os.str();

						val[2].erase(si2+1);

						g->merged = true;
						g-> exN--;

						vector<int>::iterator vi = g->exL.begin();
						vector<short int>::iterator svi = g->exF.begin();

						for(int a = 0; a < i + 1 + S_SKIP; a++)
						{
							if(val[2][a] != "-")
							{
								vi++;
								svi++;
							}
						}

						*vi = sum;

						g->exL.erase(vi+1);			
						g->exF.erase(svi+1);
						g->cdsD[2]--;

					}
				}

				else if(val[0][i+2] == "-")
                                {
                                        sum = atoi(val[2][i+1].c_str())+atoi(val[2][i+2].c_str());
                                        diff = abs(atoi(val[0][i+1].c_str()) - sum);

                                        if(diff <= MAX_DIFF && diff%3 == 0)
                                        {
                                                val[0].erase(si0+2);
						*(si1+1) = ":.";
                                                val[1].erase(si1+2);

						for(int r = 0; r < i; r++)
                                                        if(val[2][r] == "-")
                                                                gap++;

						str1 << ' ' << (i+1-gap+S_SKIP+(stg-S_SKIP-1)) << '(' << g->exL[i+1-gap+S_SKIP+(stg-S_SKIP-1)] << ')'
                                                     << '+' << (i+2-gap+S_SKIP+(stg-S_SKIP-1)) << '(' << g->exL[i+2-gap+S_SKIP+(stg-S_SKIP-1)] << ')';

                                                name_change(g,str1.str(),i+1);

                                                ostringstream os;
                                                os << sum << 'M';

                                                *(si2 + 1) = os.str();

                                                val[2].erase(si2+2);

						g->merged = true;
                                                g-> exN--;

                                                vector<int>::iterator vi = g->exL.begin();
						vector<short int>::iterator svi = g->exF.begin();

                                                for(int a = 0; a < i + 1 + S_SKIP; a++)
						{
							if(val[2][a] != "-")
							{
                                                        	vi++;
								svi++;
							}
						}

                                                *vi = sum;

                                                g->exL.erase(vi+1);
						g->exF.erase(svi+1);
						g->cdsD[2]--;
                                        }
                                }

				else if(val[2][i+1] == "-")
                                {
                                        sum = atoi(val[0][i+1].c_str())+atoi(val[0][i+2].c_str());
                                        diff = abs(atoi(val[2][i+2].c_str()) - sum);

                                        if(diff <= MAX_DIFF && diff%3 == 0)
                                        {
                                                val[2].erase(si2+1);
						*(si1+2) = ":.";
                                                val[1].erase(si1+1);

						for(int r = 0; r < i; r++)
                                                        if(val[0][r] == "-")
                                                                gap ++;
	
						str1 << ' ' << (i+1-gap+S_SKIP+(stq-S_SKIP-1)) << '(' << q->exL[i+1-gap+S_SKIP+(stq-S_SKIP-1)] << ')'
                                                     << '+' << (i+2-gap+S_SKIP+(stq-S_SKIP-1)) << '(' << q->exL[i+2-gap+S_SKIP+(stq-S_SKIP-1)] << ')';

                                                name_change(q,str1.str(),i+1);

                                                ostringstream os;
                                                os << sum << 'M';

                                                *(si0 + 2) = os.str();

                                                val[0].erase(si0+1);

						q->merged = true;
                                                q-> exN--;

                                                vector<int>::iterator vi = q->exL.begin();
						vector<short int>::iterator svi = q->exF.begin();

                                                for(int a = 0; a < i + 1 + S_SKIP; a++)
						{
							if(val[0][a] != "-")
							{
                                                        	vi++;
								svi++;
							}
						}

                                                *vi = sum;

                                                q->exL.erase(vi+1);
						q->exF.erase(svi+1);
						q->cdsD[2]--;
                                        }
                                }	

				else if(val[2][i+2] == "-")
                                {
                                        sum = atoi(val[0][i+1].c_str())+atoi(val[0][i+2].c_str());
                                        diff = abs(atoi(val[2][i+1].c_str()) - sum);

                                        if(diff <= MAX_DIFF && diff%3 == 0)
                                        {
                                                val[2].erase(si2+2);
						*(si1+1) = ":.";
                                                val[1].erase(si1+2);

						for(int r = 0; r < i; r++)
                                                        if(val[0][r] == "-")
                                                                gap ++;

						str1 << ' ' << (i+1-gap+S_SKIP+(stq-S_SKIP-1)) << '(' << q->exL[i+1-gap+S_SKIP+(stq-S_SKIP-1)] << ')'
                                                     << '+' << (i+2-gap+S_SKIP+(stq-S_SKIP-1)) << '(' << q->exL[i+2-gap+S_SKIP+(stq-S_SKIP-1)] << ')';

						name_change(q,str1.str(),i+1);

                                                ostringstream os;
                                                os << sum << 'M';

                                                *(si0 + 1) = os.str();

                                                val[0].erase(si0+2);

						q->merged = true;
                                                q-> exN--;

                                                vector<int>::iterator vi = q->exL.begin();
						vector<short int>::iterator svi = q->exF.begin();

                                                for(int a = 0; a < i + 1 + S_SKIP; a++)
						{
							if(val[0][a] != "-")
							{
                                                        	vi++;
								svi++;
							}
						}

                                                *vi = sum;

                                                q->exL.erase(vi+1);
						q->exF.erase(svi+1);
						q->cdsD[2]--;
                                        }
                                }
			}	
		}

		si0++;
		si1++;
		si2++;
	}	
	
	si0 = val[0].begin(); 
	si1 = val[1].begin(); 
	si2 = val[2].begin();

/*	for(int i = 0; i < (int)val[1].size() - 3; i++)
	{
		if(val[1][i] == "|" && val[1][i+3] == "|" && val[0][i+1].find("M") == string::npos && val[0][i+2].find("M") == string::npos
							  && val[2][i+1].find("M") == string::npos && val[2][i+2].find("M") == string::npos)
		{
			if(((val[1][i+1].find(":") != string::npos || val[1][i+1].find("*") != string::npos) 
                          && 
			  (val[1][i+2].find(":") != string::npos || val[1][i+2].find("*") != string::npos)))
			{
				int  sum,gapg = 0,gapq = 0;
				ostringstream str1,str2;

				sum = atoi(val[2][i+1].c_str())+atoi(val[2][i+2].c_str());
					
				if(sum == atoi(val[0][i+1].c_str())+atoi(val[0][i+2].c_str()))
				{
					*(si1+2) = "|";
					val[1].erase(si1+1);

					for(int r = 0; r < i; r++)
                                                if(val[0][r] == "-")
                                                        gapq++;

					for(int r = 0; r < i; r++)
						if(val[2][r] == "-")
							gapg++;
					
					str1 << ' ' << (i+1-gapq+S_SKIP+(stq-S_SKIP-1)) << '(' << q->exL[i+1-gapq+S_SKIP+(stq-S_SKIP-1)] << ')'
                                                     << '+' << (i+2-gapq+S_SKIP+(stq-S_SKIP-1)) << '(' << q->exL[i+2-gapq+S_SKIP+(stq-S_SKIP-1)] << ')';

					name_change(q,str1.str(),i+1);

					str2 << ' ' << (i+1-gapg+S_SKIP+(stg-S_SKIP-1)) << '(' << g->exL[i+1-gapg+S_SKIP+(stg-S_SKIP-1)] << ')'
                                    		     << '+' << (i+2-gapg+S_SKIP+(stg-S_SKIP-1)) << '(' << g->exL[i+2-gapg+S_SKIP+(stg-S_SKIP-1)] << ')';

					name_change(g,str2.str(),i+1);

					ostringstream os;
					os << sum << 'M';

					*(si0 + 2) = os.str();
					*(si2 + 2) = os.str();

					val[0].erase(si0+1);
					val[2].erase(si2+1);

					q->merged = true;
					q-> exN--;

					g->merged = true;
					g-> exN--;

					vector<int>::iterator vi = g->exL.begin();
					vector<short int>::iterator svi = g->exF.begin();

					for(int a = 0; a < i + 1 + S_SKIP; a++)
					{
						if(val[2][a] != "-")
						{
							vi++;
							svi++;
						}
					}

					*vi = sum;

					g->exL.erase(vi+1);			
					g->exF.erase(svi+1);
					g->cdsD[2]--;

					vi = q->exL.begin();
                                        svi = q->exF.begin();

                                        for(int a = 0; a < i + 1 + S_SKIP; a++)
                                        {
                                                if(val[0][a] != "-")
                                                {
                                                        vi++;
                                                        svi++;
                                                }
                                        }

                                        *vi = sum;

                                        q->exL.erase(vi+1);
                                        q->exF.erase(svi+1);
                                        q->cdsD[2]--;
				}
			}
		}

		si0++;
		si1++;
		si2++;
	}	*/

	if(local)
	{
		for(int i = 0; i < 3; i++)
        	{
                	ostringstream voss;

			if(i == 0)
			{
				int gcount = 0;
				for(int j = 0; j < val[i].size(); j++)
				{
					if(val[i][j] != "-")
                                		voss << val[i][j] << '[' << stq + j - gcount << ']' << ' ';
					else
					{
						voss << val[i][j] << ' ';
						gcount++;
					}
				}
			}
			
			else if(i == 1)
			{
                		for(int j = 0; j < val[i].size(); j++)
                        		voss << val[i][j] << ' ';
			}

			else if(i == 2)
                        {
				int gcount = 0;
                                for(int j = 0; j < val[i].size(); j++)
				{
					if(val[i][j] != "-")
                                        	voss << val[i][j] << '[' << stg + j - gcount << ']' << ' ';
					else
					{
						voss << val[i][j] << ' ';
						gcount++;
					}
				}
                        }

                	al[i] = voss.str();
        	}
	}
	else
	{
		for(int i = 0; i < 3; i++)
        	{
			ostringstream voss;

			for(int j = 0; j < val[i].size(); j++)
				voss << val[i][j] << ' ';

			al[i] = voss.str();
        	}
	}

//	for(int i = 0; i < al.size(); i++)
  //              cout << al[i] << endl;

//	cout << endl;

	return al;
}

string brack_clean(string s, int *st)
{
	*st = -1;
	string n;

	for(int i = 0; i < s.size(); i++)
	{
		if(s[i] == '[')
		{
			string sts;
			i++;

			while(s[i] != ']')
			{
				sts.push_back(s[i]);
				i++;
			}

			if(*st == -1)
			{
				istringstream ss(sts);
			
				ss >> *st;
			}
			else
				continue;
		}

		else
			n.push_back(s[i]);
	}

	return n;
}

void name_change(gene *g, string str1, int pos)
{
	if(g->merg_pos.size() == 0)
        {
     		g->merg_pos.push_back(g->exL.size() - 1 - S_SKIP);
                g->name += str1;
        }

        else
        {
       	        for(int bb = 0; bb < g->merg_pos.size(); bb++)
        		if(g->merg_pos[bb] >= pos - 1 - S_SKIP)
                        {
                        	istringstream ibuf(g->name);
                                ostringstream buf1;
                                string sbuf;

                                ibuf >> sbuf;

                                buf1 << sbuf << ' ' << str1 << '-';

                                ibuf >> sbuf;

				buf1 << sbuf;

				g->name = buf1.str();

				for(int gg = 0; gg < g->merg_pos.size(); gg++)
				{
					if((g->merg_pos[gg] >= g->exL.size() - 1 - S_SKIP))
						g->merg_pos[gg]--;
				}

				break;
			}
			else         //new
			{
				g->name += "-";
				g->name += str1.substr(1);
			}

		g->merg_pos.push_back(g->exL.size() - 1 - S_SKIP);
	}
	
	return;
}

bool merged_with_check(string s,vector<string> *v)
{
	for(int i = 0; i < (int)v->size(); i++)
	{
		if(v->at(i) == s)
			return true;
	}

	return false;
}
