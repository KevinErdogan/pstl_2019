#include <stdio.h>
#include <vector>
#include <iostream>
#include <utility>
#include <sstream>
#include <gmpxx.h>
#include <sys/param.h>
#include <sys/times.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <map>
#include <time.h>
#include <sys/time.h>
// compiler avec g++ unrank_bin_tree_non_plane.cpp -lgmpxx -lgmp -o btree

using namespace std;


typedef vector <mpz_class> enum_vect;
struct node
{
	node * fg;
	node * fd;
};

typedef struct tv{
	mpz_class taille;
	int nbNS;
	mpz_class PTST;
	tv *fg;
	tv *fd;
} treeVals;

typedef struct tvB{
	vector<bool> bs;
	mpz_class taille;
	int nbNS;
	mpz_class PTST; 
	tvB *fg;
	tvB *fd;
} treeValsBits;

typedef struct iso{
	int label;
	vector<int> vec;
	int nodeIsoValue;
	mpz_class taille;
	int nbNS;
	mpz_class PTST;
	int level;
	iso *pere;
	iso *fg;
	iso *fd;
} treeIso;



typedef struct{
	int maxLevel;
	vector<map<vector<int>, int>> maps;
}isoMaps;

typedef struct tree{
 tree* fg;
 tree* fd;
 int weight;
 int label;
} tree;

// make the list of numbers of trees from 1 up to n
enum_vect count(int n);
// pretty print the vector
string vect_str(enum_vect G);
// generate the tree number r with size n
node * tree_gen(int n, enum_vect G, mpz_class r);

node * leaf_init();
node * binaire(node * T, node * U);
mpz_class rand(mpz_class min, mpz_class max);
string dot_from_tree(clock_t a, node * A);
void store(int n, node * A, mpz_class r);

////////////:

//HNonPlan sans AHU
mpz_class HNonPlan(node *tree);
treeVals * tv_init(mpz_class taille, int nbNS, mpz_class PTST, treeVals *fg, treeVals *fd);
void freeTreeVals (treeVals *T);
bool isLeaf(node *n);
treeVals * HNonPlanAux(node *tree);
bool equalsTreeVals(treeVals *fg, treeVals *fd);
mpz_class factorial(mpz_class n);

//HNonPlan avec AHU
mpz_class HNonPlanIso(node *tree);
treeIso *ti_init(mpz_class taille, int nbNS, mpz_class PTST, int level, treeIso *fg, treeIso *fd);
void freeTreeIso(treeIso *ti);
treeIso *HNonPlanIsoAux(node *tree, isoMaps **m, int level);
bool equalsTreeIso(treeIso *fg, treeIso *fd);
isoMaps *init_isoMaps();
void resizeisoMaps(isoMaps **m , int levelmax);
int getNodeValueFromMaps(isoMaps **m, int level, vector<int> v);
void free_isoMaps(isoMaps *m);

//generation aleatoire uniforme
treeIso* genAleaAvecNS(node *tree);
void sousCalculAvecNS(treeIso** t, treeIso*** M, int taille, int* label);
string dot_from_treeIso(clock_t a, treeIso * A);
void storeLabeledIso(int n, treeIso * A, mpz_class m);//A doit etre etiquete
void mkCanonique(treeIso ** t);
bool sameCanonique(treeIso *t1, treeIso *t2);

int main_a()
{
	int n;
	for(n=0; n<8; n++)
	{
		enum_vect B = count(n);
		//dans B[n] : nb d'arbres à 2n+1 noeuds
	
		//cout << vect_str(B) << endl;
		mpz_class r;
		//r = rand(0, B[n]-1);
		for(int r=0; r<B[n]; r++)
		{
			node * A = tree_gen(2*n+1, B, r);
			cout << "l'arbre de taille " << 2*n+1 << " a " << HNonPlan(A) << " etiquetages possibles (sans AHU)" << endl;
			cout << "l'arbre de taille " << 2*n+1 << " a " << HNonPlanIso(A) << " etiquetages possibles (avec AHU)" << endl;
			//store(n, A, r);
		}
		cout << endl;
	}
	return 0; 
}

int main_b()//experimentation HNonPlan
{
	
	FILE *f;
	FILE *fiso;
	
	clock_t initial_time;
	clock_t final_time;
	double cpu_time=0.0;
	double cpu_time_Iso=0.0;
	int n, round;
	for(n = 1000; n <= 1010; n++){
		cpu_time = 0.0;
		cpu_time_Iso=0.0;
		initial_time = clock();
		enum_vect B = count(n);
		final_time = clock();
		cpu_time += (((double)(final_time - initial_time))/CLOCKS_PER_SEC);
		for(round = 0; round < 10; round++){
			
			mpz_class r = rand(0, B[n]-1);
			initial_time = clock();
			node * A = tree_gen(2*n+1, B, r);
			final_time = clock();
			cpu_time += (((double)(final_time - initial_time))/CLOCKS_PER_SEC);
		
			initial_time = clock();
			HNonPlan(A);
			final_time = clock();
			cpu_time += (((double)(final_time - initial_time))/CLOCKS_PER_SEC);
			
			initial_time = clock();
			HNonPlanIso(A);
			final_time = clock();
			cpu_time_Iso += (((double)(final_time - initial_time))/CLOCKS_PER_SEC);
		}
		cpu_time /= 10;
		cpu_time_Iso /= 10;
		cout << "temps moyen pour n = " << 2*n+1 << " : " << cpu_time << endl;
		cout << "temps moyen pour AHU n = " << 2*n+1 << " : " << cpu_time_Iso << endl;	
		cout << endl;	
	}
	return 0;
}

int main(){//main generation aleatoire uniforme (test avec l'arbre de taille 5)

	int n = 2;
	enum_vect B = count(n);
	mpz_class r;
	r = 0;
	node * A = tree_gen(2*n+1, B, r);
	mpz_class hook = HNonPlanIso(A);
	cout << "l'arbre de taille " << 2*n+1 << " a " << hook << " etiquetages possibles" << endl;

	int nb_distinct_gen = 0;

	cout << "Generation de 40000 arbres de taille 5" << endl;

	std::srand(time(NULL));
	treeIso *tc0 = genAleaAvecNS(A);
	tc0->fg->label = 2;
	tc0->fd->label = 3;
	tc0->fd->fg->label = 4;
	tc0->fd->fd->label = 5;
	treeIso *tc1 = genAleaAvecNS(A);
	tc1->fg->label = 3;
	tc1->fd->label = 2;
	tc1->fd->fg->label = 4;
	tc1->fd->fd->label = 5;
	treeIso *tc2 = genAleaAvecNS(A);
	tc2->fg->label = 4;
	tc2->fd->label = 2;
	tc2->fd->fg->label = 3;
	tc2->fd->fd->label = 5;
	treeIso *tc3 = genAleaAvecNS(A);
	tc3->fg->label = 5;
	tc3->fd->label = 2;
	tc3->fd->fg->label = 3;
	tc3->fd->fd->label = 4;

	int nbT0 = 0;
	int nbT1 = 0;
	int nbT2 = 0;
	int nbT3 = 0;
	for(int round = 0; round < 40000; round++){
		treeIso *T = genAleaAvecNS(A);

		if(sameCanonique(T, tc0))
			nbT0++;
		if(sameCanonique(T, tc1))
			nbT1++;
		if(sameCanonique(T, tc2))
			nbT2++;
		if(sameCanonique(T, tc3))		
			nbT3++;
	}		
	cout << "Canonique 1 : " << nbT0 << "\nCanonique 2 : " << nbT1 << "\nCanonique 3 : " << nbT2 << "\nCanonique 4 : " << nbT3 << endl;
		
	return 0;
}

enum_vect count(int n)
{
	enum_vect B;
	B.push_back(1);
	
	for(int i=1; i<n+1; i++)
	{
		// on compte les arbres de taille 2i+1
		// (sum_k=0..i-1 T_2k+1 * T_2i-2k-1)/2 + T_i/2
		// 
		mpz_class t = 0;
		int mid = floor((i-1) / 2);
		for(int k=0; k<=mid; k++)
		{
			t = t + B[k]*B[i-k-1];
		}
		if(mid == i-1-mid)
		{
			mpz_class res=2;
			mpz_class p = B[mid]*B[mid];
			mpz_fdiv_q(res.get_mpz_t(), p.get_mpz_t(), res.get_mpz_t());
			t = t - res;
			res = 2;
			mpz_fdiv_q(res.get_mpz_t(), B[mid].get_mpz_t(), res.get_mpz_t());
			t = t + res;
		}
		
		B.push_back(t);
	}
	return B;
}


string vect_str(enum_vect G)
{
	stringstream ss;
	ss.str("");
	ss << "[";
	for(int i = 0; i < G.size(); i++)
	{
  		if(i != 0)
		    ss << ",";
		ss << G[i];
	}
	ss << "]";
	
	return ss.str();
}


node * leaf_init()
{
	node * A;
	A = new node;
	
	A->fg = NULL;
	A->fd = NULL;
	return A;
}


node * binaire(node * T, node * U)
{
	node * A;
	A = new node;
	
	A->fg = T;
	A->fd = U;
	return A;
}


mpz_class rand(mpz_class min, mpz_class max)
{
	gmp_randstate_t state;
	gmp_randinit_default(state);
	
	mpz_class M = time(NULL);
	gmp_randseed(state, M.get_mpz_t());
	
	mpz_class r;
		
	M = max-min+1;
	
	mpz_urandomm(r.get_mpz_t(), state , M.get_mpz_t());
	r = r + min;
	
	return r;
}


 node * tree_gen(int n, enum_vect B, mpz_class r)
 {
 	if(n==1)
 		return leaf_init();
 	
 	mpz_class res=0;
 	int k=1;
 	int i = floor((n-1)/2);
 	while(r>=0)
 	{
		r = r - B[floor((k-1)/2)]*B[floor((n-1-k)/2)];
		if(k==i)
		{
			res=2;
			mpz_fdiv_q(res.get_mpz_t(), B[floor((i-1)/2)].get_mpz_t(), res.get_mpz_t());
			r = r + res;
		}
		k = k+2;
	}
	k = k-2;
	r = r + B[floor((k-1)/2)]*B[floor((n-1-k)/2)];
	if(k==i)
	{
		r = r - res;
	}
	
	node * f1;
	node * f2;

	mpz_class quot,rest;
	if(k<i)
	{
		mpz_fdiv_qr(quot.get_mpz_t(), rest.get_mpz_t(), r.get_mpz_t(), B[floor((k-1)/2)].get_mpz_t());
		f1 = tree_gen(k, B, rest);
		f2 = tree_gen(n-1-k, B, quot);
	} 
	else
	{
		if(r<B[floor((i-1)/2)])
		{
			mpz_fdiv_qr(quot.get_mpz_t(), rest.get_mpz_t(), r.get_mpz_t(), B[floor((i-1)/2)].get_mpz_t());
			f1 = tree_gen(i, B, rest);
			f2 = tree_gen(i, B, rest);
		}
		else
		{
			r = r-B[floor((i-1)/2)];
			mpz_fdiv_qr(quot.get_mpz_t(), rest.get_mpz_t(), r.get_mpz_t(), B[floor((i-1)/2)].get_mpz_t());
			mpz_class m, M;
			if(rest<quot)
			{
				m = rest;
				M = quot;
			}
			else
			{
				m = quot;
				M = rest;
			}
			
			f1 = tree_gen(i, B, m);
			f2 = tree_gen(i, B, M+1);
		}
	}
	return binaire(f1, f2);	
 }


string dot_from_tree(clock_t a, node * A)
{
	stringstream ss;
	
	if(A!=NULL)
	{
		ss << a << "[label=\"\",shape=point]\n";
		
		clock_t g,d;
		g = clock();
		if(A->fg !=NULL)
		{
			ss << g << "[label=\"\",shape=point];\n";
			ss << "  " << a << " -- " << g << ";\n";
			ss << dot_from_tree(g, A->fg)<< "\n";
		}
		
		if(A->fd !=NULL)
		{
			d = clock();
			while(d==g)
			{
				d = clock();
			}
			ss << d << "[label=\"\",shape=point];\n";
			ss << "  " << a << " -- " << d << ";\n";
			ss << dot_from_tree(d, A->fd) << "\n";
		}
		return ss.str();
	}
	return "";
}



void store(int n, node * A, mpz_class r)
{
	ofstream fic;
	stringstream ss;
	ss << "test/Bin_" << n << "_" << r << ".dot";
	fic.open(ss.str(), ios::out);
	
	fic << "graph {\n";
	fic << "node [shape=plaintext]\n";
	clock_t c = clock();
	fic << dot_from_tree(c, A);
	fic << "}\n";
	
	fic.close();	
}

/////////////////////BINAIRE/////////////////////////////

////Version sans AHU////

treeVals * tv_init(mpz_class taille, int nbNS, mpz_class PTST, treeVals *fg, treeVals *fd)
{
	treeVals * tv;
	tv = new treeVals;
	tv->taille = taille;
	tv->nbNS = nbNS;
	tv->PTST = PTST;
	tv->fg = fg;
	tv->fd = fd;
	return tv;
}

void freeTreeVals (treeVals *T){
	if(T->fg !=NULL){
		freeTreeVals(T->fg);
		freeTreeVals(T->fd);
	}
	delete T;
}

bool isLeaf(node *n){
	return (n->fg == NULL && n->fd == NULL);
}

mpz_class HNonPlan(node *tree){
	treeVals *rootTV = HNonPlanAux(tree);
	mpz_class fact = factorial(rootTV->taille);
	int powNS = (int) pow(2, rootTV->nbNS);
	mpz_class powMPZ(powNS);
	mpz_class PTST = rootTV->PTST;
	freeTreeVals(rootTV);
	return fact / PTST / powMPZ;
}

treeVals *HNonPlanAux(node *tree){
	if(isLeaf(tree)){
		return tv_init(1, 0, 1, NULL, NULL);
	}else{
		treeVals *fg = HNonPlanAux(tree->fg);
		treeVals *fd = HNonPlanAux(tree->fd);
		int nbNS;
		if(equalsTreeVals(fg, fd)){
			nbNS = fg->nbNS + fd->nbNS + 1;
		}else{
			nbNS = fg->nbNS + fd->nbNS;
		}
		mpz_class taille = fg->taille + fd->taille + 1;
		return tv_init(taille, nbNS, fg->PTST * fd->PTST * taille, fg, fd);
	}
}

bool equalsTreeVals(treeVals *fg, treeVals *fd){
	if(fg->taille != fd->taille){
		return false;
	}else{
		if(fg->taille == 1){
			return true;
		}else{
			return ( (equalsTreeVals(fg->fg, fd->fg) && equalsTreeVals(fg->fd, fd->fd))
					|| (equalsTreeVals(fg->fg, fd->fd) && equalsTreeVals(fg->fd, fd->fg)) );
		}
	}
}
////////////////////////

mpz_class factorial(mpz_class n)
{
	if(n <= 0) return 1;
	mpz_class res(n);
    while(n-- > 1) res *= n;
    return res;
}

////Version avec AHU////

isoMaps *init_isoMaps(){
	isoMaps *m = new isoMaps;
	m-> maxLevel = 1;
	m->maps.reserve(1);
	return m;
}

void free_isoMaps(isoMaps *m){
	if(m != NULL){
		delete m;
	}
}

void resizeisoMaps(isoMaps **m , int levelmax){//appeler sur une feuille dans HNonPlanIsoAux
	if(levelmax > ((*m)->maxLevel)){
		(*m)->maps.resize(levelmax); 
		(*m)->maxLevel = levelmax;
	}
}

int getNodeValueFromMaps(isoMaps **m, int level, vector<int> v){
	auto it = (*m)->maps[level].find(v);
	if(it != (*m)->maps[level].end()){
		return it->second;
	}else{//si le vecteur n'est pas present alors on modifie m
		int nextValue = (*m)->maps[level].size() +1;
		(*m)->maps[level].insert(pair<vector<int>, int>(v, nextValue));
		return nextValue;
	}
}

treeIso *ti_init(mpz_class taille, int nbNS, mpz_class PTST, int level, treeIso *fg, treeIso *fd){//nodeValue et vector init dans HNonPlanIsoAux
	treeIso * t;
	t = new treeIso;
	t->taille = taille;
	t->nbNS = nbNS;
	t->PTST = PTST;
	t->pere = NULL;
	t->fg = fg;
	t->fd = fd;
	t->level = level;
	t->label = -1;
	return t;
}

void freeTreeIso(treeIso *T){
	if(T->fg !=NULL){
		freeTreeIso(T->fg);
		freeTreeIso(T->fd);
	}
	delete T;
}


mpz_class HNonPlanIso(node *tree){
	isoMaps *m = init_isoMaps();
	treeIso *rootTV = HNonPlanIsoAux(tree, &m, 1);
	mpz_class fact = factorial(rootTV->taille);
	int powNS = (int) pow(2, rootTV->nbNS);
	mpz_class powMPZ(powNS);
	mpz_class PTST = rootTV->PTST;
	freeTreeIso(rootTV);
	free_isoMaps(m);
	return fact / PTST / powMPZ;
}

treeIso *HNonPlanIsoAux(node *tree, isoMaps **m, int level){
	if(isLeaf(tree)){
		treeIso *leaf = ti_init(1, 0, 1, level, NULL, NULL);
		leaf->nodeIsoValue = 0;
		resizeisoMaps(m, level);
		return leaf;
	}else{
		treeIso *fg = HNonPlanIsoAux(tree->fg, m, level+1);
		treeIso *fd = HNonPlanIsoAux(tree->fd, m, level+1);
		mpz_class taille = fg->taille + fd->taille +1;
		mpz_class PTST = fg->PTST * fd->PTST * taille;
		vector<int> v;
		if(fg->nodeIsoValue > fd->nodeIsoValue){//v[0] <= v[1]
			v.emplace_back(fd->nodeIsoValue);
			v.emplace_back(fg->nodeIsoValue);
		}else{
			v.emplace_back(fg->nodeIsoValue);
			v.emplace_back(fd->nodeIsoValue);
		}
		int nodeValue = getNodeValueFromMaps(m, level, v);
		int nbNS;
		if(equalsTreeIso(fg, fd)){
			nbNS = fg->nbNS + fd->nbNS + 1;
		}else{
			nbNS = fg->nbNS + fd->nbNS;
		}
		treeIso *n = ti_init(taille, nbNS, PTST, level, fg, fd);
		fg->pere = n;
		fd->pere = n;
		n->nodeIsoValue = nodeValue; 
		return n;
	}
}

bool equalsTreeIso(treeIso *fg, treeIso *fd){
	return fg->nodeIsoValue == fd->nodeIsoValue;
}
///////////////////////////

// generation aleatoire uniforme d'un arbre canonique

void sousCalculAvecNS(treeIso** t, treeIso*** M, int taille, int* label){
   if((*t) == NULL) return;
   int count = 0;
   for(int i=0; i< taille; i++)
	if( M[i] != NULL)
	     count++;
   // random r
   int r;
   if(count > 0){
   	r = ((std::rand()) % count) ;
   }else{
	r = 0;
   }

   // found the r-th tree* and put on his label the value of *label
   int nb=0;
   treeIso** v;
   for(int i=0; i< taille; i++){
	if( M[i] != NULL){
	     if(nb==r){
	          v = M[i];
	          break;
	     }
	     else{
	          nb++;
	     }
        }
   }
   (*v)->label = *label;
   *label = (*label)+1;	

   // construct M and recursive call
    int firstPlaceOfV = -1;
    // suppress v
    for(int i=0; i < taille; i++){
         if(M[i]==v){
	      if(firstPlaceOfV == -1)
	           firstPlaceOfV = i;
	      M[i]=NULL;
	 }
    }

    if( (*v)-> fg!= NULL && (*v)-> fd!=NULL){//t is binary strict
	 int fg_pos_in_M = (*v)->fg->taille.get_si() + firstPlaceOfV;
         for(int i=firstPlaceOfV; i < fg_pos_in_M; i++){
	      M[i]=&((*v)->fg);
         }
	 int fd_pos_in_M = (*v)->fd->taille.get_si() + firstPlaceOfV + (*v)->fg->taille.get_si();
	 for(int i = firstPlaceOfV + (*v)->fg->taille.get_si(); i < fd_pos_in_M; i++){
	      M[i]=&((*v)->fd);
         }
	 // call
	 sousCalculAvecNS(&((*v)->fg), M, taille, label);
	 sousCalculAvecNS(&((*v)->fd), M, taille, label);
    }
}


//generation aleatoire uniforme
//retourne un arbre binaire avec etiquetage canonique
treeIso* genAleaAvecNS(node *tree){
	if(tree == NULL) return NULL;
	int *label = (int *)malloc(sizeof(int));
	*label = 1;
	isoMaps *m = init_isoMaps();
	treeIso *rootTV = HNonPlanIsoAux(tree, &m, 1);
	int taille = (rootTV->taille).get_si();
	treeIso ***M = (treeIso***)malloc(sizeof(treeIso**)*taille);//attention la taille de l'arbre ne doit pas être trop grande
	for(int i = 0; i < taille; i++)
		M[i] = &rootTV;

	sousCalculAvecNS(&rootTV, M, taille, label);

	free_isoMaps(m);
	free(M);
	free(label);

	mkCanonique(&rootTV);	
	return rootTV;
}

void mkCanonique(treeIso ** t){//init une seed avant d'appeler mkCanonique
	if((*t)->fg == NULL) return;
	//Si t est un NS alors on echange les places de fg et fd de façon à ce que l'arbre ait un etiquetage canonique
	if((*t)->nbNS != ((*t)->fg->nbNS + (*t)->fd->nbNS)){
	     if((*t)->fg->label > (*t)->fd->label){//si le label du fg est superieur au label du fd
		  treeIso *tmp = (*t)->fg;
		  (*t)->fg = (*t)->fd;
		  (*t)->fd = tmp;
	     }
	}
   	mkCanonique(&((*t)->fg));
	mkCanonique(&((*t)->fd));
}

bool sameCanonique(treeIso *t1, treeIso *t2){
	if(t1->taille != t2->taille) return false;
	if(t1->label != t2 -> label) return false;
	
	if(t1->taille == 1) return true;

	return sameCanonique(t1->fg, t2->fg) && sameCanonique(t1->fd, t2->fd);

}

string dot_from_treeIso(clock_t a, treeIso * A)
{
	stringstream ss;
	
	if(A!=NULL)
	{
		
		ss << a << "[label=" << A->label << ",shape=ellipse]\n";
		
		clock_t g,d;
		g = clock();
		if(A->fg !=NULL)
		{
			ss << g << "[label=\"\",shape=ellipse];\n";
			ss << "  " << a << " -- " << g << ";\n";
			ss << dot_from_treeIso(g, A->fg)<< "\n";
		}
		
		if(A->fd !=NULL)
		{
			d = clock();
			while(d==g)
			{
				d = clock();
			}
			ss << d << "[label=\"\",shape=ellipse];\n";
			ss << "  " << a << " -- " << d << ";\n";
			ss << dot_from_treeIso(d, A->fd) << "\n";
		}
		return ss.str();
	}
	return "";
}

void storeLabeledIso(int n, treeIso * A, mpz_class m)
{
	ofstream fic;
	stringstream ss;
	ss << "testLabeled/Bin_" << n << "_" << m << ".dot";
	fic.open(ss.str(), ios::out);
	
	fic << "graph {\n";
	fic << "node [shape=plaintext]\n";
	clock_t c = clock();
	fic << dot_from_treeIso(c, A);
	fic << "}\n";
	
	fic.close();	
}

