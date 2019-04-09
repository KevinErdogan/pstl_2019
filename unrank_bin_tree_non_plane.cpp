
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

// compiler avec g++ XXX.cpp -lgmpxx -lgmp -o XXX

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
	mpz_class PTST; // produit de la taille de tous les sous arbres
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
	int vec[2];
	int nodeIsoValue;
	mpz_class taille;
	int nbNS;
	mpz_class PTST;
	unsigned int level;
	iso *root;
	iso *fg;
	iso *fd;
} treeIso;

typedef struct{
	int maxLevel;
	int ***vecs;//vecs[level][iso.vec[0]][iso.vec[1]] = iso.nodeIsoValue
	int *vecsSize;
	int *maxNodeIsoValueByLevel;
}isoVecs;

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

treeVals * tv_init(mpz_class taille, int nbNS, mpz_class PTST, unsigned int level, treeIso *fg, treeIso *fd, treeIso *root);
void freeTreeVals (treeVals *T);
bool isLeaf(node *n);
mpz_class HNonPlan(node *tree);
treeVals * HNonPlanAux(node *tree);
bool equalsTreeVals(treeVals *fg, treeVals *fd);
mpz_class factorial(mpz_class n);

//
treeValsBits * tv_initBits(mpz_class taille, int nbNS, mpz_class PTST, treeValsBits *fg, treeValsBits *fd);
void freeTreeValsBits (treeValsBits *T);
mpz_class HNonPlanBits(node *tree);
treeValsBits * HNonPlanBitsAux(node *tree);
bool equalsTreeValsBits(treeValsBits *fg, treeValsBits *fd);

//iso
treeIso *ti_init(mpz_class taille, int nbNS, mpz_class PTST, int level, treeIso *fg, treeIso *fd, treeIso *root);
void freeTreeIso(treeIso *ti);
mpz_class HNonPlanIso(node *tree);
treeIso *HNonPlanIsoAux(node *tree, node *root);
bool equalsTreeIso(treeIso *fg, treeIso *fd);
isoVecs *init_isoVecs();
void resizeIsoVecs(isoVecs **ia, int levelmax);


int main()
{
	int n;
	for(n=0; n<8; n++)
	{
		enum_vect B = count(n);
// 	dans B[n] : nb d'arbres à 2n+1 noeuds
	
		cout << vect_str(B) << endl;
		mpz_class r;
		//r = rand(0, B[n]-1);
		for(int r=0; r<B[n]; r++)
		{
			node * A = tree_gen(2*n+1, B, r);
			cout << "l'arbre de taille " << 2*n+1 << " a " << HNonPlan(A) << " etiquetages possibles" << endl;
			//store(n, A, r);
		}
		cout << endl;
	}
	n = 2;
	enum_vect B = count(n);
	mpz_class r;
	node * A = tree_gen(2*n+1, B, r);
	cout << "l'arbre de taille " << 2*n+1 << " a " << HNonPlan(A) << " etiquetages possibles" << endl;
	
	
	
	
	return 0; 
}

int main_b()//experimentation HNonPlan
{
	
	//FILE *f = fopen("HNonPlanBits.time", "w");
	
	clock_t initial_time;
	clock_t final_time;
	double cpu_time=0.0;
	int n, round;
	for(n = 1; n < 1500; n++){
		//cout << n << endl;
		cpu_time = 0.0;
		enum_vect B = count(n);
		for(round = 0; round < 10; round++){
			
			mpz_class r = rand(0, B[n]-1);
			node * A = tree_gen(2*n+1, B, r);
		
			initial_time = clock();
			HNonPlanBits(A);
			final_time = clock();
			cpu_time += (((double)(final_time - initial_time))/CLOCKS_PER_SEC);
			
		}
		cpu_time /= 10;
		cout << "temps moyen pour n = " << 2*n+1 << " : " << cpu_time << endl;
		//fprintf(f, "%d %.6f\n", 2*n+1, cpu_time);
	}
	
	//fclose(f);
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

/////////////////////////////////////////////////

//Version 1//

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
////

//Version 2 avec sequence de bits//

treeValsBits * tv_initBits(mpz_class taille, int nbNS, mpz_class PTST, treeValsBits *fg, treeValsBits *fd)
{
	treeValsBits * tv;
	tv = new treeValsBits;
	tv->taille = taille;
	tv->nbNS = nbNS;
	tv->PTST = PTST;
	tv->fg = fg;
	tv->fd = fd;
	if(taille == 1){
		tv->bs.emplace_back(0);
		tv->bs.emplace_back(1);
	}else{
		tv->bs.emplace_back(0);
		tv->bs.insert(tv->bs.end(), tv->fg->bs.begin(), tv->fg->bs.end());
		tv->bs.emplace_back(0);
		tv->bs.insert(tv->bs.end(), tv->fd->bs.begin(), tv->fd->bs.end());
		tv->bs.emplace_back(1);
	}
	return tv;
}

void freeTreeValsBits (treeValsBits *T){
	if(T->fg !=NULL){
		freeTreeValsBits(T->fg);
		freeTreeValsBits(T->fd);
	}
	delete T;
}

mpz_class HNonPlanBits(node *tree){
	treeValsBits *rootTV = HNonPlanBitsAux(tree);
	mpz_class fact = factorial(rootTV->taille);
	int powNS = (int) pow(2, rootTV->nbNS);
	mpz_class powMPZ(powNS);
	mpz_class PTST = rootTV->PTST;
	freeTreeValsBits(rootTV);
	return fact / PTST / powMPZ;
}

treeValsBits *HNonPlanBitsAux(node *tree){
	if(isLeaf(tree)){
		return tv_initBits(1, 0, 1, NULL, NULL);
	}else{
		treeValsBits *fg = HNonPlanBitsAux(tree->fg);
		treeValsBits *fd = HNonPlanBitsAux(tree->fd);
		int nbNS;
		if(equalsTreeValsBits(fg, fd)){
			nbNS = fg->nbNS + fd->nbNS + 1;
		}else{
			nbNS = fg->nbNS + fd->nbNS;
		}
		mpz_class taille = fg->taille + fd->taille + 1;
		return tv_initBits(taille, nbNS, fg->PTST * fd->PTST * taille, fg, fd);
	}
}

bool equalsTreeValsBits(treeValsBits *fg, treeValsBits *fd){
	if(fg->taille != fd->taille){
		return false;
	}else{
		if(fg->taille == 1){
			return true;
		}else{
			return fg->bs == fd->bs;//on test seulement l'égalité entre les deux séquences de bits
		}
	}
}

////

mpz_class factorial(mpz_class n)
{
	if(n <= 0) return 1;
	mpz_class res(n);
    while(n-- > 1) res *= n;
    return res;
}

//iso
/*
isoVecs *init_isoVecs(){
	isoVecs *v = new isoVecs;
	v->maxLevel = 0;
	return v;
}

void resizeIsoVecs(isoVecs **ia, int levelmax){
	int previousMax = ia->maxLevel;
	if(previousMax == 0){
		ia->vecs = (int ***)malloc(sizeof(int**)*levelmax);
		ia->vecsSize = (int *)malloc(sizeof(int)*levelmax);
		ia->maxNodeIsoValueByLevel = (int *)malloc(sizeof(int)*levelmax);
		
		int i;
		for(i = 0; i <= levelmax; i++){
			ia->vecsSize[i] = 0;
			ia->maxNodeIsoValueByLevel[i] = 0;
		}
		ia->maxLevel = levelmax;
		
	}else if(previousMax < levelmax){
		ia->vecs = (int ***)realloc(ia->vecs, sizeof(int **)*levelmax);
		ia->vecsSize = (int *)realloc(ia->vecsSize, sizeof(int)*levelmax);
		ia->maxNodeIsoValueByLevel = (int *)realloc(ia->maxNodeIsoValueByLevel, sizeof(int)*levelmax);
		
		int i;
		for(i = previousMax+1; i <= levelmax; i++){
			ia->vecsSize[i] = 0;
			ia->maxNodeIsoValueByLevel[i] = 0;
		}
		ia->maxLevel = levelmax;
	}
}

treeIso *ti_init(mpz_class taille, int nbNS, mpz_class PTST, unsigned int level, treeIso *fg, treeIso *fd, treeIso *root){
	treeIso * t;
	t = new treeIso;
	t->taille = taille;
	t->nbNS = nbNS;
	t->PTST = PTST;
	t->root = root;
	t->fg = fg;
	t->fd = fd;
	t->level = level;
	return t;
}

void freeTreeIso(treeIso *T){
	if(T->fg !=NULL){
		freeTreeValsBits(T->fg);
		freeTreeValsBits(T->fd);
	}
	delete T;
}

mpz_class HNonPlanIso(node *tree){
	treeIso *rootTV = HNonPlanIsoAux(tree, NULL, 1);
	mpz_class fact = factorial(rootTV->taille);
	int powNS = (int) pow(2, rootTV->nbNS);
	mpz_class powMPZ(powNS);
	mpz_class PTST = rootTV->PTST;
	freeTreeIso(rootTV);
	return fact / PTST / powMPZ;
}

treeIso *HNonPlanIsoAux(node *tree, node *root, unsigned int level){
	if(isLeaf(tree)){
		treeIso *leaf = ti_init(1, 0, 1, level, NULL, NULL, root);
		leaf->nodeIsoValue = 0;
		return leaf;
	}else{
		treeIso *fg = HNonPlanIsoAux(tree->fg, tree, level+1);
		treeIso *fd = HNonPlanIsoAux(tree->fd, tree, level+1);
		mpz_class taille = fg->taille + fd->taille +1;
		mpz_class PTST = fg->PTST * fd->PTST * taille;
		//int isoV = getNodeIsoValue(isoVecs, fg->nodeIsoValue, fd->nodeIsoValue);
		int nbNS;
		if(equalsTreeIso(fg, fd)){
			nbNS = fg->nbNS + fd->nbNS + 1;
		}else{
			nbNS = fg->nbNS + fd->nbNS;
		}
		//TODO
		//treeIso n = ti_init(taille, nbNS, PTST, level, fg, fd, root);
		//n->nodeIsoValue = isoV; 
	}
}

bool equalsTreeIso(treeIso *fg, treeIso *fd){
	return fg->nodeIsoValue == fd->nodeIsoValue;
}

int getNodeIsoValue(isoVecs **iv, int level, int vfg, int vfd){

	
	
}
*/
///

//generation aleatoire uniforme
void sousCalcul(tree* t, tree** M, int weight, int* label){
   if(t == NULL)
	return;
   // count number
   int count = 0;
   for(int i=0; i< weight; i++)
	if( M[i] != NULL)
	  count++;
   // random r
   //RAND_MAX = (count-1);
   int r = std::rand();
   // found the r-th tree* and put on his label the value of *label
   int nb=0;
   tree* v;
   for(int i=0; i< weight; i++){
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
   v->label = *label;
   *label++;
   // construct M and recursive call
    int firstPlaceOfV = -1;
    // supprimer v
    for(int i=0; i < weight; i++)
	if(M[i]==v){
	  if(firstPlaceOfV == -1)
	    firstPlaceOfV = i;
	  M[i]==NULL;
	}
   // constuire M et appel recursive
   if( t-> fg!= NULL & t-> fd!=NULL){
	   for(int i=firstPlaceOfV; i < t->fg->weight; i++){
	     M[i]=t->fg;
		}
	   for(int i = firstPlaceOfV+t->fg->weight; i < (firstPlaceOfV+t->fg->weight+t->fd->weight); i++){
	     M[i]=t->fd;
		}

	   // call
	    sousCalcul(v->fg, M, weight, label);
	    sousCalcul(v->fd, M, weight, label);
}else if(t->fg!=NULL){
	for(int i=firstPlaceOfV; i < t->fg->weight; i++){
	     M[i]=t->fg;
	}
        sousCalcul(v->fg, M, weight, label);
  }else{
	for(int i = firstPlaceOfV; i < (firstPlaceOfV+t->fd->weight); i++){
	     M[i]=t->fd;
	}
	sousCalcul(v->fd, M, weight, label);
  }
}


/* random labeling
 * t : struct containing fg, fd , weight and label undefined of each node
 */
tree* randomLabeling(tree* t){ 
 if(t== NULL || t->weight < 1)
  return t;

 // init
 int* label = (int*) malloc(sizeof(int));
 *label = 1;
 tree** M = (tree**) malloc(sizeof(tree*)*t->weight);
 for(int i=0; i < t->weight; i++)
  M[i] = t;
  std::srand(std::time(nullptr));
 // recursive call
 sousCalcul(t, M, t->weight, label);
 free(M);
 free(label);
 return t;
}














