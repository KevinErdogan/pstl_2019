
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

// compiler avec g++ XXX.cpp -lgmpxx -lgmp -o XXX

using namespace std;


typedef vector <mpz_class> enum_vect;
struct node
{
	node * fg;
	node * fd;
};

typedef struct tv{
				int taille;
				int nbNS;
				unsigned long long PTST;
				tv *fg;
				tv *fd;
} treeVals;

enum_vect count(int n);
string vect_str(enum_vect G);
node * leaf_init();
node * binaire(node * T, node * U);
mpz_class rand(mpz_class min, mpz_class max);
node * tree_gen(int n, enum_vect G, mpz_class r);
string dot_from_tree(clock_t a, node * A);
void store(int n, node * A, mpz_class r);

////////////:

treeVals * tv_init(int taille, int nbNS, unsigned int PTST, treeVals *fg, treeVals *fd);
bool isLeaf(node *n);
unsigned long long HNonPlan(node *tree);
treeVals * HNonPlanAux(node *tree);
bool equalsTreeVals(treeVals *fg, treeVals *fd);



int main()
{
	int n;
	for(n=0; n<8; n++)
	{
		enum_vect B = count(n);
// 	dans B[n] : nb d'arbres à 2n+1 noeuds
	
		cout << vect_str(B) << "\n";
		mpz_class r;
		//r = rand(0, B[n]-1);
		for(int r=0; r<B[n]; r++)
		{
			node * A = tree_gen(2*n+1, B, r);
			printf("l'arbre de taille %d a %llu etiquetages possibles\n", 2*n+1, HNonPlan(A));
			//store(n, A, r);
		}
		cout << endl;
	}
	/*n = 2;
	enum_vect B = count(n);
	mpz_class r;
	node * A = tree_gen(2*n+1, B, r);
	printf("l'arbre de taille %d a %llu etiquetages possibles\n", 2*n+1, HNonPlan(A));*/
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

/////////////////////////////////////////////////

treeVals * tv_init(int taille, int nbNS, unsigned int PTST, treeVals *fg, treeVals *fd)
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

bool isLeaf(node *n){
	return (n->fg == NULL && n->fd == NULL);
}

unsigned long long HNonPlan(node *tree){
	treeVals *rootTV = HNonPlanAux(tree);
	unsigned long long fact = 1;
	for(int i = 1; i <= rootTV->taille; i++){
		fact *= i;
	}
	//printf("taille = %d, fact = %llu, PTST = %llu, nbNS = %d\n", rootTV->taille, fact, rootTV->PTST, rootTV->nbNS);
	return fact / rootTV->PTST / (int) pow(2, rootTV->nbNS);
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
		int taille = fg->taille + fd->taille + 1;
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
