#include "tree.h"

hl::Tree::Tree(int deg, int val, vector<Tree> subT) : degree(deg), value(val), subTrees(subT) {}

bool hl::Tree::isIncreasing(int val, hl::Tree tree){
	if(tree.degree == 0){//Leaf
		return val < tree.value;
	}else{
		if(! (val < tree.value)){
			return false;
		}
		for(int i = 0; i < tree.degree; i++){
			if(! hl::Tree::isIncreasing(tree.value, tree.subTrees[i]))
				return false;
		}
		return true;
	}
}

bool hl::Tree::isIncreasingTree(hl::Tree tree){//return true if tree has stricly increasing subTrees
	return hl::Tree::isIncreasing(-1, tree);
}
