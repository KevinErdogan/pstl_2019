#include <vector>

namespace hl{

	using namespace std;
	class Tree {
		private:
			int degree;//number of subTrees (children)
			int value;//node's value
			vector<Tree> subTrees;
			typedef struct tv{
				int taille;
				int nbNS;
				unsigned int PTST;
				tv *fg;
				tv *fd;
			} treeVals;
		public:
			Tree(int deg, int val, vector<Tree> subT);
			bool isIncreasing(int val, hl::Tree tree);
			bool isIncreasingTree(Tree t);
			unsigned long long HNonPlan(hl::Tree tree);
			treeVals HNonPlanAux(hl::Tree tree);
			bool equalsTreeVals(treeVals fg, treeVals fd);

	};
}
