#include <vector>

namespace hl{

	using namespace std;
	class Tree {
		private:
			int degree;//number of subTrees (children)
			int value;//node's value
			vector<Tree> subTrees;
		public:
			Tree(int deg, int val, vector<Tree> subT);
			bool isIncreasing(int val, hl::Tree tree);
			bool isIncreasingTree(Tree t);
	};
}
