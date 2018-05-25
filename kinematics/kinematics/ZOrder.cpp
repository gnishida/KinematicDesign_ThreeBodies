#include "ZOrder.h"
#include "JointConnector.h"

namespace kinematics {

	std::vector<std::vector<int>> ZOrder::zorderConnectors(const std::vector<JointConnector>& connectors) {
		// type 1 constraint
		// depth of the first element has to be less than the depth of the second element
		std::vector<std::pair<int, int>> ct1;

		// type 2 constraint
		// depth of the first element has to be less than or greater than the second element
		std::vector<std::pair<int, std::vector<int>>> ct2;

		// type 3 constraint
		// depths of the first and second elements have to be different
		std::vector<std::pair<int, int>> ct3;

		for (int i = 0; i < connectors.size(); i++) {
			for (auto it = connectors[i].collisions1.begin(); it != connectors[i].collisions1.end(); it++) {
				int j = it.key();

				ct1.push_back({ j, i });
			}

			for (auto it = connectors[i].collisions2.begin(); it != connectors[i].collisions2.end(); it++) {
				//not_betweens.push_back(std::make_pair(i, it.value()));

				ct2.push_back({ i, it.value() });
			}

			for (auto it = connectors[i].neighbors.begin(); it != connectors[i].neighbors.end(); it++) {
				int j = it.key();

				ct3.push_back({ i, j });
			}
		}

		try {
			for (int max_depth = 1; max_depth < 10; max_depth++) {
				// create a model
				GecodeModel* model = new GecodeModel(connectors.size(), 1, max_depth, ct1, ct2, ct3);

				// create a search engine of DFS
				Gecode::BAB<GecodeModel> engine(model);

				// we can delete model since the search engine takes a clone of the model
				delete model;

				// get all solutions
				// the solution is actually a model.
				GecodeModel* sol = engine.next();
				if (sol) {
					std::vector<std::vector<int>> ans(3);
					int k = 0;
					for (int i = 0; i < connectors.size(); i++) {
						ans[connectors[i].type].push_back(sol->I[i].val());
					}
					delete sol;
					return ans;
				}
			}
		}
		catch (Gecode::Exception e) {
		}

		throw "No valid zorder was found.";
	}

}