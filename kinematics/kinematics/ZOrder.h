#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <unordered_map>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include <gecode/search.hh>

namespace kinematics {

	class JointConnector;

	class GecodeModel : public Gecode::Space {
	public:
		Gecode::IntVarArray I;

	public:
		GecodeModel(int num_variables, int min_value, int max_value, std::vector<std::pair<int, int>> ct1, std::vector<std::pair<int, std::vector<int>>> ct2, std::vector<std::pair<int, int>> ct3) : I(*this, num_variables, min_value, max_value) {
			// type 1 constraint
			// depth of the first element has to be less than the depth of the second element
			for (int i = 0; i < ct1.size(); i++) {
				rel(*this, I[ct1[i].first] < I[ct1[i].second]);
			}

			// type 2 constraint
			// depth of the first element has to be less than or greater than the second element
			for (int i = 0; i < ct2.size(); i++) {
				Gecode::IntVarArray ar(*this, ct2[i].second.size());
				for (int j = 0; j < ct2[i].second.size(); j++) {
					ar[j] = I[ct2[i].second[j]];
				}
				rel(*this, I[ct2[i].first] < Gecode::min(ar) || I[ct2[i].first] > Gecode::max(ar));
			}

			// type 3 constraint
			// depths of the first and second elements have to be different
			for (int i = 0; i < ct3.size(); i++) {
				rel(*this, I[ct3[i].first] != I[ct3[i].second]);
			}

			// algorithm
			// INT_VAR_SIZE_MIN(): the one with smallest domain size first
			// in this example, all have the same range [0, 9], though.
			// INT_VAL_MIN(): the smallest value of the selected variable first
			// i.e., try the values, 0, ..., 9, in this order.
			branch(*this, I, Gecode::INT_VAR_SIZE_MIN(), Gecode::INT_VAL_MIN());
		}

		// copy constructor
		GecodeModel(bool share, GecodeModel& s) : Space(share, s) {
			I.update(*this, share, s.I);
		}

		virtual Space* copy(bool share) {
			return new GecodeModel(share, *this);
		}

		virtual void constrain(const Space& _b) {
			const GecodeModel& b = static_cast<const GecodeModel&>(_b);

			rel(*this, Gecode::max(I) > Gecode::max(b.I));
		}

		void print() const {
			std::cout << I << std::endl;
		}
	};

	class ZOrder {
	protected:
		ZOrder() {}

	public:
		static std::vector<std::vector<int>> zorderConnectors(const std::vector<JointConnector>& connectors);
	};
}