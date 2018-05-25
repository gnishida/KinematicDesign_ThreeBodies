#include "LeastSquareSolver.h"
#include "KinematicUtils.h"

namespace kinematics {

	SolverForLink::SolverForLink(const std::vector<glm::dmat3x3>& poses) {
		this->poses = poses;
	}

	double SolverForLink::operator() (const column_vector& arg) const {
		glm::dvec2 A0(arg(0, 0), arg(1, 0));
		glm::dvec2 a(arg(2, 0), arg(3, 0));

		glm::dvec2 A1(poses[0] * glm::dvec3(a, 1));
		double l1_squared = glm::length(A1 - A0);
		l1_squared = l1_squared * l1_squared;

		double ans = 0.0;
		for (int i = 1; i < poses.size(); i++) {
			glm::dvec2 A(poses[i] * glm::dvec3(a, 1));
			double l_squared = glm::length(A - A0);
			l_squared = l_squared * l_squared;
			ans += (l_squared - l1_squared) * (l_squared - l1_squared);
		}

		return ans;
	}

	SolverForSlider::SolverForSlider(const std::vector<glm::dmat3x3>& poses) {
		this->poses = poses;
	}

	double SolverForSlider::operator() (const column_vector& arg) const {
		glm::dvec2 a(arg(0, 0), arg(1, 0));

		glm::dvec2 A1(poses[0] * glm::dvec3(a, 1));
		glm::dvec2 A2(poses[1] * glm::dvec3(a, 1));

		glm::dvec2 v1 = A2 - A1;
		v1 /= glm::length(v1);

		double ans = 0.0;
		for (int i = 2; i < poses.size(); i++) {
			glm::dvec2 A(poses[i] * glm::dvec3(a, 1));
			glm::dvec2 v = A - A1;
			v /= glm::length(v);

			ans += abs(crossProduct(v1, v));
		}

		return ans;
	}

	SolverForWattI::SolverForWattI(const std::vector<std::vector<glm::dmat3x3>>& poses) {
		pose_params.resize(poses.size());
		for (int i = 0; i < poses.size(); i++) {
			pose_params[i].resize(poses[i].size() - 1, std::vector<double>(3));

			for (int j = 1; j < poses[i].size(); j++) {
				glm::dmat3x3 D = poses[i][j] * glm::inverse(poses[i][0]);
				pose_params[i][j - 1][0] = std::atan2(D[0][1], D[0][0]);
				pose_params[i][j - 1][1] = D[2][0];
				pose_params[i][j - 1][2] = D[2][1];
			}
		}
	}

	double SolverForWattI::operator() (const column_vector& arg) const {
		glm::dvec2 P0(arg(0, 0), arg(1, 0));
		glm::dvec2 P1(arg(2, 0), arg(3, 0));
		glm::dvec2 P2(arg(4, 0), arg(5, 0));
		glm::dvec2 P3(arg(6, 0), arg(7, 0));
		glm::dvec2 P4(arg(8, 0), arg(9, 0));
		glm::dvec2 P5(arg(10, 0), arg(11, 0));
		glm::dvec2 P6(arg(12, 0), arg(13, 0));
		glm::dvec2 P7(arg(14, 0), arg(15, 0));
		glm::dvec2 P8(arg(16, 0), arg(17, 0));
		glm::dvec2 P9(arg(18, 0), arg(19, 0));

		double ans = 0.0;
		for (int i = 0; i < pose_params[0].size(); i++) {
			double theta = pose_params[0][i][0];
			double u = pose_params[0][i][1];
			double v = pose_params[0][i][2];
			double phi = pose_params[1][i][0];
			double s = pose_params[1][i][1];
			double t = pose_params[1][i][2];
			double psi = pose_params[2][i][0];
			double q = pose_params[2][i][1];
			double r = pose_params[2][i][2];

			double z0 = 2 * (P2.x * P0.x + P2.y * P0.y) * (1 - std::cos(theta))
				+ 2 * (P2.y * P0.x - P2.x * P0.y) * std::sin(theta)
				- 2 * u * P0.x
				- 2 * v * P0.y
				+ 2 * P2.x * (u * std::cos(theta) + v * std::sin(theta))
				+ 2 * P2.y * (-u * std::sin(theta) + v * std::cos(theta))
				+ u * u
				+ v * v;

			double z1 = 2 * (P3.x * P1.x + P3.y * P1.y) * (1 - std::cos(theta))
				+ 2 * (P3.y * P1.x - P3.x * P1.y) * std::sin(theta)
				- 2 * u * P1.x
				- 2 * v * P1.y
				+ 2 * P3.x * (u * std::cos(theta) + v * std::sin(theta))
				+ 2 * P3.y * (-u * std::sin(theta) + v * std::cos(theta))
				+ u * u
				+ v * v;

			double z2 = 2 * (P4.x * P1.x + P4.y * P1.y) * (1 - std::cos(phi))
				+ 2 * (P4.y * P1.x - P4.x * P1.y) * std::sin(phi)
				- 2 * s * P1.x
				- 2 * t * P1.y
				+ 2 * P4.x * (s * std::cos(phi) + t * std::sin(phi))
				+ 2 * P4.y * (-s * std::sin(phi) + t * std::cos(phi))
				+ s * s
				+ t * t;

			double z3 = 2 * (P4.x * P3.x + P4.y * P3.y) * (1 - std::sin(phi) * std::sin(theta) - std::cos(phi) * std::cos(theta))
				+ 2 * (P4.x * P3.y - P4.y * P3.x) * (std::cos(phi) * std::sin(theta) - std::sin(phi) * std::cos(theta))
				+ 2 * P4.x * ((t - v) * std::sin(phi) + (s - u) * std::cos(phi))
				+ 2 * P4.y * ((-s + u) * std::sin(phi) + (t - v) * std::cos(phi))
				+ 2 * P3.x * ((-t + v) * std::sin(theta) + (-s + u) * std::cos(theta))
				+ 2 * P3.y * ((s - u) * std::sin(theta) + (-t + v) * std::cos(theta))
				+ (s - u) * (s - u)
				+ (t - v) * (t - v);

			double z4 = 2 * (P6.x * P5.x + P6.y * P5.y) * (1 - std::sin(phi) * std::sin(theta) - std::cos(phi) * std::cos(theta))
				+ 2 * (P6.x * P5.y - P6.y * P5.x) * (std::cos(phi) * std::sin(theta) - std::sin(phi) * std::cos(theta))
				+ 2 * P6.x * ((t - v) * std::sin(phi) + (s - u) * std::cos(phi))
				+ 2 * P6.y * ((-s + u) * std::sin(phi) + (t - v) * std::cos(phi))
				+ 2 * P5.x * ((-t + v) * std::sin(theta) + (-s + u) * std::cos(theta))
				+ 2 * P5.y * ((s - u) * std::sin(theta) + (-t + v) * std::cos(theta))
				+ (s - u) * (s - u)
				+ (t - v) * (t - v);

			double z5 = 2 * (P8.x * P5.x + P8.y * P5.y) * (1 - std::sin(psi) * std::sin(theta) - std::cos(psi) * std::cos(theta))
				+ 2 * (P8.x * P5.y - P8.y * P5.x) * (std::cos(psi) * std::sin(theta) - std::sin(psi) * std::cos(theta))
				+ 2 * P8.x * ((r - v) * std::sin(psi) + (q - u) * std::cos(psi))
				+ 2 * P8.y * ((-q + u) * std::sin(psi) + (r - v) * std::cos(psi))
				+ 2 * P5.x * ((-r + v) * std::sin(theta) + (-q + u) * std::cos(theta))
				+ 2 * P5.y * ((q - u) * std::sin(theta) + (-r + v) * std::cos(theta))
				+ (q - u) * (q - u)
				+ (r - v) * (r - v);

			double z6 = 2 * (P8.x * P6.x + P8.y * P6.y) * (1 - std::sin(psi) * std::sin(phi) - std::cos(psi) * std::cos(phi))
				+ 2 * (P8.x * P6.y - P8.y * P6.x) * (std::cos(psi) * std::sin(phi) - std::sin(psi) * std::cos(phi))
				+ 2 * P8.x * ((r - t) * std::sin(psi) + (q - s) * std::cos(psi))
				+ 2 * P8.y * ((-q + s) * std::sin(psi) + (r - t) * std::cos(psi))
				+ 2 * P6.x * ((-r + t) * std::sin(phi) + (-q + s) * std::cos(phi))
				+ 2 * P6.y * ((q - s) * std::sin(phi) + (-r + t) * std::cos(phi))
				+ (q - s) * (q - s)
				+ (r - t) * (r - t);

			double z7 = 2 * (P9.x * P7.x + P9.y * P7.y) * (1 - std::sin(psi) * std::sin(phi) - std::cos(psi) * std::cos(phi))
				+ 2 * (P9.x * P7.y - P9.y * P7.x) * (std::cos(psi) * std::sin(phi) - std::sin(psi) * std::cos(phi))
				+ 2 * P9.x * ((r - t) * std::sin(psi) + (q - s) * std::cos(psi))
				+ 2 * P9.y * ((-q + s) * std::sin(psi) + (r - t) * std::cos(psi))
				+ 2 * P7.x * ((-r + t) * std::sin(phi) + (-q + s) * std::cos(phi))
				+ 2 * P7.y * ((q - s) * std::sin(phi) + (-r + t) * std::cos(phi))
				+ (q - s) * (q - s)
				+ (r - t) * (r - t);

			ans += z0 * z0 + z1 * z1 + z2 * z2 + z3 * z3 + z4 * z4 + z5 * z5 + z6 * z6 + z7 * z7;
		}

		return ans;
	}
	
}