#include "Link.h"
#include "Joint.h"
#include <QPolygon>

namespace kinematics {

	Link::Link(int id) {
		this->id = id;
		this->angle = 0;
		this->driver = false;
		this->actual_link = true;
		this->z = 0;
	}

	Link::Link(int id, bool driver, bool actual_link, int z) {
		this->id = id;
		this->angle = 0;
		this->driver = driver;
		this->actual_link = actual_link;
		this->z = z;
	}

	bool Link::isDetermined() {
		int count = 0;
		for (int i = 0; i < joints.size(); ++i) {
			if (joints[i]->determined) count++;
		}

		if (count >= 2) return true;
		else return false;
	}

	bool Link::isGrounded() {
		int count = 0;
		for (int i = 0; i < joints.size(); ++i) {
			if (joints[i]->ground) count++;
		}

		if (count >= 2) return true;
		else return false;
	}

	void Link::addJoint(boost::shared_ptr<Joint> joint) {
		joints.push_back(joint);
		if (joints.size() == 2) {
			// initialize the link angle
			angle = atan2(joints[1]->pos.y - joints[0]->pos.y, joints[1]->pos.x - joints[0]->pos.x);
		}
	}

	void Link::rotate(const glm::dvec2& rotation_center, double angle) {
		for (int i = 0; i < joints.size(); ++i) {
			joints[i]->rotate(rotation_center, angle);
		}
		this->angle += angle;
	}

	double Link::getLength(int joint_id1, int joint_id2) {
		return glm::length(original_shape[joint_id1] - original_shape[joint_id2]);
	}

	glm::dmat3x2 Link::getTransformMatrix() {
		std::vector<glm::dvec2> orig_p;
		std::vector<glm::dvec2> p;
		for (int i = 0; i < joints.size(); ++i) {
			if (joints[i]->determined) {
				p.push_back(joints[i]->next_pos);
				orig_p.push_back(original_shape[joints[i]->id]);
			}
		}

		if (p.size() < 2) throw "Undetermined";

		double num = (p[1].y - p[0].y) * (orig_p[1].x - orig_p[0].x) - (p[1].x - p[0].x) * (orig_p[1].y - orig_p[0].y);
		double den = (p[1].x - p[0].x) * (orig_p[1].x - orig_p[0].x) + (p[1].y - p[0].y) * (orig_p[1].y - orig_p[0].y);
		double theta = atan2(num, den);
		glm::dmat3x2 ret;
		ret[0][0] = cos(theta);
		ret[0][1] = sin(theta);
		ret[1][0] = -sin(theta);
		ret[1][1] = cos(theta);
		ret[2][0] = -cos(theta) * orig_p[0].x + sin(theta) * orig_p[0].y + p[0].x;
		ret[2][1] = -sin(theta) * orig_p[0].x - cos(theta) * orig_p[0].y + p[0].y;

		return ret;
	}

	glm::dvec2 Link::transformByDeterminedJoints(int joint_id) {
		glm::dmat3x2 mat = getTransformMatrix();

		return mat * glm::dvec3(original_shape[joint_id], 1);
	}

	void Link::draw(QPainter& painter, const QPointF& origin, float scale) {
		if (!actual_link) return;

		painter.save();

		if (driver) {
			painter.setPen(QPen(QColor(0, 0, 0), 3));
		}
		else {
			painter.setPen(QPen(QColor(90, 90, 90), 3));
		}
		if (joints.size() == 2) {
			painter.drawLine(origin.x() + joints[0]->pos.x * scale, origin.y() - joints[0]->pos.y * scale, origin.x() + joints[1]->pos.x * scale, origin.y() - joints[1]->pos.y * scale);
		}
		else {
			painter.setBrush(QBrush(QColor(192, 192, 192, 64)));
			QPolygon polygon;
			for (int i = 0; i < joints.size(); ++i) {
				polygon.append(QPoint(origin.x() + joints[i]->pos.x * scale, origin.y() - joints[i]->pos.y * scale));
			}
			painter.drawPolygon(polygon);
		}

		painter.restore();
	}

}