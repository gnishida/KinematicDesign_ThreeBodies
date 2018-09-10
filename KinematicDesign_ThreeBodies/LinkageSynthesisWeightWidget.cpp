#include "LinkageSynthesisWeightWidget.h"
#include "MainWindow.h"
#include "GLWidget3D.h"

LinkageSynthesisWeightWidget::LinkageSynthesisWeightWidget(MainWindow *parent, const std::vector<double>& weights) : QDockWidget((QWidget*)parent) {
	mainWin = parent;
	ui.setupUi(this);

	ui.horizontalSliderPositionError->setRange(1, 10);
	ui.horizontalSliderLinkageLocation->setRange(1, 10);
	ui.horizontalSliderTrajectory->setRange(1, 10);
	ui.horizontalSliderSize->setRange(1, 10);

	connect(ui.horizontalSliderPositionError, SIGNAL(valueChanged(int)), this, SLOT(onValueChanged(int)));
	connect(ui.horizontalSliderLinkageLocation, SIGNAL(valueChanged(int)), this, SLOT(onValueChanged(int)));
	connect(ui.horizontalSliderTrajectory, SIGNAL(valueChanged(int)), this, SLOT(onValueChanged(int)));
	connect(ui.horizontalSliderSize, SIGNAL(valueChanged(int)), this, SLOT(onValueChanged(int)));
	connect(ui.pushButtonApply, SIGNAL(clicked()), this, SLOT(onApply()));

	ui.horizontalSliderPositionError->setValue(weights[0]);
	ui.horizontalSliderLinkageLocation->setValue(weights[1]);
	ui.horizontalSliderTrajectory->setValue(weights[2]);
	ui.horizontalSliderSize->setValue(weights[3]);
}

LinkageSynthesisWeightWidget::~LinkageSynthesisWeightWidget() {
}

void LinkageSynthesisWeightWidget::onValueChanged(int) {
	ui.lineEditPositionError->setText(QString::number(ui.horizontalSliderPositionError->value()));
	ui.lineEditLinkageLocation->setText(QString::number(ui.horizontalSliderLinkageLocation->value()));
	ui.lineEditTrajectory->setText(QString::number(ui.horizontalSliderTrajectory->value()));
	ui.lineEditSize->setText(QString::number(ui.horizontalSliderSize->value()));

	mainWin->glWidget->weights = {
		ui.lineEditPositionError->text().toDouble(),
		ui.lineEditLinkageLocation->text().toDouble(),
		ui.lineEditTrajectory->text().toDouble(),
		ui.lineEditSize->text().toDouble()
	};
}

void LinkageSynthesisWeightWidget::onApply() {
	mainWin->glWidget->updateSolutions(GLWidget3D::LINKAGE_WATTI);
	mainWin->glWidget->setFocus();
}