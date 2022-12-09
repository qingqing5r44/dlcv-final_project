#pragma once
#ifndef MAIN_H
#define MAIN_H
#include <QWidget>
#include <QDialog>
#include "ui_main_window.h"
#include <string.h>  //memset()
#include <fstream>
namespace Ui {
	class Main;
}
class Main : public QDialog
{
	Q_OBJECT

public:
	Main(QWidget* parent = Q_NULLPTR);
	~Main();
	int argc;
	char** argv;

public slots:
	void LoadPANOMAN();
	void LoadOpenpose();
	void LoadSmplifyx();
	void UploadPhoto();

private:
	Ui::Main* ui;
	std::string photo_path;
};
void testHumanModel7(int argc, char** argv);
#endif