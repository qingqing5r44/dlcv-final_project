/********************************************************************************
** Form generated from reading UI file 'main_window.ui'
**
** Created by: Qt User Interface Compiler version 5.12.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAIN_WINDOW_H
#define UI_MAIN_WINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_Main
{
public:
    QLabel *label;
    QPushButton *btn_openpose;
    QPushButton *btn_smplifyx;
    QPushButton *btn_panoman;
    QPushButton *btn_upload;
    QLabel *label_photo;
    QLabel *label_openpose;
    QLabel *label_smplifyx;
    QLabel *label_panoman;

    void setupUi(QDialog *Main)
    {
        if (Main->objectName().isEmpty())
            Main->setObjectName(QString::fromUtf8("Main"));
        Main->resize(770, 723);
        label = new QLabel(Main);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(150, 120, 441, 421));
        label->setStyleSheet(QString::fromUtf8("background-color: rgb(255, 255, 255);"));
        btn_openpose = new QPushButton(Main);
        btn_openpose->setObjectName(QString::fromUtf8("btn_openpose"));
        btn_openpose->setGeometry(QRect(40, 610, 181, 51));
        btn_smplifyx = new QPushButton(Main);
        btn_smplifyx->setObjectName(QString::fromUtf8("btn_smplifyx"));
        btn_smplifyx->setGeometry(QRect(280, 610, 181, 51));
        btn_panoman = new QPushButton(Main);
        btn_panoman->setObjectName(QString::fromUtf8("btn_panoman"));
        btn_panoman->setGeometry(QRect(530, 610, 181, 51));
        btn_upload = new QPushButton(Main);
        btn_upload->setObjectName(QString::fromUtf8("btn_upload"));
        btn_upload->setGeometry(QRect(40, 60, 93, 28));
        label_photo = new QLabel(Main);
        label_photo->setObjectName(QString::fromUtf8("label_photo"));
        label_photo->setGeometry(QRect(170, 60, 441, 21));
        label_openpose = new QLabel(Main);
        label_openpose->setObjectName(QString::fromUtf8("label_openpose"));
        label_openpose->setGeometry(QRect(40, 680, 181, 16));
        label_smplifyx = new QLabel(Main);
        label_smplifyx->setObjectName(QString::fromUtf8("label_smplifyx"));
        label_smplifyx->setGeometry(QRect(280, 680, 181, 16));
        label_panoman = new QLabel(Main);
        label_panoman->setObjectName(QString::fromUtf8("label_panoman"));
        label_panoman->setGeometry(QRect(810, 680, 181, 16));

        retranslateUi(Main);
        QObject::connect(btn_openpose, SIGNAL(clicked()), Main, SLOT(LoadOpenpose()));
        QObject::connect(btn_panoman, SIGNAL(clicked()), Main, SLOT(LoadPANOMAN()));
        QObject::connect(btn_smplifyx, SIGNAL(clicked()), Main, SLOT(LoadSmplifyx()));
        QObject::connect(btn_upload, SIGNAL(clicked()), Main, SLOT(UploadPhoto()));

        QMetaObject::connectSlotsByName(Main);
    } // setupUi

    void retranslateUi(QDialog *Main)
    {
        Main->setWindowTitle(QApplication::translate("Main", "\345\237\272\344\272\216PANOMAN\345\217\202\346\225\260\345\214\226\346\250\241\345\236\213\347\232\204\344\270\211\347\273\264\344\272\272\344\275\223\351\207\215\345\273\272", nullptr));
        label->setText(QString());
        btn_openpose->setText(QApplication::translate("Main", "\346\217\220\345\217\226\345\205\263\351\224\256\347\202\271", nullptr));
        btn_smplifyx->setText(QApplication::translate("Main", "\344\274\260\350\256\241\345\217\202\346\225\260", nullptr));
        btn_panoman->setText(QApplication::translate("Main", "\351\207\215\345\273\272PANOMAN\346\250\241\345\236\213", nullptr));
        btn_upload->setText(QApplication::translate("Main", "\344\270\212\344\274\240\347\205\247\347\211\207", nullptr));
        label_photo->setText(QString());
        label_openpose->setText(QString());
        label_smplifyx->setText(QString());
        label_panoman->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class Main: public Ui_Main {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAIN_WINDOW_H
