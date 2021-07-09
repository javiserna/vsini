# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'rotation.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_rotation(object):
    def setupUi(self, rotation):
        rotation.setObjectName("rotation")
        rotation.resize(305, 299)
        self.centralwidget = QtWidgets.QWidget(rotation)
        self.centralwidget.setObjectName("centralwidget")
        self.buttom1 = QtWidgets.QPushButton(self.centralwidget)
        self.buttom1.setGeometry(QtCore.QRect(20, 10, 111, 31))
        self.buttom1.setStyleSheet("background-color: rgb(114, 159, 207);")
        self.buttom1.setObjectName("buttom1")
        self.buttom2 = QtWidgets.QPushButton(self.centralwidget)
        self.buttom2.setGeometry(QtCore.QRect(160, 10, 131, 41))
        self.buttom2.setObjectName("buttom2")
        self.buttom3 = QtWidgets.QPushButton(self.centralwidget)
        self.buttom3.setGeometry(QtCore.QRect(220, 110, 71, 61))
        self.buttom3.setStyleSheet("background-color: rgb(138, 226, 52);")
        self.buttom3.setObjectName("buttom3")
        self.linecenter = QtWidgets.QLineEdit(self.centralwidget)
        self.linecenter.setGeometry(QtCore.QRect(20, 50, 121, 25))
        self.linecenter.setObjectName("linecenter")
        self.linewidth = QtWidgets.QLineEdit(self.centralwidget)
        self.linewidth.setGeometry(QtCore.QRect(20, 80, 121, 25))
        self.linewidth.setObjectName("linewidth")
        self.limbdarkening = QtWidgets.QLineEdit(self.centralwidget)
        self.limbdarkening.setGeometry(QtCore.QRect(20, 140, 191, 25))
        self.limbdarkening.setObjectName("limbdarkening")
        self.resolution = QtWidgets.QLineEdit(self.centralwidget)
        self.resolution.setGeometry(QtCore.QRect(20, 110, 191, 25))
        self.resolution.setObjectName("resolution")
        self.resultado = QtWidgets.QTextEdit(self.centralwidget)
        self.resultado.setGeometry(QtCore.QRect(10, 200, 281, 41))
        self.resultado.setObjectName("resultado")
        self.buttom4 = QtWidgets.QPushButton(self.centralwidget)
        self.buttom4.setGeometry(QtCore.QRect(160, 60, 131, 41))
        self.buttom4.setObjectName("buttom4")
        self.progressBar = QtWidgets.QProgressBar(self.centralwidget)
        self.progressBar.setGeometry(QtCore.QRect(90, 170, 118, 23))
        self.progressBar.setProperty("value", 24)
        self.progressBar.setObjectName("progressBar")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(140, 240, 161, 20))
        self.label.setObjectName("label")
        rotation.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(rotation)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 305, 22))
        self.menubar.setObjectName("menubar")
        rotation.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(rotation)
        self.statusbar.setObjectName("statusbar")
        rotation.setStatusBar(self.statusbar)

        self.retranslateUi(rotation)
        QtCore.QMetaObject.connectSlotsByName(rotation)

    def retranslateUi(self, rotation):
        _translate = QtCore.QCoreApplication.translate
        rotation.setWindowTitle(_translate("rotation", "Vsini by Fourier\'s Method"))
        self.buttom1.setText(_translate("rotation", "Load"))
        self.buttom2.setText(_translate("rotation", "Plot Spectrum"))
        self.buttom3.setText(_translate("rotation", "Run"))
        self.linecenter.setText(_translate("rotation", ""))
        self.linecenter.setPlaceholderText(_translate("rotation", "Line Center"))
        self.linewidth.setText(_translate("rotation", ""))
        self.linewidth.setPlaceholderText(_translate("rotation", "Line width"))
        self.limbdarkening.setText(_translate("rotation", "0.6"))
        self.limbdarkening.setPlaceholderText(_translate("rotation", "Limb-darkening coefficient"))
        self.resolution.setText(_translate("rotation", "22500"))
        self.resolution.setPlaceholderText(_translate("rotation", "Spectral Resolution (R)"))
        self.resultado.setHtml(_translate("rotation", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p align=\"center\" style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.buttom4.setText(_translate("rotation", "Plot Line"))
        self.label.setText(_translate("rotation", "<html><head/><body><p><span style=\" font-size:8pt; color:#888a85;\">© Javier Serna / IA-UNAM/ 2020</span></p></body></html>"))

