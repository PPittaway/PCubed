'''
    Code for running the MainWindow of the PCubed GUI
    
    Created in Python 3.9.0
    
'''

from PyQt5 import QtCore, QtWidgets
import sys
import platformControl
import experimentMethod
import analysisHub

class mainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(mainWindow, self).__init__()
        global windowRun
        windowRun = False
        self.centralWidget = QtWidgets.QWidget(self)

        ####### Add Tab widget to the centralWidget of the mainWindow #######
        self.setCentralWidget(self.centralWidget)
        self.tab_main = QtWidgets.QTabWidget(self)
        self.tab_main.layout = QtWidgets.QGridLayout()
        self.tab_main.setStyleSheet("QTabBar::tab { height: 30px; width: 150px }")
        layout = QtWidgets.QGridLayout(self.centralWidget)
        layout.addWidget(self.tab_main)

        ####### Create and add custom widgets to the tab ########
        self.tab_1 = QtWidgets.QTabWidget()
        tab_1 = QtWidgets.QGridLayout(self.tab_1)
        self.methodHandler = experimentMethod.ExperimentMethod(self, main=self)
        tab_1.addWidget(self.methodHandler)
        
        self.tab_2 = QtWidgets.QTabWidget()
        tab_2 = QtWidgets.QGridLayout(self.tab_2)
        self.controller = platformControl.PlatformControl(self, main=self)
        tab_2.addWidget(self.controller)

        self.tab_3 = QtWidgets.QTabWidget()
        tab_3 = QtWidgets.QGridLayout(self.tab_3)
        self.analysisHub = analysisHub.analysisHub(self, main=self)
        tab_3.addWidget(self.analysisHub)

        ####### Add tabs and titles to the master tab widget #######
        self.tab_main.addTab(self.tab_1, "Method builder")
        self.tab_main.addTab(self.tab_2, "Platform control")
        self.tab_main.addTab(self.tab_3, "Analysis hub")

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("Platform Controller")
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("Platform Controller", "Platform Controller"))

if __name__ == "__main__":
    import sys
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps)
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = mainWindow()
    ui = mainWindow()
    ui.setupUi(MainWindow)
    MainWindow.setMinimumSize(1900, 1000)
    MainWindow.showMaximized()
    sys.exit(app.exec_())
