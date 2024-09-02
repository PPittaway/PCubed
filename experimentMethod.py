from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtWidgets import QMessageBox
import pandas as pd
import seedAmountScreen
import surfactantScreen
import nFeedsScreen
import monomerScreen
import numpy as np
import threading
import datetime
import latinHypercubeSampling
import time
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from types import SimpleNamespace
import flowrateCalculator
import DoEbuilder
import conventionalEP
import RAFT_EP
import glassTransitionPredictor
import glob
import optimiser

class ExperimentMethod(QtWidgets.QWidget):
    def __init__(self, parent, main):
        super(ExperimentMethod, self).__init__(parent)

        self.main = main

        # Initialise each possible experimental method - each method calculates parameters slightly differently
        self.surfactantScreener = surfactantScreen.surfactantScreener(self, main=self.main)
        self.seedAmountScreener = seedAmountScreen.seedAmountScreener(self, main=self.main)
        self.nFeedsScreener = nFeedsScreen.nFeedsScreener(self, main=self.main)
        self.monomerScreener = monomerScreen.monomerScreener(self, main=self.main)
        self.DoEHandler = DoEbuilder.DoEbuilder(self, main=self.main)
        self.conventionalEPHandler = conventionalEP.conventionalEP(self, main=self.main)
        self.RAFTEPHandler = RAFT_EP.RAFT_EP(self, main=self.main)
        
        # Set the master layout of the Method tab onto which all widgets are placed
        self.Layout = QtWidgets.QGridLayout()
        self.setLayout(self.Layout)
        
        ###################################### Method builder box ##############################################

        # Main groupBox for all input and calculated experimental parameters
        self.methodBuilderBox = QtWidgets.QGroupBox("Method builder")
        # self.methodBuilderBox.setMaximumSize(680, 910)
        self.methodBuilderBox.setFixedWidth(650)
        self.methodBuilderBox.setContentsMargins(10, 5, 10, 5)
        methodBuilderLayout = QtWidgets.QGridLayout()
        self.methodBuilderBox.setLayout(methodBuilderLayout)

        ################################## Polymerisation selector box #########################################

        # Create groupBox for selecting the type of polymerisation for the experiment and add grid layout
        self.polymerisationSelectorBox = QtWidgets.QGroupBox("Polymerisation selector")
        self.polymerisationSelectorBox.setMaximumSize(250, 100)
        self.polymerisationSelectorLayout = QtWidgets.QGridLayout()
        self.polymerisationSelectorBox.setLayout(self.polymerisationSelectorLayout)

        # Create a checkBox for each type of polymerisation available and setAutoExclusive so only one can be selected at a time
        self.conventionalEmulsionCheckBox = QtWidgets.QCheckBox("Seeded free-radical emulsion")
        self.conventionalEmulsionCheckBox.setAutoExclusive(True)
        self.conventionalEmulsionCheckBox.setChecked(True)
        self.RAFTEmulsionCheckBox = QtWidgets.QCheckBox("RAFT emulsion")
        self.RAFTEmulsionCheckBox.setAutoExclusive(True)
        self.RAFTDispersionCheckBox = QtWidgets.QCheckBox("RAFT dispersion")
        self.RAFTDispersionCheckBox.setAutoExclusive(True)

        # Add checkBoxes to layout
        self.polymerisationSelectorLayout.addWidget(self.conventionalEmulsionCheckBox, 0, 0)
        self.polymerisationSelectorLayout.addWidget(self.RAFTEmulsionCheckBox, 1, 0)
        self.polymerisationSelectorLayout.addWidget(self.RAFTDispersionCheckBox, 2, 0)

        # Connect checkBox clicks to widget re-formatting function
        self.conventionalEmulsionCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))
        self.RAFTEmulsionCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))

        ##################################### Method selector box ##############################################

        # Create groupBox for selecting the type of experiment to run and add grid layout
        methodSelectorBox = QtWidgets.QGroupBox("Select method type")
        methodSelectorBox.setMaximumSize(150, 120)
        methodBoxLayout = QtWidgets.QGridLayout()
        methodSelectorBox.setLayout(methodBoxLayout)

        # Create checkBoxes for each type of experiment and setAutoExclusive so only one can be selected at a time
        self.OFAATCheckBox = QtWidgets.QCheckBox("OFAAT")
        self.OFAATCheckBox.setAutoExclusive(True)
        self.DoECheckBox = QtWidgets.QCheckBox("DoE")
        self.DoECheckBox.setAutoExclusive(True)
        self.TSEMOCheckBox = QtWidgets.QCheckBox("TSEMO")
        self.TSEMOCheckBox.setAutoExclusive(True)
        self.kineticStudyCheckBox = QtWidgets.QCheckBox("Kinetics study")
        self.kineticStudyCheckBox.setAutoExclusive(True)

        # Add checkBoxes to layout
        methodBoxLayout.addWidget(self.OFAATCheckBox, 0, 0)
        methodBoxLayout.addWidget(self.DoECheckBox, 1, 0)
        methodBoxLayout.addWidget(self.TSEMOCheckBox, 2, 0)
        methodBoxLayout.addWidget(self.kineticStudyCheckBox, 3, 0)

        # Connect checkBox clicks to widget re-formatting function
        for checkBox in methodSelectorBox.findChildren(QtWidgets.QCheckBox):
            checkBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))
            checkBox.clicked.connect(lambda: self.formatWidget(groupBox=self.kineticsBox))
            checkBox.clicked.connect(lambda: self.formatWidget(groupBox=self.propertyTargetingBox))

        ###################################### Variable selector box ############################################

        # Create groupBox to contain all variable fields
        self.variableSelectorBox = QtWidgets.QGroupBox("Define variables")
        self.variableSelectorBox.resize(300, 400)
        self.variableSelectorLayout = QtWidgets.QGridLayout()
        self.variableSelectorBox.setLayout(self.variableSelectorLayout)

        ########### Common variables #############
        self.reactor1VolumeLabel = QtWidgets.QLabel("CSTR volume (mL)")
        self.reactor1VolumeText = QtWidgets.QLineEdit("1.5")

        self.nCSTRsLabel = QtWidgets.QLabel("Number of CSTRs")
        self.nCSTRsText = QtWidgets.QLineEdit("5")

        self.reactor2VolumeLabel = QtWidgets.QLabel("Reactor 2 volume (mL)")
        self.reactor2VolumeText = QtWidgets.QLineEdit("5")

        self.sampleVolumeLabel = QtWidgets.QLabel("Sample volume (mL)")
        self.sampleVolumeText = QtWidgets.QLineEdit("5")
        
        self.initiatorStockConcLabel = QtWidgets.QLabel("Initiator stock concentration (w/w)")
        self.initiatorStockConcText = QtWidgets.QLineEdit("0.025")
        self.initiatorStockConcText.setToolTip("Concentration (w/w) of initiator only in the initiator stock solution")
        self.w_initiatorStockLabel = QtWidgets.QLabel("Initiator stock total solids (w/w)")
        self.w_initiatorStockText = QtWidgets.QLineEdit("0.025")
        self.w_initiatorStockText.setToolTip("Concentration (w/w) of all components in the initiator stock solution")
        self.initiatorRatioCheckBox = QtWidgets.QCheckBox("Initiator ratio (pphm)")
        self.initiatorRatioCheckBox.setToolTip("Target concentration (w/w) of initiator in the prepared emulsion")
        self.initiatorRatioText = QtWidgets.QLineEdit("0.65")
        self.initiatorRatioText1 = QtWidgets.QLineEdit("Min")
        self.initiatorRatioText2 = QtWidgets.QLineEdit("Max")

        self.monomerRatioText = QtWidgets.QLineEdit("1.0")
        self.monomerRatioText1 = QtWidgets.QLineEdit("Min A")
        self.monomerRatioText2 = QtWidgets.QLineEdit("Max A")
        self.monomerDensityLabel = QtWidgets.QLabel("Monomer density [g/mL]")
        self.monomerDensityText = QtWidgets.QLineEdit("1")
        self.monomerDensityLabel1 = QtWidgets.QLabel("Monomer A density [g/mL]")
        self.monomerDensityText1 = QtWidgets.QLineEdit("Density A")
        self.monomerDensityLabel2 = QtWidgets.QLabel("Monomer B density [g/mL]")
        self.monomerDensityText2 = QtWidgets.QLineEdit("Density B")

        self.residenceTimeCheckBox = QtWidgets.QCheckBox("Residence time (min)")
        self.residenceTimeText = QtWidgets.QLineEdit("100")
        self.residenceTimeText1 = QtWidgets.QLineEdit("Min")
        self.residenceTimeText2 = QtWidgets.QLineEdit("Max")

        self.productConcCheckBox = QtWidgets.QCheckBox("Product concentration [g/g]")
        self.productConcentrationText = QtWidgets.QLineEdit("0.3")
        self.productConcentrationText1 = QtWidgets.QLineEdit("Min")
        self.productConcentrationText2 = QtWidgets.QLineEdit("Max")

        self.w_makeupWaterLabel = QtWidgets.QLabel("Makeup water total solids [w/w]")
        self.w_makeupWaterText = QtWidgets.QLineEdit("0.001")
        self.w_makeupWaterText.setToolTip("Concentration [w/w] of additional components added to the makeup water")

        self.feedRateCheckBox = QtWidgets.QCheckBox("Number of feeds")
        self.numFeedsText = QtWidgets.QLineEdit("5")
        self.numFeedsText1 = QtWidgets.QLineEdit("Min")
        self.numFeedsText2 = QtWidgets.QLineEdit("Max")

        self.numExpsLabel = QtWidgets.QLabel("Number of experiments")
        self.numExpsText = QtWidgets.QLineEdit("1")        

        ########### Conventional emulsion polymerisation variables ##############
        # Core variables #
        self.seedConcentrationLabel = QtWidgets.QLabel("Seed concentration [w/w]")
        self.seedConcentrationText = QtWidgets.QLineEdit("0.0935")

        self.w_eMaxLabel = QtWidgets.QLabel("Max emulsion concentration [w/w]")
        self.w_eMaxText = QtWidgets.QLineEdit("0.7")

        # Optional variables #
        self.surfactantStockConcLabel = QtWidgets.QLabel("Surfactant stock concentration [w/w]")
        self.surfactantStockConcText = QtWidgets.QLineEdit("0.05")
        self.surfactantStockConcText.setToolTip("Concentration [w/w] of surfactant only in the surfactant stock solution")
        self.w_surfactantStockLabel = QtWidgets.QLabel("Surfactant stock total solids [w/w]")
        self.w_surfactantStockText = QtWidgets.QLineEdit("0.05")
        self.w_surfactantStockText.setToolTip("Concentration (w/w) of all components in the surfactant stock solution")
        self.surfactantRatioCheckBox = QtWidgets.QCheckBox("Surfactant ratio (pphm)")
        self.surfactantRatioCheckBox = QtWidgets.QCheckBox("Surfactant ratio [pphm]")
        self.surfactantRatioText = QtWidgets.QLineEdit("1.7")
        self.surfactantRatioText1 = QtWidgets.QLineEdit("Min")
        self.surfactantRatioText2 = QtWidgets.QLineEdit("Max")

        self.seedAmountText = QtWidgets.QLineEdit("0.05")
        self.seedAmountText1 = QtWidgets.QLineEdit("Min R value")
        self.seedAmountText2 = QtWidgets.QLineEdit("Max R value")

        self.monomerRatioCheckBox = QtWidgets.QCheckBox("Monomer ratio [g/g]")
        self.seedAmountCheckBox = QtWidgets.QCheckBox("Seed fraction [g/g]")

        ########### RAFT emulsion polymerisation variables ##############
        # Core variables #
        self.seededQueryCheckBox = QtWidgets.QCheckBox("Using seed?")
        self.macroCTAStockConcLabel = QtWidgets.QLabel("MacroCTA stock concentration (w/w)")
        self.macroCTAStockConcText = QtWidgets.QLineEdit("0.1")
        self.macroCTAStockConcText.setToolTip("Concentration (w/w) of macroCTA only in the macroCTA stock solution")
        self.macroCTAPolymerConcLabel = QtWidgets.QLabel("MacroCTA polymer concentration (w/w)")
        self.macroCTAPolymerConcText = QtWidgets.QLineEdit("0.117")
        self.macroCTAMonomerConcLabel = QtWidgets.QLabel("MacroCTA monomer concentration (w/w)")
        self.macroCTAMonomerConcText = QtWidgets.QLineEdit("0.011022")
        self.macroCTAPolymerConcLabel.setHidden(True)
        self.macroCTAPolymerConcText.setHidden(True)
        self.macroCTAMonomerConcLabel.setHidden(True)
        self.macroCTAMonomerConcText.setHidden(True)
        self.molarMassMacroCTALabel = QtWidgets.QLabel("MacroCTA molar mass (g/mol)")
        self.molarMassMacroCTAText = QtWidgets.QLineEdit("23123.75")

        self.molarMassMonomerLabel = QtWidgets.QLabel("Monomer molar mass (g/mol)")
        self.molarMassMonomerText = QtWidgets.QLineEdit("128.171")

        self.molarMassInitiatorLabel = QtWidgets.QLabel("Initiator molar mass (g/mol)")
        self.molarMassInitiatorText = QtWidgets.QLineEdit("280.28")

        # Optional variables
        self.DPCheckBox = QtWidgets.QCheckBox("DP")
        self.DPText = QtWidgets.QLineEdit("300")
        self.DPText1 = QtWidgets.QLineEdit("Min")
        self.DPText2 = QtWidgets.QLineEdit("Max")

        self.ctaToInitiatorCheckBox = QtWidgets.QCheckBox("CTA:I")
        self.ctaToInitiatorText = QtWidgets.QLineEdit("5")
        self.ctaToInitiatorText1 = QtWidgets.QLineEdit("Min")
        self.ctaToInitiatorText2 = QtWidgets.QLineEdit("Max")
        
        self.calculateParamsBtn = QtWidgets.QPushButton("\nCalculate parameters\n")
        self.calculateParamsBtn.setStyleSheet(
            "background-color: #e9e9e9;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #6e6e6e;"
            "font-size: 10pt"
        )

        self.variableSelectorLayout.addWidget(self.reactor1VolumeLabel, 0, 0)
        self.variableSelectorLayout.addWidget(self.reactor1VolumeText, 0, 1)
        self.variableSelectorLayout.addWidget(self.nCSTRsLabel, 1, 0)
        self.variableSelectorLayout.addWidget(self.nCSTRsText, 1, 1)
        self.variableSelectorLayout.addWidget(self.reactor2VolumeLabel, 2, 0)
        self.variableSelectorLayout.addWidget(self.reactor2VolumeText, 2, 1)
        self.variableSelectorLayout.addWidget(self.sampleVolumeLabel, 3, 0)
        self.variableSelectorLayout.addWidget(self.sampleVolumeText, 3, 1)
        self.variableSelectorLayout.addWidget(self.w_eMaxLabel, 4, 0)
        self.variableSelectorLayout.addWidget(self.w_eMaxText, 4, 1)
        self.variableSelectorLayout.addWidget(self.monomerDensityLabel, 5, 0)
        self.variableSelectorLayout.addWidget(self.monomerDensityText, 5, 1)
        self.variableSelectorLayout.addWidget(self.residenceTimeCheckBox, 6, 0)
        self.variableSelectorLayout.addWidget(self.residenceTimeText, 6, 1)
        self.variableSelectorLayout.addWidget(self.residenceTimeText1, 6, 1)
        self.variableSelectorLayout.addWidget(self.residenceTimeText2, 6, 2)
        self.variableSelectorLayout.addWidget(self.productConcCheckBox, 7, 0)
        self.variableSelectorLayout.addWidget(self.productConcentrationText, 7, 1)
        self.variableSelectorLayout.addWidget(self.productConcentrationText1, 7, 1)
        self.variableSelectorLayout.addWidget(self.productConcentrationText2, 7, 2)
        self.variableSelectorLayout.addWidget(self.initiatorStockConcLabel, 8, 0)
        self.variableSelectorLayout.addWidget(self.initiatorStockConcText, 8, 1)
        self.variableSelectorLayout.addWidget(self.w_initiatorStockLabel, 9, 0)
        self.variableSelectorLayout.addWidget(self.w_initiatorStockText, 9, 1)
        self.variableSelectorLayout.addWidget(self.initiatorRatioCheckBox, 10, 0)
        self.variableSelectorLayout.addWidget(self.initiatorRatioText, 10, 1)
        self.variableSelectorLayout.addWidget(self.initiatorRatioText1, 10, 1)
        self.variableSelectorLayout.addWidget(self.initiatorRatioText2, 10, 2)
        self.variableSelectorLayout.addWidget(self.w_makeupWaterLabel, 11, 0)
        self.variableSelectorLayout.addWidget(self.w_makeupWaterText, 11, 1)
        self.variableSelectorLayout.addWidget(self.seedConcentrationLabel, 12, 0)
        self.variableSelectorLayout.addWidget(self.seedConcentrationText, 12, 1)
        self.variableSelectorLayout.addWidget(self.surfactantStockConcLabel, 13, 0)
        self.variableSelectorLayout.addWidget(self.surfactantStockConcText, 13, 1)
        self.variableSelectorLayout.addWidget(self.w_surfactantStockLabel, 14, 0)
        self.variableSelectorLayout.addWidget(self.w_surfactantStockText, 14, 1)
        self.variableSelectorLayout.addWidget(self.surfactantRatioCheckBox, 15, 0)
        self.variableSelectorLayout.addWidget(self.surfactantRatioText, 15, 1)
        self.variableSelectorLayout.addWidget(self.surfactantRatioText1, 15, 1)
        self.variableSelectorLayout.addWidget(self.surfactantRatioText2, 15, 2)
        self.variableSelectorLayout.addWidget(self.monomerRatioCheckBox, 16, 0)
        self.variableSelectorLayout.addWidget(self.monomerRatioText, 16, 1)
        self.variableSelectorLayout.addWidget(self.monomerRatioText1, 16, 1)
        self.variableSelectorLayout.addWidget(self.monomerDensityText1, 17, 1)
        self.variableSelectorLayout.addWidget(self.monomerRatioText2, 16, 2)
        self.variableSelectorLayout.addWidget(self.monomerDensityText2, 17, 2)
        self.variableSelectorLayout.addWidget(self.seedAmountCheckBox, 18, 0)
        self.variableSelectorLayout.addWidget(self.seedAmountText, 18, 1)
        self.variableSelectorLayout.addWidget(self.seedAmountText1, 18, 1)
        self.variableSelectorLayout.addWidget(self.seedAmountText2, 18, 2)
        self.variableSelectorLayout.addWidget(self.feedRateCheckBox, 19, 0)
        self.variableSelectorLayout.addWidget(self.numFeedsText, 19, 1)
        self.variableSelectorLayout.addWidget(self.numFeedsText1, 19, 1)
        self.variableSelectorLayout.addWidget(self.numFeedsText2, 19, 2)
        self.variableSelectorLayout.addWidget(self.seededQueryCheckBox, 20, 0)
        self.variableSelectorLayout.addWidget(self.macroCTAStockConcLabel, 21, 0)
        self.variableSelectorLayout.addWidget(self.macroCTAStockConcText, 21, 1)
        self.variableSelectorLayout.addWidget(self.macroCTAPolymerConcLabel, 22, 0)
        self.variableSelectorLayout.addWidget(self.macroCTAPolymerConcText, 22, 1)
        self.variableSelectorLayout.addWidget(self.macroCTAMonomerConcLabel, 23, 0)
        self.variableSelectorLayout.addWidget(self.macroCTAMonomerConcText, 23, 1)
        self.variableSelectorLayout.addWidget(self.molarMassMacroCTALabel, 24, 0)
        self.variableSelectorLayout.addWidget(self.molarMassMacroCTAText, 24, 1)
        self.variableSelectorLayout.addWidget(self.molarMassInitiatorLabel, 25, 0)
        self.variableSelectorLayout.addWidget(self.molarMassInitiatorText, 25, 1)
        self.variableSelectorLayout.addWidget(self.molarMassMonomerLabel, 26, 0)
        self.variableSelectorLayout.addWidget(self.molarMassMonomerText, 26, 1)
        self.variableSelectorLayout.addWidget(self.DPCheckBox, 27, 0)
        self.variableSelectorLayout.addWidget(self.DPText, 27, 1)
        self.variableSelectorLayout.addWidget(self.DPText1, 27, 1)
        self.variableSelectorLayout.addWidget(self.DPText2, 27, 2)
        self.variableSelectorLayout.addWidget(self.ctaToInitiatorCheckBox, 28, 0)
        self.variableSelectorLayout.addWidget(self.ctaToInitiatorText, 28, 1)
        self.variableSelectorLayout.addWidget(self.ctaToInitiatorText1, 28, 1)
        self.variableSelectorLayout.addWidget(self.ctaToInitiatorText2, 28, 2)
        self.variableSelectorLayout.addWidget(self.numExpsLabel, 29, 0)
        self.variableSelectorLayout.addWidget(self.numExpsText, 29, 1)
        self.variableSelectorLayout.addWidget(self.calculateParamsBtn, 30, 0, 1, 2)

        self.residenceTimeCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))
        self.productConcCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))
        self.initiatorRatioCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))

        self.surfactantRatioCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))
        self.monomerRatioCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))
        self.seedAmountCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))
        self.feedRateCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))

        self.seededQueryCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))
        self.DPCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))
        self.ctaToInitiatorCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.variableSelectorBox))

        self.productConcentrationText.textEdited.connect(self.maxSeedFracCalculator) # Actively sets maximum possible seed fraction based on the fixed product conc selected
        self.productConcentrationText2.textEdited.connect(self.maxSeedFracCalculator) # Actively sets the maximum possible seed fraction based on max product conc selected
        self.w_eMaxText.textEdited.connect(self.maxSeedFracCalculator) # Actively sets the maximum possible seed fraction based on maximum emulsion concentration input
        self.seedConcentrationText.textEdited.connect(self.maxSeedFracCalculator) # Actively sets the maximum possible seed fraction based on seed concentration input
        self.seedAmountText2.textEdited.connect(self.maxSeedFracValidator) # Actively validates the upper bound of seed fraction according to the maximum specified by w_e, w_s, and w_f

        self.calculateParamsBtn.clicked.connect(self.buildExperiment)
        
        # Reset widget by default until method selections are made
        for label in self.variableSelectorBox.findChildren(QtWidgets.QLabel):
            label.setHidden(True)
        for lineEdit in self.variableSelectorBox.findChildren(QtWidgets.QLineEdit):
            lineEdit.setHidden(True)
        for checkBox in self.variableSelectorBox.findChildren(QtWidgets.QCheckBox):
            checkBox.setHidden(True)

        ################################# Property targeting box #########################################
            
        # Create groupBox to contain all 'optimisable' properties
        self.propertyTargetingBox = QtWidgets.QGroupBox("Select optimisation targets")
        self.propertyTargetingBox.resize(300, 400)
        self.propertyTargetingLayout = QtWidgets.QGridLayout()
        self.propertyTargetingBox.setLayout(self.propertyTargetingLayout)

        self.particleSizeCheckBox = QtWidgets.QCheckBox("Particle size [nm]")
        self.particleSizeText = QtWidgets.QLineEdit("Size target")
        self.particleSizeMappingCheckBox = QtWidgets.QCheckBox("Enable size mapping (Flip-Flop)")

        # Need conversion monitoring to optimise for conversion, and GPC for dispersity optimisation...

        self.surfactantConcentrationObjectiveCheckBox = QtWidgets.QCheckBox("Minimise surfactant [pphm]")
        self.seedFractionObjectiveCheckBox = QtWidgets.QCheckBox("Minimise seed [pphm]")
        self.initiatorRatioObjectiveCheckBox = QtWidgets.QCheckBox("Minimise initiator ratio [pphm]")
        self.ctaToInitiatorObjectiveCheckBox = QtWidgets.QCheckBox("Maximise CTA : initiator [mol/mol]")

        self.TgCheckBox = QtWidgets.QCheckBox("Tg (\u00b0C)")
        self.TgText = QtWidgets.QLineEdit("Tg target")
        self.TgPolymer1Text = QtWidgets.QLineEdit("Tg polymer 1")
        self.TgPolymer2Text = QtWidgets.QLineEdit("Tg polymer 2")
        self.fractionFixedCopolymer1Label = QtWidgets.QLabel("Fixed copolymer 1 conc [w/w]")
        self.fractionFixedCopolymer1Text = QtWidgets.QLineEdit("0.01")
        self.fractionFixedCopolymer2Label = QtWidgets.QLabel("Fixed copolymer 2 conc [w/w]")
        self.fractionFixedCopolymer2Text = QtWidgets.QLineEdit("0.02")
        self.TgFixedCopolymer1Label = QtWidgets.QLabel("Fixed copolymer 1 Tg [\u00b0C]")
        self.TgFixedCopolymer1Text = QtWidgets.QLineEdit("106")
        self.TgFixedCopolymer2Label = QtWidgets.QLabel("Fixed copolymer 2 Tg [\u00b0C]")
        self.TgFixedCopolymer2Text = QtWidgets.QLineEdit("228")
        self.calculateMonomerCompositionBtn = QtWidgets.QPushButton("Calculate monomer composition")
        self.monomerCompositionLabel1 = QtWidgets.QLabel("Monomer composition: ")
        self.monomerCompositionLabel2 = QtWidgets.QLabel()
        self.calculateMonomerCompositionBtn.setStyleSheet("background-color: #e9e9e9;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #6e6e6e;"
            "font-size: 10pt")
        self.calculateMonomerCompositionBtn.setFixedSize(200, 50)

        self.trainingDataBrowseButton = QtWidgets.QPushButton('Select training data folder')
        self.trainingDataBrowseButton.setStyleSheet(
            "background-color: #e9e9e9;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #6e6e6e"
        )
        self.trainingDataBrowseButton.clicked.connect(self.searchTrainingData)
        self.trainingDataBrowseButton.setFixedSize(200, 50)
        self.trainingDataPathText = QtWidgets.QLineEdit()

        self.maxIterationsLabel = QtWidgets.QLabel("Max iterations")
        self.maxIterationsText = QtWidgets.QLineEdit("30")

        self.sampleStrategyLabel = QtWidgets.QLabel("Sample every __ experiments")
        self.sampleStrategyText = QtWidgets.QLineEdit("3")

        self.propertyTargetingLayout.addWidget(self.particleSizeCheckBox, 0, 0)
        self.propertyTargetingLayout.addWidget(self.particleSizeText, 0, 1)
        self.propertyTargetingLayout.addWidget(self.particleSizeMappingCheckBox, 1, 0)
        self.propertyTargetingLayout.addWidget(self.surfactantConcentrationObjectiveCheckBox, 2, 0)
        self.propertyTargetingLayout.addWidget(self.seedFractionObjectiveCheckBox, 3, 0)
        self.propertyTargetingLayout.addWidget(self.initiatorRatioObjectiveCheckBox, 4, 0)
        self.propertyTargetingLayout.addWidget(self.ctaToInitiatorObjectiveCheckBox, 5, 0)
        self.propertyTargetingLayout.addWidget(self.TgCheckBox, 6, 0)
        self.propertyTargetingLayout.addWidget(self.TgText, 6, 1)
        self.propertyTargetingLayout.addWidget(self.TgPolymer1Text, 7, 1)
        self.propertyTargetingLayout.addWidget(self.TgPolymer2Text, 8, 1)
        self.propertyTargetingLayout.addWidget(self.fractionFixedCopolymer1Label, 9, 0)
        self.propertyTargetingLayout.addWidget(self.fractionFixedCopolymer1Text, 9, 1)
        self.propertyTargetingLayout.addWidget(self.fractionFixedCopolymer2Label, 10, 0)
        self.propertyTargetingLayout.addWidget(self.fractionFixedCopolymer2Text, 10, 1)
        self.propertyTargetingLayout.addWidget(self.TgFixedCopolymer1Label, 11, 0)
        self.propertyTargetingLayout.addWidget(self.TgFixedCopolymer1Text, 11, 1)
        self.propertyTargetingLayout.addWidget(self.TgFixedCopolymer2Label, 12, 0)
        self.propertyTargetingLayout.addWidget(self.TgFixedCopolymer2Text, 12, 1)
        self.propertyTargetingLayout.addWidget(self.calculateMonomerCompositionBtn, 13, 0, 1, 2)
        self.propertyTargetingLayout.addWidget(self.monomerCompositionLabel1, 14, 0)
        self.propertyTargetingLayout.addWidget(self.monomerCompositionLabel2, 15, 0)
        self.propertyTargetingLayout.addWidget(self.trainingDataBrowseButton, 16, 0, 1, 2)
        self.propertyTargetingLayout.addWidget(self.trainingDataPathText, 17, 0, 1, 2)
        self.propertyTargetingLayout.addWidget(self.maxIterationsLabel, 18, 0)
        self.propertyTargetingLayout.addWidget(self.maxIterationsText, 18, 1)
        self.propertyTargetingLayout.addWidget(self.sampleStrategyLabel, 19, 0)
        self.propertyTargetingLayout.addWidget(self.sampleStrategyText, 19, 1)

        self.propertyTargetingBox.setHidden(True)

        for label in self.propertyTargetingBox.findChildren(QtWidgets.QLabel):
            label.setHidden(True)
        for lineEdit in self.propertyTargetingBox.findChildren(QtWidgets.QLineEdit):
            lineEdit.setHidden(True)
        for checkBox in self.propertyTargetingBox.findChildren(QtWidgets.QCheckBox):
            checkBox.setHidden(True)
        self.calculateMonomerCompositionBtn.setHidden(True)
        self.trainingDataBrowseButton.setHidden(True)

        self.particleSizeCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.propertyTargetingBox))
        self.particleSizeMappingCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.propertyTargetingBox))
        self.TgCheckBox.clicked.connect(lambda: self.formatWidget(groupBox=self.propertyTargetingBox))
        self.calculateMonomerCompositionBtn.clicked.connect(self.calculateMonomerComposition)

        ###################################### Kinetic screening box #####################################

        self.kineticsBox = QtWidgets.QGroupBox("Kinetic screening")
        self.kineticsBox.resize(300, 400)
        self.kineticsLayout = QtWidgets.QGridLayout()
        self.kineticsBox.setLayout(self.kineticsLayout)

        self.numberKineticsSamplesLabel = QtWidgets.QLabel("Number of kinetics samples")
        self.numberKineticsSamplesLineEdit = QtWidgets.QLineEdit("6")

        self.kineticsLayout.addWidget(self.numberKineticsSamplesLabel, 0, 0)
        self.kineticsLayout.addWidget(self.numberKineticsSamplesLineEdit, 0, 1)

        self.kineticsBox.setHidden(True)

        ###################################### Cleaning box ##############################################
        cleaningBox = QtWidgets.QGroupBox("Cleaning")
        cleaningBox.resize(300, 150)
        cleaningBox.setMaximumWidth(300) 
        self.cleaningBoxLayout = QtWidgets.QGridLayout()
        cleaningBox.setLayout(self.cleaningBoxLayout)

        self.cleaningRateLabel = QtWidgets.QLabel("Cleaning flow rate [mL/min]")
        self.cleaningRateText = QtWidgets.QLineEdit("5")
        self.cleaningRateText.setMaximumWidth(50)

        self.cleaningVolumesLabel = QtWidgets.QLabel("Number of reactor volumes")
        self.cleaningVolumesText = QtWidgets.QLineEdit("2")
        self.cleaningVolumesText.setMaximumWidth(50)

        self.samplesNumberLabel = QtWidgets.QLabel("Number of samples")
        self.samplesNumberText = QtWidgets.QLineEdit("7")
        self.samplesNumberText.setMaximumWidth(50)

        self.cleanSampleTubesButton = QtWidgets.QPushButton("\nClean autosampler\n")
        self.cleanSampleTubesButton.setStyleSheet(
            "background-color: #e9e9e9;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #6e6e6e;"
            "font-size: 10pt"
        )

        self.emptyWasteButton = QtWidgets.QPushButton("\nWaste bottle emptied\n")
        self.emptyWasteButton.setStyleSheet(
            "background-color: #e9e9e9;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #6e6e6e;"
            "font-size: 10pt"
        )

        self.wasteVolumeCounterLabel = QtWidgets.QLabel("Waste container volume\navailable (mL):")
        self.wasteVolumeCounterLabel.setStyleSheet(
            "font-size: 10pt"
        )
        self.wasteVolumeCounter = QtWidgets.QLabel()
        self.wasteVolumeCounter.setStyleSheet(
            "font-size: 10pt"
        )

        self.seedFillButton = QtWidgets.QPushButton("\nSeed syringe refilled\n")
        self.seedFillButton.setStyleSheet(
            "background-color: #e9e9e9;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #6e6e6e;"
            "font-size: 10pt"
        )

        self.seedVolumeCounterLabel = QtWidgets.QLabel("Seed syringe volume\navailable (mL):")
        self.seedVolumeCounterLabel.setStyleSheet(
            "font-size: 10pt"
        )
        self.seedVolumeCounter = QtWidgets.QLabel()
        self.seedVolumeCounter.setStyleSheet(
            "font-size: 10pt"
        )

        self.cleaningBoxLayout.addWidget(self.cleaningRateLabel, 0, 0)
        self.cleaningBoxLayout.addWidget(self.cleaningRateText, 0, 1)
        self.cleaningBoxLayout.addWidget(self.cleaningVolumesLabel, 1, 0)
        self.cleaningBoxLayout.addWidget(self.cleaningVolumesText, 1, 1)
        self.cleaningBoxLayout.addWidget(self.samplesNumberLabel, 2, 0)
        self.cleaningBoxLayout.addWidget(self.samplesNumberText, 2, 1)
        self.cleaningBoxLayout.addWidget(self.cleanSampleTubesButton, 3, 0, 1, 2)
        self.cleaningBoxLayout.addWidget(self.emptyWasteButton, 4, 0, 1, 2)
        self.cleaningBoxLayout.addWidget(self.wasteVolumeCounterLabel, 5, 0, QtCore.Qt.AlignTop)
        self.cleaningBoxLayout.addWidget(self.wasteVolumeCounter, 6, 0, 1, 2, QtCore.Qt.AlignTop)
        self.cleaningBoxLayout.addWidget(self.seedFillButton, 7, 0, 1, 2)
        self.cleaningBoxLayout.addWidget(self.seedVolumeCounterLabel, 8, 0, QtCore.Qt.AlignTop)
        self.cleaningBoxLayout.addWidget(self.seedVolumeCounter, 9, 0, 1, 2, QtCore.Qt.AlignTop)

        self.Layout.addWidget(cleaningBox, 0, 2, 1, 1)

        self.cleanSampleTubesButton.clicked.connect(self.startSamplerCleaning)
        self.emptyWasteButton.clicked.connect(self.wasteEmptied)
        self.seedFillButton.clicked.connect(self.seedFilled)

        ###################################### Calculated variables ############################################
        calcVarBox = QtWidgets.QGroupBox("Calculated variables")
        calcVarBox.setMaximumSize(680, 250)
        calcVarBoxLayout = QtWidgets.QGridLayout()
        calcVarBox.setLayout(calcVarBoxLayout)

        self.parametersTable = QtWidgets.QTableWidget()
        self.parametersTable.setMinimumWidth(600)

        # calcVarBoxLayout.addWidget(self.parametersTable, 0, 0)

        ####################################### Save directory box #############################################
        self.directoryLineEdit = QtWidgets.QLineEdit()
        self.directoryLineEdit.setMaximumSize(300, 50)
        self.directoryBrowseButton = QtWidgets.QPushButton('Browse file locations')
        self.directoryBrowseButton.setMaximumSize(120, 50)
        self.directoryBrowseButton.setStyleSheet(
            "background-color: #e9e9e9;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #6e6e6e"
        )
        self.directoryBrowseButton.clicked.connect(self.searchDirectories)
        self.experimentIDLineEdit = QtWidgets.QLineEdit()
        self.experimentIDLineEdit.setMaximumSize(120, 50)
        self.experimentIDLineEdit.setText(str(datetime.date.today()))
        self.setDataPathButton = QtWidgets.QPushButton('Set data path')
        self.setDataPathButton.setMaximumSize(100, 50)
        self.setDataPathButton.setStyleSheet(
            "background-color: #e9e9e9;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #6e6e6e"
        )
        self.setDataPathButton.clicked.connect(self.setDataPath)
        self.directoryBox = QtWidgets.QGroupBox('Set save location')
        self.directoryBox.setMaximumSize(680, 50)
        directoryBoxLayout = QtWidgets.QHBoxLayout()
        self.directoryBox.setLayout(directoryBoxLayout)
        directoryBoxLayout.addWidget(self.directoryLineEdit, QtCore.Qt.AlignTop)
        directoryBoxLayout.addWidget(self.directoryBrowseButton, QtCore.Qt.AlignTop)
        directoryBoxLayout.addWidget(self.experimentIDLineEdit, QtCore.Qt.AlignTop)
        directoryBoxLayout.addWidget(self.setDataPathButton, QtCore.Qt.AlignTop)

        self.savePathCreated = 0 # Used to track if a save file has been created to store metadata and results

        ####################################### Dialogue widget ################################################
        # Used for displaying messages to the user

        self.dialogueLabel = QtWidgets.QLabel()
        self.ackButton = QtWidgets.QPushButton("OK")
        self.ackButton.setFixedSize(60, 30)
        self.ackButton.setStyleSheet(
            "background-color: #e9e9e9;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #6e6e6e;"
            "font-size: 10pt"
        )
        self.dialogueLabel.setHidden(True)
        self.dialogueLabel.setStyleSheet(
            "background-color: #fcfc03;"
            "font-size: 10pt;"
            "border-radius: 5px"
        )
        self.ackButton.setHidden(True)        

        ####################################### Add methodBuilder widgets #############################################
        self.runButton = QtWidgets.QPushButton("\nStart experiment\n")
        self.runButton.setStyleSheet(
            "background-color: #0996D4;" 
            "color: white;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #2d5e86;"
            "font-size: 10pt")
        self.runButton.setFixedSize(200, 50)
        self.stopButton = QtWidgets.QPushButton("\nStop experiment\n")
        self.stopButton.setStyleSheet(
            "background-color: #d27b79;"
            "color: white;"
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #862f2d;"
            "font-size: 10pt")
        self.stopButton.setFixedSize(200, 50)
        methodBuilderLayout.addWidget(self.directoryBox, 0, 0, 1, 6, QtCore.Qt.AlignTop)
        methodBuilderLayout.addWidget(methodSelectorBox, 1, 0, 1, 2, QtCore.Qt.AlignTop)
        methodBuilderLayout.addWidget(self.polymerisationSelectorBox, 1, 2, 1, 2, QtCore.Qt.AlignTop)
        methodBuilderLayout.addWidget(self.variableSelectorBox, 2, 0, 2, 4, QtCore.Qt.AlignTop)
        methodBuilderLayout.addWidget(self.propertyTargetingBox, 2, 4, QtCore.Qt.AlignTop)
        methodBuilderLayout.addWidget(self.kineticsBox, 2, 4, QtCore.Qt.AlignTop)
        methodBuilderLayout.addWidget(self.runButton, 5, 1, 1, 2, QtCore.Qt.AlignTop)
        methodBuilderLayout.addWidget(self.stopButton, 5, 3, 1, 2, QtCore.Qt.AlignTop)
        methodBuilderLayout.addWidget(self.dialogueLabel, 6, 1, 1, 3, QtCore.Qt.AlignTop)
        methodBuilderLayout.addWidget(self.ackButton, 6, 4, 1, 2, QtCore.Qt.AlignTop)
        self.Layout.addWidget(self.methodBuilderBox, 0, 0, 2, 1, QtCore.Qt.AlignTop)

        self.runButton.clicked.connect(self.startExperiment)
        self.stopButton.clicked.connect(self.stopExperimentManual)
        self.ackButton.clicked.connect(self.unpauseExperiment)

        ######################################## Method graph widget ############################################
        self.methodGraphBox = QtWidgets.QGroupBox("Experimental design")
        self.methodGraphLayout = QtWidgets.QGridLayout()
        self.methodGraphBox.setLayout(self.methodGraphLayout)
        self.Layout.addWidget(self.methodGraphBox, 0, 1, 1, 1)
        self.methodGraphLayout.addWidget(self.parametersTable, 0, 0, 2, 2, QtCore.Qt.AlignLeft)
        self.nGraph = 0 # To track number of times experimental graph has been constructed - if first time running the app, n = 0. Prevents multiple graphs being displayed

        ######################################## Platform diagram widget #########################################

        self.diaGroupBox = QtWidgets.QGroupBox("Experimental set-up")
        self.diaGroupBoxLayout = QtWidgets.QGridLayout()
        self.diaGroupBox.setLayout(self.diaGroupBoxLayout)
        diagramLoc = r'C:\Users\pm15pmp\Miniconda3\envs\MPSR\polymer_platform\PCubed\Emulsion_FRP.png'
        self.diagram = QtGui.QPixmap(diagramLoc)
        self.diaLabel = QtWidgets.QLabel(self)
        self.diaLabel.setPixmap(self.diagram)
        self.diaLabel.setScaledContents(True)
        self.diaLabel.setFixedSize(1050, 430)

        self.pumpAq1FlowRateText = QtWidgets.QLineEdit()
        self.pumpAq1FlowRateText.setFixedSize(35, 20)
        self.pumpAq2FlowRateText = QtWidgets.QLineEdit()
        self.pumpAq2FlowRateText.setFixedSize(35, 20)
        self.pumpM1FlowRateText = QtWidgets.QLineEdit()
        self.pumpM1FlowRateText.setFixedSize(35, 20)
        self.pumpM2FlowRateText = QtWidgets.QLineEdit()
        self.pumpM2FlowRateText.setFixedSize(35, 20)
        self.pumpSolventFlowRateText = QtWidgets.QLineEdit()
        self.pumpSolventFlowRateText.setFixedSize(35, 20)
        self.pumpSeedFlowRateText = QtWidgets.QLineEdit()
        self.pumpSeedFlowRateText.setFixedSize(35, 20)
        self.Aq1StartBtn = QtWidgets.QPushButton("")
        self.Aq1StartBtn.setFixedSize(16, 16)
        self.Aq1StartBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #13ad2a")
        self.Aq1StopBtn = QtWidgets.QPushButton("")
        self.Aq1StopBtn.setFixedSize(16, 16)
        self.Aq1StopBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #ad1813")
        self.Aq2StartBtn = QtWidgets.QPushButton("")
        self.Aq2StartBtn.setFixedSize(16, 16)
        self.Aq2StartBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #13ad2a")
        self.Aq2StopBtn = QtWidgets.QPushButton("")
        self.Aq2StopBtn.setFixedSize(16, 16)
        self.Aq2StopBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #ad1813")
        self.M1StartBtn = QtWidgets.QPushButton("")
        self.M1StartBtn.setFixedSize(16, 16)
        self.M1StartBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #13ad2a")
        self.M1StopBtn = QtWidgets.QPushButton("")
        self.M1StopBtn.setFixedSize(16, 16)
        self.M1StopBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #ad1813")
        self.M2StartBtn = QtWidgets.QPushButton("")
        self.M2StartBtn.setFixedSize(16, 16)
        self.M2StartBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #13ad2a")
        self.M2StopBtn = QtWidgets.QPushButton("")
        self.M2StopBtn.setFixedSize(16, 16)
        self.M2StopBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #ad1813")
        self.solventStartBtn = QtWidgets.QPushButton("")
        self.solventStartBtn.setFixedSize(16, 16)
        self.solventStartBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #13ad2a")
        self.solventStopBtn = QtWidgets.QPushButton("")
        self.solventStopBtn.setFixedSize(16, 16)
        self.solventStopBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #ad1813")
        self.seedStartBtn = QtWidgets.QPushButton("")
        self.seedStartBtn.setFixedSize(16, 16)
        self.seedStartBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #13ad2a")
        self.seedStopBtn = QtWidgets.QPushButton("")
        self.seedStopBtn.setFixedSize(16, 16)
        self.seedStopBtn.setStyleSheet(
            "border-radius : 8;"
            "border : 0.5px solid black;"
            "background-color : #ad1813")
        self.solventValveText = QtWidgets.QLineEdit()
        self.solventValveText.setFixedSize(20, 20)
        self.emulsionValveText = QtWidgets.QLineEdit()
        self.emulsionValveText.setFixedSize(20, 20)
        self.outletValveText = QtWidgets.QLineEdit()
        self.outletValveText.setFixedSize(20, 20)
        self.diaGroupBoxLayout.addWidget(self.diaLabel, 0, 0, 1000, 1000)
        # self.diaGroupBoxLayout.addWidget(self.pumpAq1FlowRateText, 50, 3)
        # self.diaGroupBoxLayout.addWidget(self.pumpAq2FlowRateText, 230, 3)
        # self.diaGroupBoxLayout.addWidget(self.pumpM1FlowRateText, 400, 3)
        # self.diaGroupBoxLayout.addWidget(self.pumpM2FlowRateText, 569, 3)
        # self.diaGroupBoxLayout.addWidget(self.pumpSolventFlowRateText, 875, 72, 1, 2)
        # self.diaGroupBoxLayout.addWidget(self.Aq1StartBtn, 110, 22)
        # self.diaGroupBoxLayout.addWidget(self.Aq1StopBtn, 110, 23)
        # self.diaGroupBoxLayout.addWidget(self.Aq2StartBtn, 280, 22)
        # self.diaGroupBoxLayout.addWidget(self.Aq2StopBtn, 280, 23)
        # self.diaGroupBoxLayout.addWidget(self.M1StartBtn, 454, 22)
        # self.diaGroupBoxLayout.addWidget(self.M1StopBtn, 454, 23)
        # self.diaGroupBoxLayout.addWidget(self.M2StartBtn, 630, 22)
        # self.diaGroupBoxLayout.addWidget(self.M2StopBtn, 630, 23)
        # self.diaGroupBoxLayout.addWidget(self.solventValveText, 760, 50)
        # self.diaGroupBoxLayout.addWidget(self.emulsionValveText, 435, 215)
        # self.diaGroupBoxLayout.addWidget(self.solventStartBtn, 920, 72)
        # self.diaGroupBoxLayout.addWidget(self.solventStopBtn, 920, 73)
        # self.diaGroupBoxLayout.addWidget(self.pumpSeedFlowRateText, 830, 310, 1, 2)
        # self.diaGroupBoxLayout.addWidget(self.seedStartBtn,790, 310)
        # self.diaGroupBoxLayout.addWidget(self.seedStopBtn, 790, 311)
        # self.diaGroupBoxLayout.addWidget(self.outletValveText, 607, 850)

        self.Layout.addWidget(self.diaGroupBox, 1, 1, 1, 2)

        # self.Aq1StartBtn.clicked.connect(self.main.controller.pump4.start())
        # self.Aq1StopBtn.clicked.connect(self.main.controller.pump4.stop())
        # self.Aq2StartBtn.clicked.connect(self.main.controller.pump5.start())
        # self.Aq2StopBtn.clicked.connect(self.main.controller.pump5.stop())
        # self.M1StartBtn.clicked.connect(self.main.controller.pump2.start())
        # self.M1StopBtn.clicked.connect(self.main.controller.pump2.stop())
        # self.M2StartBtn.clicked.connect(self.main.controller.pump3.start())
        # self.M2StopBtn.clicked.connect(self.main.controller.pump3.stop())
        # self.seedStartBtn.clicked.connect(self.main.controller.pump1.start())
        # self.seedStopBtn.clicked.connect(self.main.controller.pump1.stop())
        # self.solventStartBtn.clicked.connect(self.main.controller.pump6.start())
        # self.solventStopBtn.clicked.connect(self.main.controller.pump6.stop())

        # self.pumpAq1FlowRateText.returnPressed.connect(lambda: self.setFlowrate(pump=self.main.controller.pump4, value=self.pumpAq1FlowRateText.text()))
        # self.pumpAq2FlowRateText.returnPressed.connect(lambda: self.setFlowrate(pump=self.main.controller.pump5, value=self.pumpAq2FlowRateText.text()))
        # self.pumpM1FlowRateText.returnPressed.connect(lambda: self.setFlowrate(pump=self.main.controller.pump2, value=self.pumpM1FlowRateText.text()))
        # self.pumpM2FlowRateText.returnPressed.connect(lambda: self.setFlowrate(pump=self.main.controller.pump3, value=self.pumpM2FlowRateText.text()))
        # self.pumpSeedFlowRateText.returnPressed.connect(lambda: self.setFlowrate(pump=self.main.controller.pump1, value=self.pumpSeedFlowRateText.text()))
        # self.pumpSolventFlowRateText.returnPressed.connect(lambda: self.setFlowrate(pump=self.main.controller.pump6, value=self.pumpSolventFlowRateText.text()))

        # self.solventValveText.returnPressed.connect(lambda: self.valveSwitch(valve=self.main.controller.valve1, position=int(self.solventValveText.text())))
        # self.emulsionValveText.returnPressed.connect(lambda: self.valveSwitch(valve=self.main.controller.valve2, position=int(self.emulsionValveText.text())))
        # self.outletValveText.returnPressed.connect(lambda: self.valveSwitch(valve=self.main.controller.valve5, position=int(self.outletValveText.text())))

    def setFlowrate(self, pump, value):
        pump.setFlowrateText.setText(value) 
        pump.setFlowrate()  

    def valveSwitch(self, valve, position):
        if position == 1:
            valve.valveHome()
        else:
            valve.valveSwitch(port=position)

    def maxSeedFracCalculator(self):
        w_eMax = float(self.w_eMaxText.text())
        try: w_fMax = float(self.productConcentrationText.text())
        except(ValueError): pass
        if self.productConcCheckBox.isChecked():
            try: w_fMax = float(self.productConcentrationText2.text())
            except(ValueError): pass
        w_s = float(self.seedConcentrationText.text())
        try: 
            self.R_sMax = round((w_s*(w_eMax - w_fMax))/(w_fMax*(w_eMax - w_s)), 4)
            self.seedAmountText2.setText(str(self.R_sMax))
        except(ZeroDivisionError, UnboundLocalError): pass
        return
    
    def maxSeedFracValidator(self):
        w_eMax = float(self.w_eMaxText.text())
        try: w_fMax = float(self.productConcentrationText.text())
        except(ValueError): pass
        if self.productConcCheckBox.isChecked():
            try: w_fMax = float(self.productConcentrationText2.text())
            except(ValueError): pass
        w_s = float(self.seedConcentrationText.text())
        try: 
            self.R_sMax = round((w_s*(w_eMax - w_fMax))/(w_fMax*(w_eMax - w_s)), 4)
        except(ZeroDivisionError, UnboundLocalError): pass
        R_sValueEntered = float(self.seedAmountText2.text())
        if R_sValueEntered > self.R_sMax:
            self.seedAmountText2.setText(str(self.R_sMax))
        return
    
    def calculateMonomerComposition(self):
        targetTg = float(self.TgText.text()) + 273
        polymer1Tg = float(self.TgPolymer1Text.text()) + 273
        polymer2Tg = float(self.TgPolymer2Text.text()) + 273
        fixedPolymer1Fraction = float(self.fractionFixedCopolymer1Text.text())
        fixedPolymer2Fraction = float(self.fractionFixedCopolymer2Text.text())
        fixedPolymer1Tg = float(self.TgFixedCopolymer1Text.text()) + 273
        fixedPolymer2Tg = float(self.TgFixedCopolymer2Text.text()) + 273

        transitionTemps = [polymer1Tg, polymer2Tg, fixedPolymer1Tg, fixedPolymer2Tg]
        fixedComponentFractions = [fixedPolymer1Fraction, fixedPolymer2Fraction]

        monomerComposition = glassTransitionPredictor.calculateComposition(targetTg, fixedComponentFractions, transitionTemps)
        self.monomerCompositionLabel2.setText(f'A: {monomerComposition[0]}, B: {monomerComposition[1]}, C: {monomerComposition[2]}, D: {monomerComposition[3]}')

    def formatWidget(self, groupBox):
        self.resetWidget(groupBox)
        self.runButton.setText("Start experiment")
        self.reactor1VolumeLabel.setHidden(False)
        self.reactor1VolumeText.setHidden(False)
        self.nCSTRsLabel.setHidden(False)
        self.nCSTRsText.setHidden(False)
        self.reactor2VolumeLabel.setHidden(False)
        self.reactor2VolumeText.setHidden(False)
        self.sampleVolumeLabel.setHidden(False)
        self.sampleVolumeText.setHidden(False)
        self.initiatorStockConcLabel.setHidden(False)
        self.initiatorStockConcText.setHidden(False)
        self.numFeedsText.setHidden(False)
        self.numExpsLabel.setHidden(False)
        self.numExpsText.setHidden(False)
        self.kineticsBox.setHidden(True)
        self.propertyTargetingBox.setHidden(True)
        self.calculateMonomerCompositionBtn.setHidden(True)
        self.numExpsLabel.setText("Number of experiments")

        if self.DoECheckBox.isChecked():
            self.numExpsLabel.setHidden(True)
            self.numExpsText.setHidden(True)

        if self.TSEMOCheckBox.isChecked():
            self.propertyTargetingBox.setHidden(False)
            self.trainingDataBrowseButton.setHidden(False)
            self.trainingDataPathText.setHidden(False)
            self.maxIterationsLabel.setHidden(False)
            self.maxIterationsText.setHidden(False)
            self.sampleStrategyLabel.setHidden(False)
            self.sampleStrategyText.setHidden(False)
            self.numExpsLabel.setText("Number of LHC experiments")
            self.runButton.setText("Start optimisation")
            for checkBox in self.propertyTargetingBox.findChildren(QtWidgets.QCheckBox):
                checkBox.setHidden(False)
            self.particleSizeCheckBox.setHidden(False)
            # for checkBox in self.propertyTargetingBox.findChildren(QtWidgets.QCheckBox):
            #     checkBox.setHidden(False)
            # for label in self.propertyTargetingBox.findChildren(QtWidgets.QLabel):
            #     label.setHidden(False)
            # for lineEdit in self.propertyTargetingBox.findChildren(QtWidgets.QLineEdit):
            #     lineEdit.setHidden(False)

            if self.TgCheckBox.isChecked():
                self.TgText.setHidden(False)
                self.TgPolymer1Text.setHidden(False)
                self.TgPolymer2Text.setHidden(False)
                self.TgFixedCopolymer1Label.setHidden(False)
                self.TgFixedCopolymer1Text.setHidden(False)
                self.TgFixedCopolymer2Label.setHidden(False)
                self.TgFixedCopolymer2Text.setHidden(False)
                self.fractionFixedCopolymer1Label.setHidden(False)
                self.fractionFixedCopolymer1Text.setHidden(False)
                self.fractionFixedCopolymer2Label.setHidden(False)
                self.fractionFixedCopolymer2Text.setHidden(False)
                self.calculateMonomerCompositionBtn.setHidden(False)
                self.monomerCompositionLabel1.setHidden(False)
                self.monomerCompositionLabel2.setHidden(False)
            
            if self.particleSizeCheckBox.isChecked():
                self.particleSizeText.setHidden(False)
                self.particleSizeMappingCheckBox.setHidden(False)

            if self.particleSizeMappingCheckBox.isChecked():
                self.particleSizeText.setHidden(True)

            if self.conventionalEmulsionCheckBox.isChecked():
                self.surfactantConcentrationObjectiveCheckBox.setHidden(False)
                self.seedFractionObjectiveCheckBox.setHidden(False)
                self.initiatorRatioObjectiveCheckBox.setHidden(False)
            
            if self.RAFTEmulsionCheckBox.isChecked():
                self.seedFractionObjectiveCheckBox.setHidden(False)
                self.ctaToInitiatorObjectiveCheckBox.setHidden(False)
                
        if self.OFAATCheckBox.isChecked():
            self.initiatorRatioCheckBox.setAutoExclusive(True)
            self.surfactantRatioCheckBox.setAutoExclusive(True)
            self.monomerRatioCheckBox.setAutoExclusive(True)
            self.seedAmountCheckBox.setAutoExclusive(True)
            self.feedRateCheckBox.setAutoExclusive(True)
            self.residenceTimeCheckBox.setAutoExclusive(True)
            self.DPCheckBox.setAutoExclusive(True)
            self.ctaToInitiatorCheckBox.setAutoExclusive(True)
            self.productConcCheckBox.setAutoExclusive(True)
        else:
            self.surfactantRatioCheckBox.setAutoExclusive(False)
            self.monomerRatioCheckBox.setAutoExclusive(False)
            self.seedAmountCheckBox.setAutoExclusive(False)
            self.feedRateCheckBox.setAutoExclusive(False)
            self.residenceTimeCheckBox.setAutoExclusive(False)
            self.DPCheckBox.setAutoExclusive(False)
            self.ctaToInitiatorCheckBox.setAutoExclusive(False)
            self.productConcCheckBox.setAutoExclusive(False)

        if self.productConcCheckBox.isChecked():
            self.productConcentrationText1.setHidden(False)
            self.productConcentrationText2.setHidden(False)
            self.productConcentrationText.setHidden(True)

        if self.residenceTimeCheckBox.isChecked():
            self.residenceTimeText1.setHidden(False)
            self.residenceTimeText2.setHidden(False)
            self.residenceTimeText.setHidden(True)

        if self.feedRateCheckBox.isChecked():
            self.numFeedsText1.setHidden(False)
            self.numFeedsText2.setHidden(False)
            self.numFeedsText.setHidden(True)

        if self.kineticStudyCheckBox.isChecked():
            self.kineticsBox.setHidden(False)
            for label in self.kineticsBox.findChildren(QtWidgets.QLabel):
                label.setHidden(False)
            for lineEdit in self.kineticsBox.findChildren(QtWidgets.QLineEdit):
                lineEdit.setHidden(False)
            
    ######################################### Conventional emulsion formatting #########################################

        if self.conventionalEmulsionCheckBox.isChecked():
            diagramLoc = r'C:\Users\pm15pmp\Miniconda3\envs\MPSR\polymer_platform\PCubed\Emulsion_FRP.png'
            self.diagram = QtGui.QPixmap(diagramLoc)
            self.diaLabel.setPixmap(self.diagram)
            self.diaLabel.setScaledContents(True)
            self.diaLabel.setFixedSize(1050, 430)
            self.residenceTimeCheckBox.setHidden(False)
            self.residenceTimeText.setHidden(False)
            self.w_eMaxLabel.setHidden(False)
            self.w_eMaxText.setHidden(False)
            self.productConcCheckBox.setHidden(False)
            self.productConcentrationText.setHidden(False)
            self.seedConcentrationLabel.setHidden(False)
            self.seedConcentrationText.setHidden(False)
            self.initiatorRatioCheckBox.setHidden(False)
            self.initiatorRatioText.setHidden(False)
            self.w_initiatorStockLabel.setHidden(False)
            self.w_initiatorStockText.setHidden(False)
            self.monomerDensityLabel.setHidden(False)
            self.monomerDensityText.setHidden(False)
            self.surfactantStockConcLabel.setHidden(False)
            self.surfactantStockConcText.setHidden(False)
            self.surfactantRatioCheckBox.setHidden(False)
            self.surfactantRatioText.setHidden(False)
            self.w_surfactantStockLabel.setHidden(False)
            self.w_surfactantStockText.setHidden(False)
            self.monomerRatioCheckBox.setHidden(False)
            self.monomerRatioText.setHidden(False)
            self.monomerDensityLabel.setHidden(False)
            self.monomerDensityText.setHidden(False)
            self.seedAmountCheckBox.setHidden(False)
            self.seedAmountText.setHidden(False)
            self.feedRateCheckBox.setHidden(False)
            self.numFeedsText.setHidden(False)
        else:
            self.initiatorRatioCheckBox.setChecked(False)
            self.surfactantRatioCheckBox.setChecked(False)
            self.monomerRatioCheckBox.setChecked(False)
            self.seedAmountCheckBox.setChecked(False)

        if self.initiatorRatioCheckBox.isChecked():
            self.initiatorRatioText1.setHidden(False)
            self.initiatorRatioText2.setHidden(False)
            self.initiatorRatioText.setHidden(True)

        if self.surfactantRatioCheckBox.isChecked():
            self.surfactantRatioText1.setHidden(False)
            self.surfactantRatioText2.setHidden(False)
            self.surfactantRatioText.setHidden(True)

        if self.monomerRatioCheckBox.isChecked():
            self.monomerRatioText1.setHidden(False)
            self.monomerRatioText2.setHidden(False)
            self.monomerDensityText1.setHidden(False)
            self.monomerDensityText2.setHidden(False)
            self.monomerRatioText.setHidden(True)
            self.monomerDensityLabel.setHidden(True)
            self.monomerDensityText.setHidden(True)

        if self.seedAmountCheckBox.isChecked():
            self.seedAmountText1.setHidden(False)
            self.seedAmountText2.setHidden(False)
            self.w_eMaxText.setHidden(False)
            self.seedAmountText.setHidden(True)

    ####################################### RAFT emulsion formatting ##########################################

        if self.RAFTEmulsionCheckBox.isChecked():
            diagramLoc = r'C:\Users\pm15pmp\Miniconda3\envs\MPSR\polymer_platform\PCubed\RAFT_platform.png'
            self.diagram = QtGui.QPixmap(diagramLoc)
            self.diaLabel.setPixmap(self.diagram)
            self.diaLabel.setScaledContents(True)
            self.diaLabel.setFixedSize(1100, 450)
            self.residenceTimeCheckBox.setHidden(False)
            self.residenceTimeText.setHidden(False)
            self.macroCTAStockConcLabel.setHidden(False)
            self.macroCTAStockConcText.setHidden(False)
            self.seededQueryCheckBox.setHidden(False)
            if self.seededQueryCheckBox.isChecked():
                self.macroCTAPolymerConcLabel.setHidden(False)
                self.macroCTAPolymerConcText.setHidden(False)
                self.macroCTAMonomerConcLabel.setHidden(False)
                self.macroCTAMonomerConcText.setHidden(False)
                self.macroCTAStockConcLabel.setHidden(True)
                self.macroCTAStockConcText.setHidden(True)
            self.DPCheckBox.setHidden(False)
            self.DPText.setHidden(False)
            self.ctaToInitiatorCheckBox.setHidden(False)
            self.ctaToInitiatorText.setHidden(False)
            self.productConcCheckBox.setHidden(False)
            self.productConcentrationText.setHidden(False)
            self.molarMassMonomerLabel.setHidden(False)
            self.molarMassMonomerText.setHidden(False)
            self.molarMassInitiatorLabel.setHidden(False)
            self.molarMassInitiatorText.setHidden(False)
            self.molarMassMacroCTALabel.setHidden(False)
            self.molarMassMacroCTAText.setHidden(False)
            self.monomerDensityLabel.setHidden(False)
            self.monomerDensityText.setHidden(False)
            self.feedRateCheckBox.setHidden(False)
            self.numFeedsText.setHidden(False)
            
        else:
            self.DPCheckBox.setChecked(False)
            self.ctaToInitiatorCheckBox.setChecked(False)

        if self.DPCheckBox.isChecked():
            self.DPText1.setHidden(False)
            self.DPText2.setHidden(False)
            self.DPText.setHidden(True)

        if self.ctaToInitiatorCheckBox.isChecked():
            self.ctaToInitiatorText1.setHidden(False)
            self.ctaToInitiatorText2.setHidden(False)
            self.ctaToInitiatorText.setHidden(True)

    def resetWidget(self, groupBox):
        '''Function used to hide all widgets in a specified groupbox for user convenience'''
        for label in groupBox.findChildren(QtWidgets.QLabel):
            label.setHidden(True)
        for lineEdit in groupBox.findChildren(QtWidgets.QLineEdit):
            lineEdit.setHidden(True)
        for checkBox in groupBox.findChildren(QtWidgets.QCheckBox):
            checkBox.setHidden(True)

    def buildExperiment(self):

        '''
        This is the main function for defining the experiment to run, which gathers all the values from user inputs
        and sends them to the appropriate calculator depending on the type of experiment and type of polymerisation selected.
        
        Where values are not specified by the user, a default value is usually provided - check your inputs carefully.
        
        In a nutshell, most of this function just organises the experimental parameters (reactor volumes, residence times etc.)
        into a dictionary of {variable : [values]} which then gets a dictionary of {pump : [flowrates]} in an experimental design (expDesign).

        The length of the [values] lists are defined by the number of pre-programmed experiments requested (i.e. number of reactions in an
        automated screen such as an OFAAT or DoE). For optimisation experiments, a set of pre-programmed experiments will be prepared based
        on latin hypercube sampling before the optimisation sequence begins. Alternatively, training data (from a model or previous experiment) 
        can be specified, and this function will use the optimiser to 'pre-programme' the first iteration of the optimisation from it.

        Calculators are given lists of every variable necessary for making the flow rate calculations. Where these variables are fixed, the
        lists will contain the same value repeated for the number of experiments specified.
        
        '''

        # Make sure a method has been selected:
        if not (self.OFAATCheckBox.isChecked() or self.DoECheckBox.isChecked() or self.TSEMOCheckBox.isChecked() or self.kineticStudyCheckBox.isChecked()):
            print('Please select a method')

        ### Catch issues with more than one variable being selected for an OFAAT experiment
        variableCounter = 0
        for checkBox in self.variableSelectorBox.findChildren(QtWidgets.QCheckBox):
            if checkBox.isChecked():
                variableCounter += 1
        if self.seededQueryCheckBox.isChecked():
            variableCounter += -1
        if self.OFAATCheckBox.isChecked() and variableCounter > 1:
            print("Too many variables selected for OFAAT experiment - please just choose one")
            return

        '''First collect all the experimental parameters from user inputs - error handling to
        ignore problems with trying to read default text like "min" and "max" in hidden lineEdits'''

        # Get common fixed parameters from user inputs #
        volumeCSTR = [float(self.reactor1VolumeText.text())] # Volume of a single CSTR
        numberCSTRs = [int(self.nCSTRsText.text())] # Number of CSTRs in cascade
        volumeReactor2 = [float(self.reactor2VolumeText.text())]
        deadVolume = [1] # Volume of tubing between reactor outlet and sampling valve [mL]
        volumeSample = [float(self.sampleVolumeText.text())] # Specify amount of sample to collect
        initiatorStockConc = [float(self.initiatorStockConcText.text())] # Concentration of only initiator in initiator stock solution [w/w]
        w_initiatorStock = [float(self.w_initiatorStockText.text())] # Solids content of all components in initiator stock solution (buffer etc.) [w/w]
        densityAq = [1] # Density of aqueous solution (assuming water)
        densityMonomer = [float(self.monomerDensityText.text())] # Density of monomer if composition is kept constant [g/mL]
        densityProduct = [1] # Density of final product (assuming density of water for now)
        numExp = int(self.numExpsText.text()) # Number of experiments to run/samples to take
        cleanRate = [float(self.cleaningRateText.text())] # Flow rate for solvent washing
        cleanTime = [float(self.cleaningVolumesText.text())*((volumeCSTR[0]*numberCSTRs[0]) + volumeReactor2[0])/cleanRate[0]] # Time taken to clean between reactions (min)
        cleanVolume = [cleanRate[0]*cleanTime[0]] # Volume of material used for each cleaning step (mL)
        flushTime = 2*((volumeCSTR[0]*numberCSTRs[0]) + volumeReactor2[0])/cleanRate[0] # Time taken to flush reactor after experiment

        # Get common variable parameters from user inputs #
        w_f = [float(self.productConcentrationText.text())] # Specify final product concentration [w/w]
        try:
            w_fMin = [float(self.productConcentrationText1.text())] # Minimum product concentration [w/w]
            w_fMax = [float(self.productConcentrationText2.text())] # Maximum product concentration [w/w]
        except(ValueError): pass

        tau = [float(self.residenceTimeText.text())] # Residence time - for seeded processes this is the average RT of a single seed particle [min]
        try:
            tauMin = [float(self.residenceTimeText1.text())] # Minimum residence time [min]
            tauMax = [float(self.residenceTimeText2.text())] # Maximum residence time [min]
        except(ValueError): pass

        initiatorRatio = [float(self.initiatorRatioText.text())] # Initiator concentration in emulsion if amount is kept constant [w/w]
        try:
            initiatorRatioMin = [float(self.initiatorRatioText1.text())] # Minimum bound of initiator concentration in emulsion [w/w]
            initiatorRatioMax = [float(self.initiatorRatioText2.text())] # Maximum bound of initiator concentration in emulsion [w/w]
        except(ValueError): pass
        densityInitiator = [1] # Density of initiator solution [g/mL]
        
        monomerAFraction = [float(self.monomerRatioText.text())]
        try:
            monomerAFractionMin = [float(self.monomerRatioText1.text())] # Minimum fraction [w/w] of monomer 'A' (value between 0.0 and 1.0)
            monomerAFractionMax = [float(self.monomerRatioText2.text())] # Maximum fraction [w/w] of monomer 'A' (value between 0.0 and 1.0)
            densityMonomerA = [float(self.monomerDensityText1.text())] # Density of monomer 'A' [g/mL]
            densityMonomerB = [float(self.monomerDensityText2.text())] # Density of monomer 'B' [g/mL]
        except(ValueError): 
            densityMonomerA = [float(self.monomerDensityText.text())]
            densityMonomerB = [0]

        # Get fixed parameters specific to conventional emulsion polymerisation #
        w_s = [float(self.seedConcentrationText.text())] # Specify seed concentration [w/w]
        w_eMax = [float(self.w_eMaxText.text())] # Maximum permissible emulsion concentration
        surfactantStockConc = [float(self.surfactantStockConcText.text())] # Concentration of only surfactant in surfactant stock solution [w/w]
        w_surfactantStock = [float(self.w_surfactantStockText.text())] # Solids content of all components in surfactant stock solution (surfactant, buffer etc.) [w/w]
        w_makeupWater = [float(self.w_makeupWaterText.text())] # Solids content of all components in make-up water (buffer etc.) [w/w]
        densitySeed = [1] # Density of latex seed (assuming density of water, can be measured before exp)

        # Get variable parameters specific to conventional emulsion polymerisation #
        surfactantRatio = [float(self.surfactantRatioText.text())] # Surfactant concentration in emulsion if amount is kept constant [w/w]
        
        try:
            surfactantRatioMin = [float(self.surfactantRatioText1.text())] # Minimum bound of surfactant concentration in emulsion [w/w]
            surfactantRatioMax = [float(self.surfactantRatioText2.text())] # Maximum bound of surfactant concentration in emulsion [w/w]
        except(ValueError): pass

        densitySurfactant = [1] # Density of surfactant solution [g/mL]
        seedFrac = [float(self.seedAmountText.text())] # Seed fraction if kept constant (ratio of seed solids to total solids content of product)
        try:
            seedFracMin = [float(self.seedAmountText1.text())] # Minimum seed fraction (ratio of seed solids to total solids content of product)
            seedFracMax = [float(self.seedAmountText2.text())] # Maximum seed fracion (ratio of seed solids to total solids content of product)
        except(ValueError): pass    

        numFeeds = [int(self.numFeedsText.text())] # Number of emulsion feeds if kept constant
        try:
            numFeedsMin = [int(self.numFeedsText1.text())] # Minimum number of emulsion feeds (proxy for emulsion feed rate)
            numFeedsMax = [int(self.numFeedsText2.text())] # Maximum number of emulsion feeds (proxy for emulsion feed rate)
        except(ValueError): pass

        # Get  parameters specific to RAFT emulsion polymerisation #
        densityMacroCTA = [1]
        try:
            macroCTAStockConc = [float(self.macroCTAStockConcText.text())] # Concentrtion of macro-CTA in stock solution [w/w]
            macroCTAPolymerConc = [float(self.macroCTAPolymerConcText.text())] # Concentration of polymer only in the macro-CTA solution (usually for seeded reactions) [w/w]
            macroCTAMonomerConc = [float(self.macroCTAMonomerConcText.text())] # Concentration of the residual monomer in the macro-CTA feed (usually for seeded reactions) [w/w]
            molarMassMacroCTA = [float(self.molarMassMacroCTAText.text())]# Molar mass of macro-CTA
            molarMassMonomer = [float(self.molarMassMonomerText.text())]# Molar mass of monomer
            molarMassInitiator = [float(self.molarMassInitiatorText.text())]# Molar mass of initiator
        except(ValueError): pass

        # Get variable parameters specific to RAFT emulsion polymerisation #
        DP = [float(self.DPText.text())]# Target degree of polymerisation if kept constant [moles of monomer/moles of macro-CTA]
        try:
            DP_Min = [float(self.DPText1.text())]# Minimum bound of block copolymer degree of polymerisation (DP) [moles of monomer/moles of macro-CTA]
            DP_Max = [float(self.DPText2.text())]# Maximum bound of block copolymer degree of polymerisation (DP) [moles of monomer/moles of macro-CTA]
        except(ValueError): pass
        ctaToInitiator = [float(self.ctaToInitiatorText.text())]# CTA:initiator ratio if kept constant [moles mcaro-CTA/moles initiator]
        try:
            ctaToInitiatorMin = [float(self.ctaToInitiatorText1.text())]# Minimum bound of CTA:initiator ratio [moles macro-CTA/moles initiator]
            ctaToInitiatorMax = [float(self.ctaToInitiatorText2.text())]# Maximum bound of CTA:initiator ratio [moles macro-CTA/moles initiator]
        except(ValueError): pass

        '''
        Next part of the function collects all the selected parameters (names 
        and bounds) into lists ready to be handled by the chosen method
        '''

        self.parameterNames = [] # Blank list to hold the names (str) of the varied parameters involved in the experiment
        self.parameterBounds = [] # Blank list to hold the min/max of each selected parameter [min, max]

        # If OFAAT method selected, build an index list to loop over based on the number of experiments required
        if self.OFAATCheckBox.isChecked():
            OFAAT_dict = {}
            expIndex = np.linspace(1, numExp, num=numExp)
            OFAAT_dict["Experiment"] = expIndex

        # kinetics study currently set up the same as OFAAT, so can experiment with different feed rates and residence times
        elif self.kineticStudyCheckBox.isChecked():
            kinetics_dict = {}
            expIndex = np.linspace(1, numExp, num=numExp)
            kinetics_dict["Experiment"] = expIndex

        # For each possible variable, check if it has been selected and either a) run a OFAAT calculation or b) add the name 
        # and bounds to the experiment builder lists above

        if self.productConcCheckBox.isChecked():
            productConcBounds = w_fMin + w_fMax
            self.parameterNames.append("Product_concentration")
            self.parameterBounds.append(productConcBounds)
            if self.OFAATCheckBox.isChecked():
                OFAAT_dict["Product_concentration"] = np.linspace(w_fMin[0], w_fMax[0], numExp)
        
        if self.residenceTimeCheckBox.isChecked():
            residenceTimeBounds = tauMin + tauMax
            self.parameterNames.append("Residence_time")
            self.parameterBounds.append(residenceTimeBounds)
            if self.OFAATCheckBox.isChecked():
                OFAAT_dict["Residence_time"] = np.linspace(tauMin[0], tauMax[0], numExp)
            if self.kineticStudyCheckBox.isChecked():
                kinetics_dict["Residence_time"] = np.linspace(tauMin[0], tauMax[0], numExp)

        if self.initiatorRatioCheckBox.isChecked():
            initiatorRatioBounds = initiatorRatioMin + initiatorRatioMax
            self.parameterNames.append("Initiator_ratio")
            self.parameterBounds.append(initiatorRatioBounds)
            if self.OFAATCheckBox.isChecked():
                OFAAT_dict["Initiator_ratio"] = np.linspace(initiatorRatioMin[0], initiatorRatioMax[0], numExp)

        if self.surfactantRatioCheckBox.isChecked():
            surfactantRatioBounds = surfactantRatioMin + surfactantRatioMax
            self.parameterNames.append("Surfactant_ratio")
            self.parameterBounds.append(surfactantRatioBounds)
            if self.OFAATCheckBox.isChecked():
                OFAAT_dict["Surfactant_ratio"] = np.linspace(surfactantRatioMin[0], surfactantRatioMax[0], numExp)

        if self.monomerRatioCheckBox.isChecked():
            monomerRatioBounds = monomerAFractionMin + monomerAFractionMax
            self.parameterNames.append("Monomer_A_fraction")
            self.parameterBounds.append(monomerRatioBounds)
            if self.OFAATCheckBox.isChecked():
                OFAAT_dict["Monomer_A_fraction"] = np.linspace(monomerAFractionMin[0], monomerAFractionMax[0], numExp)

        if self.seedAmountCheckBox.isChecked():
            seedAmountBounds = seedFracMin + seedFracMax
            self.parameterNames.append("Seed_fraction")
            self.parameterBounds.append(seedAmountBounds)
            if self.OFAATCheckBox.isChecked():
                OFAAT_dict["Seed_fraction"] = np.linspace(seedFracMin[0], seedFracMax[0], numExp)

        if self.feedRateCheckBox.isChecked():
            emulsionFeedsBounds = numFeedsMin + numFeedsMax
            self.parameterNames.append("Num_feeds")
            self.parameterBounds.append(emulsionFeedsBounds)
            if self.OFAATCheckBox.isChecked():
                OFAAT_dict["Num_feeds"] = np.linspace(numFeedsMin[0], numFeedsMax[0], numExp).astype(int)
            if self.kineticStudyCheckBox.isChecked():
                kinetics_dict["Num_feeds"] = np.linspace(numFeedsMin[0], numFeedsMax[0], numExp).astype(int)
            
            print(f'kinetics_dict: {kinetics_dict}')
        
        if self.DPCheckBox.isChecked():
            DP_bounds = DP_Min + DP_Max
            self.parameterNames.append("DP")
            self.parameterBounds.append(DP_bounds)
            if self.OFAATCheckBox.isChecked():
                OFAAT_dict["DP"] = np.linspace(DP_Min[0], DP_Max[0], numExp)

        if self.ctaToInitiatorCheckBox.isChecked():
            ctaToInitiatorBounds = ctaToInitiatorMin + ctaToInitiatorMax
            self.parameterNames.append("CTA_to_initiator")
            self.parameterBounds.append(ctaToInitiatorBounds)
            if self.OFAATCheckBox.isChecked():
                OFAAT_dict["CTA_to_initiator"] = np.linspace(ctaToInitiatorMin[0], ctaToInitiatorMax[0], numExp)

        '''
        Now pass the lists of varied parameters to the relevant experiment builder (DoE or LHS) to return a
        dataFrame containing experimental conditions for a series of pre-programmed experiments
        '''
        if self.OFAATCheckBox.isChecked():
            expDesign = pd.DataFrame.from_dict(OFAAT_dict)

        elif self.DoECheckBox.isChecked():
            if len(self.parameterNames) == 1:
                print("Please select more variables for DoE experiment")
                return
            expDesign = self.DoEHandler.buildCC(bounds=self.parameterBounds, variables=self.parameterNames)

        elif self.TSEMOCheckBox.isChecked():
            '''
            Builds the optimisation experiment by collecting selected input variables and bounds in a dictionary {variable : bounds}
            and then takes the selected optimisation objectives also in a dictionary, but with max/min also specified - {(objective, strategy) : bounds}.
            'Bounds' don't really matter for the objectives, but the optimiser function requires it, so just use the bounds specified in the inputs, or
            some arbitrary broad bounds (as for the particle size objective)
            '''
            self.optIterations = 1 # Tracks the total number of experimental iterations (including training data)
            self.maxIterations = int(self.maxIterationsText.text())
            self.inputsDict = {}
            for variable, bounds in zip(self.parameterNames, self.parameterBounds):
                self.inputsDict[variable] = bounds

            self.objectivesDict = {} # Build a dictionary for passing objective variables to the algorithm
            if self.particleSizeCheckBox.isChecked():
                particleSizeBounds = [0, 1000] # Arbitrary broad bounds of particle size
                if self.particleSizeMappingCheckBox.isChecked():
                    strategy = "FLIPFLOP" # 'FLIPFLOP' strategy is employed when we want to map particle sizes whilst min/max-ing other objectives, this alternates between minimising and maximising the particle size objective
                    self.sizeKey = "Particle_size"
                    self.sizeTarget=None
                else:
                    strategy = "MIN"
                    self.sizeTarget = self.particleSizeText.text()
                    self.sizeKey = "Size_function" # For targeting a specific size we need to convert particle size to a function we can minimise, not sure exactly how to handle this for now
                self.objectivesDict[(self.sizeKey, strategy)] = particleSizeBounds
            
            if self.surfactantConcentrationObjectiveCheckBox.isChecked():
                surfactantConcBounds = [surfactantRatioMax[0], surfactantRatioMax[0]]
                self.objectivesDict[("Surfactant_ratio_min", "MIN")] = surfactantConcBounds # Always want to minimise surfactant concentration by default

            if self.seedFractionObjectiveCheckBox.isChecked():
                seedFractionBounds = [seedFracMin[0], seedFracMax[0]]
                self.objectivesDict[("Seed_fraction_min", "MIN")] = seedFractionBounds # Always want to minimise seed fraction by default

            if self.initiatorRatioObjectiveCheckBox.isChecked():
                initiatorRatioBounds = [initiatorRatioMin[0], initiatorRatioMax[0]]
                self.objectivesDict[("Initiator_ratio_min", "MIN")] = initiatorRatioBounds # Always want to minimise seed fraction by default
        
            allVariables=[] # Put all inputs and objectives in a list, used for checking that the training data provides all the required information
            for var in list(self.inputsDict.keys()):
                allVariables.append(var)
            for obj in list(self.objectivesDict.keys()):
                allVariables.append(obj[0])
            
            self.trainingDataQ=False
            if not self.trainingDataPathText.text()=='': # If training data lineEdit is not blank (i.e. training data has been given)
                self.trainingDataQ=True # Update trainingDataQuery to true so the optimiser knows to use training data
                self.trainingDataPath = self.trainingDataPathText.text()
                countFiles = 0
                for file in glob.glob(self.trainingDataPath + r"\*.xlsx"): # Counts all files in the training data directory, must be exactly 1 excel sheet of training data
                    countFiles+=1
                self.nonExclusiveTrainingDataMsgBox = QMessageBox()
                self.nonExclusiveTrainingDataMsgBox.setStandardButtons(QMessageBox.Ok)
                if countFiles > 1:
                    self.nonExclusiveTrainingDataMsgBox.setWindowTitle("Too many files in specified directory")
                    self.nonExclusiveTrainingDataMsgBox.setText("Tried to read more than one file in the specified directory, please\nput the training data in its own folder and reselect the path")
                    self.nonExclusiveTrainingDataMsgBox.exec()
                elif countFiles < 1:
                    self.nonExclusiveTrainingDataMsgBox.setWindowTitle("No files found in current directory")
                    self.nonExclusiveTrainingDataMsgBox.setText("No files found in the specified directory, please check the location\nof your training data and that it is .xlsx format")
                    self.nonExclusiveTrainingDataMsgBox.exec()
                if self.nonExclusiveTrainingDataMsgBox.clickedButton() is self.nonExclusiveTrainingDataMsgBox.button(QMessageBox.Ok):
                    return
                else:
                    self.trainingDatadf = pd.read_excel(file)
                    self.trainingDataDict = self.trainingDatadf.to_dict(orient='list') # convert the training data excel file to a dictionary
                self.optIterations = self.trainingDatadf.shape[0] # Update the total number of experimental iterations to include the training data

                # Check that all the variables in the training data match the selected inputs and objectives of the current experiment
                for key_to_check in allVariables:
                    try:
                        check = self.trainingDataDict[key_to_check]
                    except KeyError:
                        self.missingVariableMsgBox = QMessageBox()
                        self.missingVariableMsgBox.setStandardButtons(QMessageBox.Ok)
                        self.missingVariableMsgBox.setWindowTitle("Selected variable does not exist in training data")
                        self.missingVariableMsgBox.setText("Either an input or objective variable is missing from the training data, please match\nthe training data with the inputs and objectives you are using")
                        self.missingVariableMsgBox.exec()
                        if self.missingVariableMsgBox.clickedButton() is self.missingVariableMsgBox.button(QMessageBox.Ok):
                            return
                expDesign = optimiser.getNextExperiment(self.inputsDict, 
                                                        self.objectivesDict,
                                                        self.optIterations,
                                                        current_dataSet=self.trainingDataDict)
            
            # Get LHC sampling experiments or request training data
            if int(self.numExpsText.text()) > 0:
                self.LHCQ = True # Set LHCQ true if the user has specified to run LHC experiments
                self.expDesignDict = latinHypercubeSampling.getHypercubeSamplingParams(numExp, self.inputsDict)
                expDesign = pd.DataFrame(self.expDesignDict)
            elif not self.trainingDataQ:
                self.LHCQ = False # Set LHCQ false if the user has NOT specified to run LHC experiments
                self.noLHCMsgBox = QMessageBox()
                self.noLHCMsgBox.setWindowTitle("No LHC experiments planned")
                self.noLHCMsgBox.setText("You have selected to run zero LHC experiments, either provide your own training data\nor continue and the algorithm will run its own LHC")
                self.noLHCMsgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
                self.noLHCMsgBox.setDefaultButton(QMessageBox.Cancel)
                self.noLHCMsgBox.exec()
                if self.noLHCMsgBox.standardButton(self.noLHCMsgBox.clickedButton()) == QMessageBox.Ok:
                    expDesign = optimiser.getNextExperiment(self.inputsDict, 
                                                            self.objectivesDict, 
                                                            iteration=self.optIterations)
                if self.noLHCMsgBox.clickedButton() is self.noLHCMsgBox.button(QMessageBox.Cancel):
                    return

            if self.trainingDataQ: # If training data is provided, add the expDesign (from initial conditions or LHC) to it
                self.optimisationSummary = self.trainingDataDict # Store the training data in a new 'optimisation summmary' which holds the input and objective variables
                for key in expDesign:
                    self.optimisationSummary[key].extend(expDesign[key]) # Append the new input variables from the expDesign to the optimisation summary
            else: # Otherwise, create a blank optimisation summary and fill it in with the expDesign
                self.optimisationSummary = {}
                for var in allVariables:
                    self.optimisationSummary[var] = []
                for key in expDesign:
                    if key in self.optimisationSummary:
                        self.optimisationSummary[key].extend(expDesign[key])
                    
            expDesign = pd.DataFrame(expDesign)

        elif self.kineticStudyCheckBox.isChecked():
            expDesign = pd.DataFrame.from_dict(kinetics_dict)

        '''Use the experimental design to build lists for each variable with n elements based on the
        'shape' of expDesign. Having each variable as a list (even a list of the same repeating value)
        allows the pump flow rates to be calculated in the calculateFlowrates function for any future
        combination of varied experimental parameters'''

        numExp = pd.DataFrame.from_dict(expDesign).shape[0]
        self.expIndex = np.linspace(0, numExp - 1, numExp)
        self.variablesDict = {}
        for n in range(0, numExp-1):
            ### Core variables ###
            volumeCSTR.append(volumeCSTR[0])
            numberCSTRs.append(numberCSTRs[0])
            volumeReactor2.append(volumeReactor2[0])
            deadVolume.append(deadVolume[0])
            volumeSample.append(volumeSample[0])
            initiatorStockConc.append(initiatorStockConc[0])
            w_initiatorStock.append(w_initiatorStock[0])
            densityInitiator.append(densityInitiator[0])
            densityAq.append(densityAq[0])
            densityProduct.append(densityProduct[0])
            tau.append(tau[0])
            w_f.append(w_f[0])
            initiatorRatio.append(initiatorRatio[0])
            monomerAFraction.append(monomerAFraction[0])
            densityMonomer.append(densityMonomer[0])
            densityMonomerA.append(densityMonomerA[0])
            densityMonomerB.append(densityMonomerB[0])
            cleanRate.append(cleanRate[0])
            cleanTime.append(cleanTime[0])
            cleanVolume.append(cleanVolume[0])
            
            ### Conventional emulsion polymerisation variables ###
            w_s.append(w_s[0])
            surfactantStockConc.append(surfactantStockConc[0])
            w_surfactantStock.append(w_surfactantStock[0])
            densitySurfactant.append(densitySurfactant[0])
            w_makeupWater.append(w_makeupWater[0])
            densitySeed.append(densitySeed[0])
            surfactantRatio.append(surfactantRatio[0])
            seedFrac.append(seedFrac[0])
            numFeeds.append(numFeeds[0])

            ### RAFT emulsion polymerisation variables ###
            densityMacroCTA.append(densityMacroCTA[0])
            macroCTAStockConc.append(macroCTAStockConc[0])
            macroCTAPolymerConc.append(macroCTAPolymerConc[0])
            macroCTAMonomerConc.append(macroCTAMonomerConc[0])
            molarMassMacroCTA.append(molarMassMacroCTA[0])
            molarMassMonomer.append(molarMassMonomer[0])
            molarMassInitiator.append(molarMassInitiator[0])
            DP.append(DP[0])
            ctaToInitiator.append(ctaToInitiator[0])

        '''Now build a dict of parameter key-words and their list of values, and recreate the lists of 
        those parameters that are varied by substituting them to the values specified by the experimental design'''
        
        self.variablesDict = {
            'volumeCSTR' : volumeCSTR,
            'numberCSTRs' : numberCSTRs,
            'volumeReactor2' : volumeReactor2,
            'deadVolume' : deadVolume,
            'volumeSample' : volumeSample,
            'initiatorStockConc' : initiatorStockConc,
            'w_initiatorStock' : w_initiatorStock,
            'densityInitiator' : densityInitiator,
            'densityAq' : densityAq,
            'densityProduct' : densityProduct,
            'tau' : tau,
            'w_f' : w_f,
            'Initiator_ratio' : initiatorRatio,
            'monomerAFraction' : monomerAFraction,
            'densityMonomer' : densityMonomer,
            'densityMonomerA' : densityMonomerA,
            'densityMonomerB' : densityMonomerB,
            'w_s' : w_s,
            'surfactantStockConc' : surfactantStockConc,
            'w_surfactantStock' : w_surfactantStock,
            'densitySurfactant' : densitySurfactant,
            'w_makeupWater' : w_makeupWater,
            'densitySeed' : densitySeed,
            'Surfactant_ratio' : surfactantRatio,
            'Seed_fraction' : seedFrac,
            'numFeeds' : numFeeds,
            'macroCTAStockConc' : macroCTAStockConc,
            'macroCTAPolymerConc' : macroCTAPolymerConc,
            'macroCTAMonomerConc' : macroCTAMonomerConc,
            'densityMacroCTA' : densityMacroCTA,
            'molarMassMacroCTA' : molarMassMacroCTA,
            'molarMassMonomer' : molarMassMonomer,
            'molarMassInitiator' : molarMassInitiator,
            'DP' : DP,
            'ctaToInitiator' : ctaToInitiator,
            'cleanRate' : cleanRate,
            'cleanTime' : cleanTime,
            'cleanVolume' : cleanVolume
        }

        '''Substitute the 'default' values of the manipulated variables for the lists defined by the experimental design'''

        for columnName in expDesign:

            if columnName == "Random state value": # Collect the random state value for building the latin hypercube and then drop it from the DF. This will make it accessible in the experiment summary but not interfere with plotting/formatting
                randomStateValue = expDesign[columnName].to_list()
                self.variablesDict["Random state value"] = randomStateValue
                expDesign.drop(columns="Random state value", inplace=True, axis=1)

            if columnName == 'Product_concentration':
                w_f = expDesign[columnName].to_list()
                self.variablesDict["w_f"] = list(np.round(w_f, 3))

            if columnName == 'Residence_time':
                tau = expDesign[columnName].to_list()
                self.variablesDict["tau"] = list(np.round(tau, 3))

            if columnName == 'Initiator_ratio':
                initiatorRatio = expDesign[columnName].to_list()
                self.variablesDict["Initiator_ratio"] = list(np.round(initiatorRatio, 3))

            if columnName == 'Surfactant_ratio':
                surfactantRatio = expDesign[columnName].to_list()
                self.variablesDict["Surfactant_ratio"] = list(np.round(surfactantRatio, 3))

            if columnName == 'Monomer_A_fraction':
                monomerAFraction = expDesign[columnName].to_list()
                self.variablesDict["monomerAFraction"] = list(np.round(monomerAFraction, 3))

                densityMonomer = [] # Calculate sequence of monomer densities as the weighted average of the two monomers in their defined ratio
                for i in monomerAFraction:
                    densityMonomer_temp = np.round((i*densityMonomerA[0] + (1 - i)*densityMonomerB[0]), 3)
                    densityMonomer.append(densityMonomer_temp)
                
                self.variablesDict["densityMonomerA"] = densityMonomerA
                self.variablesDict["densityMonomerB"] = densityMonomerB
                self.variablesDict["densityMonomer"] = densityMonomer

            if columnName == 'Seed_fraction':
                seedFrac = expDesign[columnName].to_list()
                self.variablesDict["Seed_fraction"] = list(np.round(seedFrac, 3))

            if columnName == 'Num_feeds':
                numFeeds = expDesign[columnName].to_list()
                self.variablesDict["numFeeds"] = numFeeds

            if columnName == 'DP':
                DP = expDesign[columnName].to_list()
                self.variablesDict["DP"] = list(np.round(DP, 0))

            if columnName == 'CTA_to_initiator':
                ctaToInitiator = expDesign[columnName].to_list()
                self.variablesDict["ctaToInitiator"] = list(np.round(ctaToInitiator, 3))

        '''Now pass the dict to the appropriate calculator'''

        if self.conventionalEmulsionCheckBox.isChecked():
            self.flowRates = flowrateCalculator.calculateFlowrates(
                'conventional',
                self.variablesDict,
                numExp, 
                )
            self.main.controller.pump1.setHidden(False)
            self.main.controller.pump2.setHidden(False)
            self.main.controller.pump6.setHidden(False)
            self.main.controller.pump9.setHidden(False)
            self.main.controller.pump10.setHidden(False)
            self.main.controller.pump11.setHidden(False)
            self.main.controller.pump12.setHidden(False)

        if self.RAFTEmulsionCheckBox.isChecked():
            if self.seededQueryCheckBox.isChecked():
                self.flowRates = flowrateCalculator.calculateFlowrates(
                    'RAFT',
                    self.variablesDict,
                    numExp,
                    seeded=True
                )
            else:
                self.flowRates = flowrateCalculator.calculateFlowrates(
                    'RAFT',
                    self.variablesDict,
                    numExp
                    )
            self.main.controller.pump2.setHidden(False)
            self.main.controller.pump6.setHidden(False)
            self.main.controller.pump9.setHidden(False)
            self.main.controller.pump10.setHidden(False)
            self.main.controller.pump11.setHidden(False)
            self.main.controller.pump13.setHidden(False)

        self.variablesDict.update(self.flowRates)
        
        '''After creating the experimental design, the following code is used to display the table and graph in the GUI'''

        if self.conventionalEmulsionCheckBox.isChecked():
            self.parametersTable.setRowCount(11)
            self.parametersTable.setColumnCount(numExp)
            self.parametersTable.setVerticalHeaderLabels([
                'TSC (g/g)',
                'Residence time (min)',
                'Initiator ratio (pphm)',
                'Surfactant ratio (pphm)', 
                'Seed fraction (g/g)',
                'Emulsion feeds',
                'Seed flow rate (mL/min)',
                'Surfactant flow rate (mL/min)',
                'Initiator flow rate (mL/min)',
                'Monomer flow rate (mL/min)',
                'Makeup water flow rate (mL/min)'
                ])
            
            for i in self.expIndex.astype(int):
                tableTSC = QtWidgets.QTableWidgetItem(str(round(w_f[i], 4)))
                tableRT = QtWidgets.QTableWidgetItem(str(round(tau[i], 4)))
                tableinitRatio = QtWidgets.QTableWidgetItem(str(round(initiatorRatio[i], 4)))
                tablesurfactantRatio = QtWidgets.QTableWidgetItem(str(round(surfactantRatio[i], 4)))
                tableSeedFrac = QtWidgets.QTableWidgetItem(str(round(seedFrac[i], 4)))
                tableNumFeeds = QtWidgets.QTableWidgetItem(str(round(numFeeds[i], 4)))
                tablev_seed = QtWidgets.QTableWidgetItem(str(round(self.flowRates.get("v_seed")[i], 4)))
                tablev_surfactant = QtWidgets.QTableWidgetItem(str(round(self.flowRates.get("v_surfactant")[i], 4)))
                tablev_initiator = QtWidgets.QTableWidgetItem(str(round(self.flowRates.get("v_initiator")[i], 4)))
                tablev_monomer = QtWidgets.QTableWidgetItem(str(round(self.flowRates.get("v_monomer")[i], 4)))
                tablev_water = QtWidgets.QTableWidgetItem(str(round(self.flowRates.get("v_water")[i], 4)))
                self.parametersTable.setItem(0, i, tableTSC)
                self.parametersTable.setItem(1, i, tableRT)
                self.parametersTable.setItem(2, i, tableinitRatio)
                self.parametersTable.setItem(3, i, tablesurfactantRatio)
                self.parametersTable.setItem(4, i, tableSeedFrac)
                self.parametersTable.setItem(5, i, tableNumFeeds)
                self.parametersTable.setItem(6, i, tablev_seed)
                self.parametersTable.setItem(7, i, tablev_surfactant)
                self.parametersTable.setItem(8, i, tablev_initiator)
                self.parametersTable.setItem(9, i, tablev_monomer)
                self.parametersTable.setItem(10, i, tablev_water)

        if self.RAFTEmulsionCheckBox.isChecked():

            self.parametersTable.setRowCount(8)
            self.parametersTable.setColumnCount(numExp)
            self.parametersTable.setVerticalHeaderLabels([
                'TSC (g/g)',
                'Residence time (min)',
                'DP (mol/mol)',
                'CTA:I (mol/mol)', 
                'Macro-CTA flow rate (mL/min)',
                'Initiator flow rate (mL/min)',
                'Monomer flow rate (mL/min)',
                'Makeup water flow rate (mL/min)'
                ])
            
            for i in self.expIndex.astype(int):
                tableTSC = QtWidgets.QTableWidgetItem(str(round(w_f[i], 4)))
                tableRT = QtWidgets.QTableWidgetItem(str(round(tau[i], 4)))
                tableDP = QtWidgets.QTableWidgetItem(str(round(DP[i], 4)))
                tableCTAToInitiator = QtWidgets.QTableWidgetItem(str(round(ctaToInitiator[i], 4)))
                tablev_macroCTA = QtWidgets.QTableWidgetItem(str(round(self.flowRates.get("v_macroCTA")[i], 4)))
                tablev_initiator = QtWidgets.QTableWidgetItem(str(round(self.flowRates.get("v_initiator")[i], 4)))
                tablev_monomer = QtWidgets.QTableWidgetItem(str(round(self.flowRates.get("v_monomer")[i], 4)))
                tablev_water = QtWidgets.QTableWidgetItem(str(round(self.flowRates.get("v_water")[i], 4)))
                self.parametersTable.setItem(0, i, tableTSC)
                self.parametersTable.setItem(1, i, tableRT)
                self.parametersTable.setItem(2, i, tableDP)
                self.parametersTable.setItem(3, i, tableCTAToInitiator)
                self.parametersTable.setItem(4, i, tablev_macroCTA)
                self.parametersTable.setItem(5, i, tablev_initiator)
                self.parametersTable.setItem(6, i, tablev_monomer)
                self.parametersTable.setItem(7, i, tablev_water)

        self.parametersTable.resizeColumnsToContents()

        if not self.nGraph == 0: # Clear old graph if one has already been displayed
            plt.close(self.fig)
            self.methodGraphLayout.removeWidget(self.canvas)
            self.methodGraphLayout.removeWidget(self.toolbar)

        self.nGraph = 1 # Lets the function know that this isn't the first time displaying a graph to clear the old graph
        self.fig = plt.figure()
        self.fig.clear()
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(canvas=self.canvas)
        self.methodGraphLayout.addWidget(self.toolbar, 0, 3, QtCore.Qt.AlignTop)
        self.methodGraphLayout.addWidget(self.canvas, 1, 3, QtCore.Qt.AlignTop)

        # Create a dictionary of 'var(i) : parameterName' used for plotting the design in GUI #
        expDesign_vars = {}
        for i, name in zip(range(len(self.parameterNames)), self.parameterNames):
            var_name = "var%d" % i
            expDesign_vars[var_name] = name
        designColumns = list(expDesign.columns)

        if len(designColumns) == 1:
            self.ax = self.fig.add_subplot(111)
        if len(designColumns) == 2:
            self.ax = self.fig.add_subplot(111)
            self.fig.subplots_adjust(left=0.2, bottom=0.15)
        elif len(designColumns) == 3:
            self.ax = self.fig.add_subplot(111, projection='3d')
            self.fig.subplots_adjust(bottom=0.2)

        plotArgs = []
        for i in designColumns:
            newPlotArg = [expDesign[i]]
            plotArgs.append(newPlotArg)

        if len(designColumns) < 4:
            self.ax.scatter(*plotArgs, c='#0996D4', marker='o')
            if len(designColumns) == 1:
                x_variable = np.linspace(1, numExp+1, 1)
                y_variable = expDesign_vars.get("var0")
            x_variable = expDesign_vars.get("var0")
            y_variable = expDesign_vars.get("var1")
            self.ax.set_xlabel(x_variable)
            self.ax.set_ylabel(y_variable)
            if len(designColumns) == 3:
                z_variable = expDesign_vars.get("var2")
                self.ax.set_zlabel(z_variable)

    def startExperiment(self):
        self.startMsgBox = QMessageBox()
        self.startMsgBox.setWindowTitle("Confirm waste bottle emptied")
        self.startMsgBox.setText("Confirm you have emptied the waste bottle - do not continue the experiment without doing so")
        self.startMsgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        self.startMsgBox.setDefaultButton(QMessageBox.Cancel)
        self.startMsgBox.exec()
        if self.startMsgBox.standardButton(self.startMsgBox.clickedButton()) == QMessageBox.Ok:
            self.wasteVolumeAvailable = 850 # Resets the volume available in the waste bottle to 970 mL (3% margin of safety)
            self.wasteVolumeCounter.setText(str(round(self.wasteVolumeAvailable, 1)))
        if self.startMsgBox.clickedButton() is self.startMsgBox.button(QMessageBox.Cancel):
            return
        if self.startMsgBox.clickedButton() is self.startMsgBox.button(QMessageBox.Cancel):
            return
        self.startMsgBox2 = QMessageBox()
        self.startMsgBox2.setWindowTitle("Confirm seed syringe filled")
        self.startMsgBox2.setText("Confirm you have filled the seed syringe - do not continue the experiment without doing so")
        self.startMsgBox2.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        self.startMsgBox2.setDefaultButton(QMessageBox.Cancel)
        self.startMsgBox2.exec()
        if self.startMsgBox2.standardButton(self.startMsgBox2.clickedButton()) == QMessageBox.Ok:
            self.seedVolumeAvailable = 95 # Resets the volume available in the waste bottle to 980 mL (2% margin of safety)
            self.seedVolumeCounter.setText(str(round(self.seedVolumeAvailable, 2)))
        if self.startMsgBox2.clickedButton() is self.startMsgBox2.button(QMessageBox.Cancel):
            return
        self.saveSummaryData()
        self.stopThread = False
        self.manualStop = False
        self.runButton.setStyleSheet(
            "background-color: #8fd9a7;" 
            "color: black;" 
            "border-radius:5px;"
            "border: 1px solid;"
            "border-color: #002b0e;"
            "font-size: 10pt"
        )
        self.runButton.setText("\nExperiment in progress...\n")

        optimisationFlag=0
        if self.TSEMOCheckBox.isChecked():
            optimisationFlag=1
            strategy="TSEMO"

                #
                #
                # Need to start a TSEMO experiment but with an experimental design specified by the LHC earlier. These initial experiments should be stored in 'variablesDict' alongside other 'standard' variables
                #
                #

        if self.conventionalEmulsionCheckBox.isChecked():
            if self.kineticStudyCheckBox.isChecked():
                self.conventionalEPHandler.start(method="kinetics") # If kinetics study selected, pass method type as kinetic study to run the correct script
            else:
                self.conventionalEPHandler.start()

        if self.RAFTEmulsionCheckBox.isChecked():
            if self.kineticStudyCheckBox.isChecked():
                self.RAFTEPHandler.start(method="kinetics") # If kinetics study selected, pass method type as kinetic study to run the correct script
            else:
                self.RAFTEPHandler.start()
                
    def stopExperiment(self):
        self.stopThread = True
        self.manualStop = False
        self.main.methodHandler.runButton.setChecked(False)
        self.main.methodHandler.runButton.setText('\nStart experiment\n')
        self.main.methodHandler.runButton.setStyleSheet(
        "background-color: #0996D4;" 
        "color: white;" 
        "border-radius:5px;"
        "border: 1px solid;"
        "border-color: #2d5e86;"
        "font-size: 10pt"
        )

    def stopExperimentManual(self):
        # Stop an experiment with manual flag true (doesn't immediately start cleaning reactor)
        self.stopThread = True
        self.manualStop = True
        # self.experimentThread.join()
        print('Experiment aborted by user at ' + str(datetime.datetime.now()))
        self.main.methodHandler.runButton.setChecked(False)
        self.main.methodHandler.runButton.setText('\nStart experiment\n')
        self.main.methodHandler.runButton.setStyleSheet(
        "background-color: #0996D4;" 
        "color: white;" 
        "border-radius:5px;"
        "border: 1px solid;"
        "border-color: #2d5e86;"
        "font-size: 10pt"
        )

    def unpauseExperiment(self):
        self.pauseExperiment = False
        self.dialogueLabel.setHidden(True)
        self.ackButton.setHidden(True)

    def startSamplerCleaning(self):
        self.sampleCleanThread = threading.Thread(target=self.cleanSampler)
        self.sampleCleanThread.start()
        self.stopCleanThread = False

    def cleanSampler(self):
        # Script for conveniently cleaning all sample ports
        solventPump = self.main.controller.pump6
        solventValve = self.main.controller.valve1
        emulsionValve = self.main.controller.valve2
        outletValve = self.main.controller.valve5
        nSamples = int(self.samplesNumberText.text())

        while not self.stopCleanThread:
            emulsionValve.valveSwitch(5)
            solventPump.start()

            for i in range(2, nSamples + 2):
                outletValve.valveSwitch(i)
                portTimerStart = datetime.datetime.now()
                while not self.stopCleanThread and (datetime.datetime.now() - portTimerStart).seconds < 45:
                    time.sleep(1)
                if i == (nSamples - 1):
                    solventValve.valveSwitch(position='B')
            
            for i in reversed(range(2, nSamples + 2)):
                outletValve.valveSwitch(i)
                portTimerStart = datetime.datetime.now()
                while not self.stopCleanThread and (datetime.datetime.now() - portTimerStart).seconds < 45:
                    time.sleep(1)
            break
        print("Sampler cleaning finished")
        return
    
    def wasteEmptied(self):
        self.wasteMsgBox = QMessageBox()
        self.wasteMsgBox.setWindowTitle("Confirm waste bottle emptied")
        self.wasteMsgBox.setText("Confirm you have emptied the waste bottle - do not continue the experiment without doing so")
        self.wasteMsgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        self.wasteMsgBox.setDefaultButton(QMessageBox.Cancel)
        self.wasteMsgBox.exec()
        if self.wasteMsgBox.standardButton(self.wasteMsgBox.clickedButton()) == QMessageBox.Ok:
            print("Waste bottle volume reset")
            self.wasteVolumeAvailable = 980 # Resets the volume available in the waste bottle to 980 mL (2% margin of safety)
            return
        if self.wasteMsgBox.clickedButton() is self.wasteMsgBox.button(QMessageBox.Cancel):
            return

    def seedFilled(self):
        self.seedMsgBox = QMessageBox()
        self.seedMsgBox.setWindowTitle("Confirm seed syringe filled")
        self.seedMsgBox.setText("Confirm you have filled the seed syringe - do not continue the experiment without doing so")
        self.seedMsgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        self.seedMsgBox.setDefaultButton(QMessageBox.Cancel)
        self.seedMsgBox.exec()
        if self.seedMsgBox.standardButton(self.seedMsgBox.clickedButton()) == QMessageBox.Ok:
            print("Seed syringe volume reset")
            self.seedVolumeAvailable = 95 # Resets the volume available in the seed syringe to 95 mL (5% margin of safety)
            return
        if self.seedMsgBox.clickedButton() is self.seedMsgBox.button(QMessageBox.Cancel):
            return

    def searchDirectories(self):
        self.activeFolder = QtWidgets.QFileDialog.getExistingDirectory(self, caption='Select save folder')
        print(self.activeFolder)
        self.directoryLineEdit.setText(self.activeFolder)

    def searchTrainingData(self):
        self.trainingDataFolder = QtWidgets.QFileDialog.getExistingDirectory(self, caption='Select training data folder')
        self.trainingDataPathText.setText(self.trainingDataFolder)

    def setDataPath(self):
        self.saveFolder = (str(self.directoryLineEdit.text()) + '/')
        self.savePath = (self.saveFolder + self.experimentIDLineEdit.text())
        print(f"savePath: {self.savePath}")
        self.savePathCreated = 1

    def saveSummaryData(self):
        # Saves a summary of the experimental conditions in pre-programmed experiments (OFAAT and DoE) into the save path #
        self.expSummary = pd.DataFrame(self.variablesDict)
        if self.conventionalEmulsionCheckBox.isChecked(): # If conventional emulsion polymerisation selected - drop unrelated (default) variables from the summary 
            self.expSummary.drop(columns=["macroCTAStockConc",
                                            "macroCTAPolymerConc",
                                            "macroCTAMonomerConc",
                                            "densityMacroCTA",
                                            "molarMassMacroCTA",
                                            "molarMassMonomer",
                                            "molarMassInitiator",
                                            "DP",
                                            "ctaToInitiator"], inplace=True, axis=1)
        if self.RAFTEmulsionCheckBox.isChecked(): # If RAFT emulsion polymerisation selected - drop unrelated (default) variables from the summary
            self.expSummary.drop(columns=["Initiator_ratio",
                                            "w_s",
                                            "surfactantStockConc",
                                            "w_surfactantStock",
                                            "densitySurfactant",
                                            "densitySeed",
                                            "Surfactant_ratio"
                                            ], inplace=True, axis=1)
        summarySavePath = (self.savePath + '-summary.csv')
        pd.DataFrame(self.expSummary).to_csv(summarySavePath)
        print("Summary data saved")

    def getNewConditions(self, current_dataSet=None, strategy="TSEMO"):
        if strategy=="TSEMO":
            newConditions = optimiser.getNextExperiment(self.inputsDict, # inputs dictionary is defined by the method handler when parameters are calculated
                                                             self.objectivesDict, # objectives dictionary is defined by the method handler when parameters are calculated
                                                             self.optIterations, # number of total experimental iterations is held in the method handler
                                                             current_dataSet=current_dataSet) # current_dataSet can be updated externally and new conditions returned here
            print(f"new conditions = {newConditions}")
            for key in newConditions:
                newExpCount = len(newConditions[key]) # Count how many new experiments were added (usually 1 but gives room to request multiple suggested experiments)
                pass
            for key in self.variablesDict:
                if key in newConditions:
                    self.variablesDict[key].extend(newConditions[key]) # If the keyword (input variable) exists in the newConditions dict, append the new value(s) to the corresponding list in
                    self.optimisationSummary[key].extend(newConditions[key])
                elif not key in self.flowRates:
                    for i in range(newExpCount): # Else, count how many new conditions were added, and extend the constant variables (V_CSTR, nCSTRs etc.) with their first/default value
                        self.variablesDict[key].append(self.variablesDict[key][0])

            if self.conventionalEmulsionCheckBox.isChecked():
                newFlowRates = flowrateCalculator.calculateFlowrates('conventional',
                                                                        self.variablesDict,
                                                                        strategy="optimisation")                
                for key in newFlowRates:
                    self.variablesDict[key].extend(newFlowRates[key])
            return self.variablesDict
        
        
