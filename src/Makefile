CC = g++
CFLAGS = `root-config --cflags`  
LFLAGS = -O3 -Wl,--no-as-needed `root-config --libs`
CORE = CorrelationUtils.o  PIDUtils.o AnalysisFW.o AnalysisBinning.o SetupCustomTrackTree.o EvtSelection.o EvtAnalyzer.o SpectrumUtils.o
EXT = CorrelationUtils_FigureTools.o PhyConst.o GraphTools.o ContMatrix.o
FULL = $(CORE) $(EXT)
CMN =  commonSetup.h

all : ../bin/preprocess ../bin/process ../bin/preprocess_minimal

###################################################################################

### Proton, kaon and pion - charged hadron correlation
# process
../bin/process : process.o $(CMN) $(FULL)
	$(CC) $(LFLAGS) -o $@ $< $(FULL)

process.o : process.cpp 
	$(CC) $(CFLAGS) -o $@ -c $<

# process_light
../bin/process_light : process_light.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE)

process_light.o : process_light.cpp 
	$(CC) $(CFLAGS) -o $@ -c $<

# preprocess
../bin/preprocess : preprocess.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE)

preprocess.o : preprocess.cpp $(CMN) $(CORE)
	$(CC) $(CFLAGS) -o $@ -c $<

# preprocess_genlevel
../bin/preprocess_genlevel : preprocess_genlevel.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE)

preprocess_genlevel.o : preprocess_genlevel.cpp $(CMN) $(CORE)
	$(CC) $(CFLAGS) -o $@ -c $<

# SimplePreProc
../bin/SimplePreProc : SimplePreProc.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE)

SimplePreProc.o : SimplePreProc.cpp $(CMN) $(CORE)
	$(CC) $(CFLAGS) -o $@ -c $<


####################################################################################
### dEdxMap
../bin/makedEdxMap : makedEdxMap.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE) 

makedEdxMap.o : makedEdxMap.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


####################################################################################
### PIDSetup
../bin/PIDSetup : PIDSetup.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE) 

PIDSetup.o : PIDSetup.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

####################################################################################
### PIDStudy
../bin/PIDStudy : PIDStudy.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE) 

PIDStudy.o : PIDStudy.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

####################################################################################
### EtaPhiDistr
../bin/EtaPhiDistr : EtaPhiDistr.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE) 

EtaPhiDistr.o : EtaPhiDistr.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

####################################################################################
### MC_Contamination_Matrix
../bin/MC_Contamination_Matrix : MC_Contamination_Matrix.o $(CMN) $(FULL)
	$(CC) $(LFLAGS) -o $@ $< $(FULL) 

MC_Contamination_Matrix.o : MC_Contamination_Matrix.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

# MC_Contamination_Matrix_viewer
../bin/MC_Contamination_Matrix_viewer : MC_Contamination_Matrix_viewer.o $(CMN) $(FULL)
	$(CC) $(LFLAGS) -o $@ $< $(FULL) 

MC_Contamination_Matrix_viewer.o : MC_Contamination_Matrix_viewer.cpp
	$(CC) $(CFLAGS)-o $@ -c $<  

####################################################################################
### TrackCorrection
../bin/TrackCorrection : TrackCorrection.o $(CMN) PIDUtils.o AnalysisFW.o AnalysisBinning.o PhyConst.o SetupCustomTrackTree.o EvtSelection.o EvtAnalyzer.o
	$(CC) $(LFLAGS) -o $@ $< PIDUtils.o AnalysisBinning.o PhyConst.o AnalysisFW.o SetupCustomTrackTree.o EvtSelection.o EvtAnalyzer.o

TrackCorrection.o : TrackCorrection.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

###########################################
###  *** Framework and utility test *** ###
# AnalysisFWTest
../bin/AnaFwTest : AnaFwTest.o $(CMN) $(FULL)
	$(CC) $(LFLAGS) -o $@ $< $(FULL) 

AnaFwTest.o : AnaFwTest.cpp
	$(CC) $(CFLAGS)-o $@ -c $<  

# TrackCorrection_viewer
../bin/TrackCorrection_viewer : TrackCorrection_viewer.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE) 

TrackCorrection_viewer.o : TrackCorrection_viewer.cpp
	$(CC) $(CFLAGS)-o $@ -c $<  

# UtilityTest
../bin/utilitytest : UtilityTest.o $(CMN) $(CORE)
	$(CC) $(LFLAGS) -o $@ $< $(CORE) 

UtilityTest.o : UtilityTest.cpp
	$(CC) $(CFLAGS)-o $@ -c $<  

# MC_CrosscheckFig
../bin/MC_CrosscheckFig : MC_CrosscheckFig.o $(CMN) AnalysisBinning.o
	$(CC) $(LFLAGS) -o $@ $<  AnalysisBinning.o

MC_CrosscheckFig.o : MC_CrosscheckFig.cpp
	$(CC) $(CFLAGS)-o $@ -c $<  AnalysisBinning.o

# Compare_vns
../bin/Compare_vns : Compare_vns.o $(CMN) AnalysisBinning.o
	$(CC) $(LFLAGS) -o $@ $< AnalysisBinning.o

Compare_vns.o : Compare_vns.cpp
	$(CC) $(CFLAGS)-o $@ -c $< AnalysisBinning.o

# Compare_vns_multi
../bin/Compare_vns_multi : Compare_vns_multi.o $(CMN) AnalysisBinning.o
	$(CC) $(LFLAGS) -o $@ $< AnalysisBinning.o

Compare_vns_multi.o : Compare_vns_multi.cpp
	$(CC) $(CFLAGS)-o $@ -c $< AnalysisBinning.o

###########################################
### *** Framework classes
# AnalysisBinning
AnalysisBinning.o : AnalysisBinning.cpp AnalysisBinning.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# AnalysisFW
AnalysisFW.o : AnalysisFW.cpp AnalysisFW.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# GraphTools
GraphTools.o : GraphTools.cpp GraphTools.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# PhyConst
PhyConst.o : PhyConst.cpp PhyConst.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# CorrelationUtils
CorrelationUtils.o : CorrelationUtils.cpp CorrelationUtils.h
	$(CC) $(CFLAGS) -o $@ -c $<

# CorrelationUtils_FigureTools
CorrelationUtils_FigureTools.o : CorrelationUtils_FigureTools.cpp CorrelationUtils.h
	$(CC) $(CFLAGS) -o $@ -c $<

# ContMatrix
ContMatrix.o : ContMatrix.cpp ContMatrix.h
	$(CC) $(CFLAGS) -o $@ -c $<

# PIDUtils
PIDUtils.o : PIDUtils.cpp PIDUtils.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# SpectrumUtils
SpectrumUtils.o : SpectrumUtils.cpp SpectrumUtils.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# SetupCustomTrackTree
SetupCustomTrackTree.o : SetupCustomTrackTree.cpp SetupCustomTrackTree.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# SetupSkimTree
SetupSkimTree.o : SetupSkimTree.cpp SetupSkimTree.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# SetupEvtTree
SetupEvtTree.o : SetupEvtTree.cpp SetupEvtTree.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# EvtSelection
EvtSelection.o : EvtSelection.cpp EvtSelection.h
	$(CC) $(CFLAGS) -c $< -o $@ 

# EvtAnalyzer
EvtAnalyzer.o : EvtAnalyzer.cpp EvtAnalyzer.h
	$(CC) $(CFLAGS) -c $< -o $@ 
