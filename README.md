###LBNL Run14 D0 Spectra and Nuclear Modification Factor Analysis.  
  
Principal Authors:  
	Guannan Xie (guannanxie@lbl.gov)  
	Mustafa Mustafa (mmustafa@lbl.gov)  

- - -
###How to build this code:  
```bash
mkdir myAnalysis  
cd myAnalysis  

# Replace address below with your own fork if you have one  
git clone git@github.com:MustafaMustafa/auau200GeVRun14Ana.git

# Clone LBNL PicoHFLib  
git clone git@github.com:rnc-lbl/auau200GeVRun14.git

# Now you need to get StPicoDstMaker  
# if compiling at PDSF you need to get a klog token as below. You don't need this step at RCF.  
klog ­principal YOURRCFUSERNAME  
cvs co -r Run14_AuAu200_physics offline/users/dongx/pico/source/StPicoDstMaker  

mkdir StRoot  
cd StRoot  
ln -s ../auau200GeVRun14Ana/StPicoD0AnaMaker  
ln -s ../auau200GeVRun14Ana/StRoot/StPicoD0EventMaker  
ln -s ../auau200GeVRun14Ana/StRoot/StPicoPrescales  
ln -s ../auau200GeVRun14Ana/StRoot/StPicoHFMaker  
ln -s ../offline/users/dongx/pico/source/StPicoDstMaker  
starver SL15c  
cons  
```

###How to get a list of files:  
```bash
git clone git@github.com:rnc-lbl/fileLists.git  
# The list of daily D0 production will be under:  
# fileLists/Run14/AuAu/200GeV/physics/picoD0Lists/daily  
```

###How to run this code:  
```bash
cd myAnalysis
ln -s auau200GeVRun14Ana/StRoot/macros/runPicoD0AnaMaker.C  
root4star -l -b -q -x 'runPicoD0AnaMaker.C(“d0Trees.list”,”outputfile.root”)'  
```

###How to submit jobs:
```bash
# You cah find STAR Scheduler XML file under:  
cp -p auau200GeVRun14/starSubmit/submitPicoD0AnaMaker.xml  
# The README file contains a how to use.  
```
