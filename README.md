Everything starts here
====

    ssh -Y lxplus.cern.ch -o ServerAliveInterval=240

    setenv SCRAM_ARCH slc6_amd64_gcc481
    cmsrel CMSSW_7_2_0
    cd CMSSW_7_2_0/src
    cmsenv


Get the material and compile it
====

    git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720
    git clone https://github.com/piedraj/TopAnalysis.git TopAnalysis
    mkdir TopAnalysis/TopTreeProducer/interface

    scram b -j 10


Do a test run
====

    cmsenv
    voms-proxy-init

    cd TopTreeProducer/test
    cmsRun skimToTreeSUSYMCtfs.py


CRAB3
====

    cmsenv
    source /cvmfs/cms.cern.ch/crab3/crab.csh
    voms-proxy-init

    cd TopTreeProducer/test/
    crab submit -c crabConfig.py
    crab status --dir crab_TTJets_PU30bx50


It is commit time
====

First get the latest changes in the repository, if any.

    git pull

And then commit your changes.

    git status
    git add <filepattern>
    git commit -m 'Modified'
    git push origin master

