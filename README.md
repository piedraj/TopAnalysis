Everything starts here
====

The most common machines to run are lxplus.

    ssh -Y lxplus.cern.ch -o ServerAliveInterval=240

Once logged in we need to set the architecture.

    setenv SCRAM_ARCH slc6_amd64_gcc491

If not running from lxplus, the CMS environment has to be set.

    source /cvmfs/cms.cern.ch/cmsset_default.csh

Now we can choose our favorite CMSSW release.

    cmsrel CMSSW_7_4_0_pre9
    cd CMSSW_7_4_0_pre9/src
    cmsenv


It is time to get the material
====

Go to the master repository (https://github.com/piedraj/TopAnalysis) and click **Fork** in the top-right corner of the page. Now get the code in your working area.

    git clone https://github.com/YOUR-GIT-USER-NAME/TopAnalysis.git TopAnalysis
    mkdir TopAnalysis/TopTreeProducer/interface


Do a test run
====

    cd CMSSW_7_4_0_pre9/src
    cmsenv
    scram b -j 10
    voms-proxy-init -voms cms
    cd TopTreeProducer/test
    cmsRun skimToTreeSUSYMCtfs.py


CRAB3
====

    cd CMSSW_7_4_0_pre9/src
    cmsenv
    source /cvmfs/cms.cern.ch/crab3/crab.csh
    voms-proxy-init -voms cms
    cd TopTreeProducer/test/
    crab submit -c crabConfig.py
    crab status --dir crab_DYJetsToLL_PU20bx25_PHYS14


It is commit time
====

First get the latest changes in the repository, if any.

    git pull https://github.com/piedraj/TopAnalysis.git

And then commit your changes.

    git status
    git add <filepattern>
    git commit -m 'Modified'
    git push origin master

If the changes have been made in a fork of the master, go to https://github.com/YOUR-USERNAME/TopAnalysis.git and click **Pull Request**.

