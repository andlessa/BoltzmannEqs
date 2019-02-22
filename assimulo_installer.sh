#!/bin/sh

homeDIR="$( pwd )"

#Get Sundials tarball
sundials="sundials-2.7.0.tar.gz"
URL=https://computation.llnl.gov/projects/sundials/download/$sundials
mkdir sundials_src;
echo "[installer] getting $sundials"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $sundials -C sundials_src --strip-components=1;
#Check if file was downloaded and extracted:
if [ ! "$(ls -A $homeDIR/sundials_src)" ]; then
  echo "SUNDIALS was not properly downloaded. Check the tarball path"
  exit
fi

#Make installation dir
mkdir sundials;
#Make build dir (sundials can not be built in the source or installation folders)
mkdir sundials_build;
cd sundials_build;
#Install SUNDIALS
echo "[installer] installing sundials";
cmake -DCMAKE_INSTALL_PREFIX=$homeDIR/sundials $homeDIR/sundials_src; make; make install;
cd ..;
#Check if installation was successful:
if [ ! "$(ls -A $homeDIR/sundials)" ]; then
  echo "SUNDIALS was not installed"
  exit
fi

rm -r sundials_src sundials_build $sundials


#Add sundials folder to path:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$homeDIR/sundials/lib

#Install Assimulo:
echo "[installer] installing Assimulo";
pip install --user assimulo

echo "[installer] Checking installation";
./check_installation.py

printf "\n\n\n Please, add sundials folder ($homeDIR/sundials) to the library path. Include the following line in your .bashrc:\n\n  LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$homeDIR/sundials/lib\n"

