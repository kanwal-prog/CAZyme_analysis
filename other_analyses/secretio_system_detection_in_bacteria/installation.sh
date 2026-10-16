#Bacterial Secretion Systems Detection using TXSScan / MacSyFinder

# create conda envrionemnt
conda create -n txsscan python=3.11 -y
# activate the environment
conda activate txsscan
# install the macsyfinder package latest version 2.1.6
conda install -c bioconda macsyfinder -y


#check available msf_data
msf_data available 
msf_data install --user TXSScan
# to show
msf_data show TXSScan	


# for TFFscan
msf_data install --user TFFscan
# for CasFinder
msf_data install --user CasFinder
# for CONJScan
msf_data install --user CONJScan

#TXSScan (1.1.4)            - TXSScan - Models for 14 types of bacterial protein secretion systems and appendages (T1SS, T2SS, T3SS, T4aP, T4bP, pT4SS, T5SSs, T6SSs, Flagellum, Tad pilus, T9SS, Com pilus, MSH, Archaeal TFF).
#TFFscan (1.0.1)            - TFF-SF - Models for the type IV filament (TFF) superfamily (including T2SS, T4P, Tad pili and Competence machinery).
#CasFinder (3.1.1)          - CasFinder - Models for detection of CRISPR-Cas systems.
#CONJScan (2.1.0)           - CONJScan - Models for detection of conjugative and mobilisable elements.

TXSScan
    ├-archaea
    │   └-Archaeal-T4P
    └-bacteria
         ├-diderm
         │   ├-Flagellum
         │   ├-MSH
         │   ├-T1SS
         │   ├-T2SS
         │   ├-T3SS
         │   ├-T4aP
         │   ├-T4bP
         │   ├-T5aSS
         │   ├-T5bSS
         │   ├-T5cSS
         │   ├-T6SSi
         │   ├-T6SSii
         │   ├-T6SSiii
         │   ├-T9SS
         │   ├-Tad
         │   ├-pT4SSi
         │   └-pT4SSt
         └-monoderm
              └-ComM

TXSScan (1.1.4) : 19 models