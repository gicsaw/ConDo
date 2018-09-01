# ConDo
Contact based protein Domain boundary prediction method

# Pre-requisite:

PSIBLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/ 
! (do not use blast+) with sqeucne database such as UniRef or NR 

HHblitz: https://github.com/soedinglab/hh-suite.git

PSIPRED: http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/

SANN: https://github.com/newtonjoo/sann  

or https://lee.kias.re.kr

Jackhmmer:http://hmmer.org/  
!must install hmmer/easel   (enclosed in hmmer)

UniRef90: https://www.uniprot.org/downloads  
!Use UniRef90 (recommended)

CCMPRED:git clone --recursive https://github.com/soedinglab/CCMpred.git

python2

numpy

KERAS with TensorFlow or theano 

gcc or icc

#Installation:

git clone https://github.com/gicsaw/ConDo

cd ConDo

gcc src/feature.c -o bin/feature -lm -fopenmp -O2

or icc src/feature.c -o bin/feature -qopenmp -O2

Edit ConDodir variable in bin/ConDo.sh 

Edit hhpath, condodir, database, and jackhmmerbin variables in bin/run_jackhmmer.sh 

Edit blastbin, dbname, psipred, condodir, sann, NNDB_HOME variables in bin/gen_features.sh 

Edit ccmpredbindir variables in bin/run_ccmpred.sh

#Run a example:
We prepared two targets such as 1c7cA and 1sxjH in examples dir

cd examples/$target   !replace $target to 1c7cA or 1sxjH

Condo.sh $target.fasta $ncpu 


Output file is $target.ConDo

First and second columns of the output file are residue index and domain boundary score, respectively

The cut-off of score is 1.4 

In gnuplot, plot "$target.ComDo" u 1:2 w lp, 1.4


#etc: 
bin/feature # generate input features of Machine learning and some other output files such as PAS, contact mat, modularity of contact. 

Input files are: 

$target.fasta    ! target sequence

$target.ss2      ! Secondary Structure predicted by PSIPRED

$target.a22      ! Solvent Accessibility predicted by SANN

$target.a3       ! Solvent Accessibility predicted by SANN

$target.ck2      ! sequence profile converted from chk of blast 

$target.msa      ! multiple sequence alignment converted by Jackhammer 

$target.ccmpred  ! predicted contact by CCMPRED

Outout files are
$target_feature.txt  ! input features for machine

$target_PAS3.txt     ! PAS information 

result_ccm2.txt      ! predicted contact after filtering 

community_ccm2.txt   ! Modularity of predicted contact

How to show 

In gnuplot,

set size square

plot $target_PAS3.txt u 1:2:3 w image

plot result_ccm2.txt u 1:2:3 w image

plot community_ccm2.txt u 1:2 w lp

