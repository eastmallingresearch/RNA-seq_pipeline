# Function annotation for transcriptomes

## Interproscan 5

### Installation instructions
Interproscan can be a bit of a pain to get working...

These instructions were made using interproscan 5.27-66.0 

```shell
# download interproscan files
mkdir ~/build
cd build
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.27-66.0/interproscan-5.27-66.0-64-bit.tar.gz
tar -zxf interproscan-5.27-66.0-64-bit.tar.gz

# interproscan 5 requires gcc 4.8 or later (4.7 installed on cluster)
# See https://github.com/eastmallingresearch/seq_tools/blob/master/install_progs/instructions.md for instructions for getting gcc upgraded
# alternatively use gcc installed in Andy's home directory

### install prereq libraries: sfld ###
cd ~/build/interproscan-5.27-66.0/src/sfld/
sed -i -e 's/^EASEL_DIR=/EASEL_DIR=~\/build\/hmmer-3.1b2-linux-intel-x86_64\/easel/' Makefile # might need to specify full path
make 
mv sfld_preprocess ../../../bin/sfld/.
mv sfld_postprocess ../../../bin/sfld/. 

### install prereq libraries: cath-resolve-hits ###
cd ~/build/interproscan-5.27-66.0/bin/gene3d/4.1.0
wget "https://github.com/UCLOrengoGroup/cath-tools/releases/download/v0.16.2/cath-refine-align.centos6"
cd ../../../
echo -e "\n#cath-resolve settings\ncath.resolve.hits.switches.gene3d=--input-format=hmmer_domtblout" >>interproscan.properties

### install prereq libraries: rpsbproc ###
# rpsbproc is part of BLAST+. The below uses a prebuilt version which may or may not last for long
cd ~/build/interproscan-5.27-66.0/bin/blast/ncbi-blast-2.6.0+
wget "https://github.com/ebi-pf-team/interproscan/tree/master/core/jms-implementation/support-mini-x86-32/bin/blast/rpsbproc"
chmod rpsbproc 755

# You also need to create a rpsbproc.ini in the same location:
cat <<EOT>> ~/build/interproscan-5.27-66.0/bin/blast/ncbi-blast-2.6.0+/rpsbproc.ini
cdd.signature.library.release=3.16
cdd.signature.list.path=data/cdd/3.16/data/cddid.tbl
cdd.library.path=data/cdd/3.16/db/Cdd_NCBI
cdd.data.path=data/cdd/3.16/data
EOT

# the below is the "correct" way to get rpsbproc to work (blast+ install) - requires gcc 4.8+
mkdir cddblast
cd cddblast
wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz
tar xvzf ncbi-blast-2.6.0+-src.tar.gz
cd ncbi-blast-2.6.0+-src/c++/src/app/
wget -r --no-parent -l 1 -np -nd -nH -P rpsbproc ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/rpsbproc-src/
#edit Makefile.in and make sure SUB_PROJ is assigned two applications as follows: SUB_PROJ = blast rpsbproc
cd ../../
./compilers/unix/GCC.sh
make && make install 
#after compilation is complete (or update interprocscan.profile to installed file location)
cp ReleaseMT/bin/rpsblast ../../../bin/blast/ncbi-blast-2.6.0+/
cp ReleaseMT/bin/rpsbproc ../../../bin/blast/ncbi-blast-2.6.0+/


# Compiling interproscan (borrowing various other precompiled bits from Andy's home directory)
cd ~/build/interproscan-5.27-66.0
./configure \
 --prefix=/$HOME/usr/local  \
 --with-gmp="/home/armita/prog/gmp/gmp/6.1.2" \
 --with-mpfr="/home/armita/prog/mpfr/mpfr-3.1.5" \
 --with-mpc="/home/armita/prog/mpc/mpc-0.8.2" \
 --disable-multilib
```

#### Issues

Gene3D is still not working correctly
The translation option (-t n) is not working correctly (uses EMBOSS getorf which I couldn't get to work).
The included get_all_cds.pl (or get_longest_cds) and translate.pl can convert a nucleotide sequence(s) to protein:
```
get_all_cds.pl <mytranscriptome.fa|translate.pl|sed 's/\*//g' >myproteins.fa
```
It outputs all proteins at least 25 aa long - it's a bit slow mind. The sed gets rid of any leading or trailing stops (\*) as interproscan doesn't like them.   

### Running on cluster
Interproscan has a built in option to run on an SGE cluster - but I couldn't get it to work.  
Therefore I've written my own version which splits a transcriptome into multiple files and runs multiple simultaneous jobs (hardcoded to 10). Compared to the Interproscan inbuilt clustering this will take an additional core for each simultaneous job.  

To run interproscan on the cluster use interproscan.sh. This takes two compulsory + any iterproscan arguments.
```shell
interproscan.sh FOLDER FILE <options>
```
FOLDER is the location of the transcriptome file(s) - all output files will be copied here.  
FILE is the name of the transcriptome file (accepts shell expandable wild cards to specify multiple files).  
At least one of the options *must* specify which analyses to perform (-appl) - otherwise the Gene3D analysis will cause the JVM to crash. Below is the list of possible analyses:  
 -appl CDD,COILS,HAMAP,MobiDBLite,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM  
 
 Other interesting options:  
 -iprlookup (provides mappings from matched member database signatures to the InterPro entries that they are integrated into)  
 -goterms  (provides mappings to the Gene Ontology - based on InterPro (IPR) matches)  
 -pa (KEGG, MetaCyc and Reactome pathway lookup - based on InterPro (IPR) matches)  
 -dra (switch off residue level annotations - speeds up analysis)  


The cluster job will create a tempory folder <TEMPFOLDER> in <FOLDER> which will contain the transcriptome split into multiple file containing 1000 sequences, and a file called <FILE>_splitfiles.txt which contains the path to each of the multiple files.
 
Interproscan can take a while to run, if the job needs to be restarted, edit <FILE>_splitfiles.txt removing all files which have already been processed, then run restart_interproscan.sh.  

```shell
restart_interproscan.sh FOLDER FILE TEMPFOLDER <options>
```
