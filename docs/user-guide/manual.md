## Warning: All the following steps are in the Linux system

## Download and unzip the main package

Download the CancerVar zip package at [here](https://github.com/WGLab/CancerVar/archive/master.zip) using `wget https://github.com/WGLab/CancerVar/archive/master.zip -O CancerVar.zip`:

```
qli@login1|:~> wget https://github.com/WGLab/CancerVar/archive/master.zip -O CancerVar.zip
--2022-01-24 20:50:53--  https://github.com/WGLab/CancerVar/archive/master.zip
Resolving github.com (github.com)... 140.82.114.3
Connecting to github.com (github.com)|140.82.114.3|:443... connected.
HTTP request sent, awaiting response... 302 Found
Location: https://codeload.github.com/WGLab/CancerVar/zip/master [following]
--2022-01-24 20:50:53--  https://codeload.github.com/WGLab/CancerVar/zip/master
Resolving codeload.github.com (codeload.github.com)... 140.82.114.10
Connecting to codeload.github.com (codeload.github.com)|140.82.114.10|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: unspecified [application/zip]
Saving to: 'CancerVar.zip'

    [      <=>                                                                                     ] 5,448,305   2.78MB/s   in 1.9s

2022-01-24 20:50:56 (2.78 MB/s) - 'CancerVar.zip' saved [5448305]

```

 

Assume that we have successfully  downloaded the InterVar package as `CancerVar.zip` and used `unzip CancerVar.zip` to unpack the package.

```
qli@login1|:~> unzip CancerVar.zip
Archive:  CancerVar.zip
b34bc913eba64c9b72817841f871d2ed5797e077
   creating: CancerVar-master/
  inflating: CancerVar-master/.gitignore
  inflating: CancerVar-master/CancerVar.py
   creating: CancerVar-master/OPAI/
   creating: CancerVar-master/OPAI/saves/
  inflating: CancerVar-master/OPAI/saves/ensemble.pt
  inflating: CancerVar-master/OPAI/saves/evs.pt
  inflating: CancerVar-master/OPAI/saves/nonmissing_db.npy
   creating: CancerVar-master/OPAI/scripts/
  inflating: CancerVar-master/OPAI/scripts/feature_preprocess.py
  inflating: CancerVar-master/OPAI/scripts/myDis.py
  inflating: CancerVar-master/OPAI/scripts/opai_predictor.py
  inflating: CancerVar-master/README.md
   creating: CancerVar-master/cancervardb/
  inflating: CancerVar-master/cancervardb/LOF.genes.exac_me_cancers
  inflating: CancerVar-master/cancervardb/cancer_census.genes
  inflating: CancerVar-master/cancervardb/cancers_genes.list_kegg.txt
  inflating: CancerVar-master/cancervardb/cancervar.cancer.types
  inflating: CancerVar-master/cancervardb/cancervar.out.txt
  inflating: CancerVar-master/cancervardb/knownGeneCanonical.txt
  inflating: CancerVar-master/cancervardb/mim2gene.txt
  inflating: CancerVar-master/cancervardb/mim_orpha.txt
  inflating: CancerVar-master/cancervardb/mim_pheno.txt
  inflating: CancerVar-master/cancervardb/orpha.txt
  inflating: CancerVar-master/config.ini
   creating: CancerVar-master/docs/
  inflating: CancerVar-master/docs/favicon.ico
   creating: CancerVar-master/docs/img/
 extracting: CancerVar-master/docs/img/new.png
   creating: CancerVar-master/docs/misc/
  inflating: CancerVar-master/docs/misc/contributing.md
  inflating: CancerVar-master/docs/misc/credit.md
  inflating: CancerVar-master/docs/misc/faq.md
  inflating: CancerVar-master/docs/misc/whatsnew.md
  inflating: CancerVar-master/docs/readme.md
   creating: CancerVar-master/docs/user-guide/
  inflating: CancerVar-master/docs/user-guide/download.md
  inflating: CancerVar-master/docs/user-guide/startup.md
   creating: CancerVar-master/example/
  inflating: CancerVar-master/example/FDA.hg19_multianno.txt
  inflating: CancerVar-master/example/FDA.hg19_multianno.txt.cancervar
  inflating: CancerVar-master/example/FDA.hg19_multianno.txt.cancervar.ensemble.csv
  inflating: CancerVar-master/example/FDA.hg19_multianno.txt.cancervar.ensemble.pred
  inflating: CancerVar-master/example/FDA.hg19_multianno.txt.cancervar.evs.csv
  inflating: CancerVar-master/example/FDA.hg19_multianno.txt.cancervar.evs.pred
  inflating: CancerVar-master/example/FDA.hg19_multianno.txt.grl_p
  inflating: CancerVar-master/example/FDA_hg19.av

```
Go to the folder `cd CancerVar-master` and check the files using `ls -alrt .`

```
qli@login1|:~> cd CancerVar-master
qli@login1|:~/CancerVar-master> ls -alrt .
total 198
-rw-r-----  1 qli qli 11833 Jan 24 01:23 README.md
drwxr-x---  4 qli qli  4096 Jan 24 01:23 OPAI
-rw-r-----  1 qli qli  1045 Jan 24 01:23 .gitignore
drwxr-x---  2 qli qli  4096 Jan 24 01:23 example
drwxr-x---  5 qli qli  4096 Jan 24 01:23 docs
-rw-r-----  1 qli qli  2588 Jan 24 01:23 config.ini
-rwxr-xr-x  1 qli qli 87114 Jan 24 01:23 CancerVar.py
drwxr-x---  2 qli qli  4096 Jan 24 01:23 cancervardb
drwxr-x---  6 qli qli  4096 Jan 24 01:23 .
drwx------ 49 qli qli 32768 Jan 24 20:52 ..

```
Now you can find the main python program as `CancerVar.py`, and test the main program of CancerVar  can run properly or not by `python CancerVar.py`;

Also the main OPAI scripts in folder "OPAI", the OPAI need to install 4 python modules, These are **numpy** https://numpy.org, **pandas** https://pandas.pydata.org , **scikit-learn** https://scikit-learn.org and **pytorch** https://pytorch.org.

There are two ways to install these modules:

- Using CONDA and manage the environment.
```
     conda create  -n opai python=3.6
     conda activate opai
     conda install -c anaconda numpy pandas scikit-learn
     conda install -c pytorch pytorch=1.9
```

- Using pip
```
    python3.6 -m pip install numpy --user
    python3.6 -m pip install pandas --user
    python3.6 -m pip install scikit-learn --user
    python3.6 -m pip install pytorch --user
```
after install the modules, we can test CancerVar and OPAI :


```
@login1|:~/CancerVar-master> python CancerVar.py
Usage: CancerVar.py [OPTION] -i  INPUT -o  OUTPUT ...
       CancerVar.py  --config=config.ini ...


=============================================================================
CancerVar Interpretation of Pathogenic/Benign for cancer variants (v.1.1)
=============================================================================
........................................................................
..%%%%....%%%%...%%..%%...%%%%...%%%%%%..%%%%%...%%..%%...%%%%...%%%%%..
.%%..%%..%%..%%..%%%.%%..%%..%%..%%......%%..%%..%%..%%..%%..%%..%%..%%.
.%%......%%%%%%..%%.%%%..%%......%%%%....%%%%%...%%..%%..%%%%%%..%%%%%..
.%%..%%..%%..%%..%%..%%..%%..%%..%%......%%..%%...%%%%...%%..%%..%%..%%.
..%%%%...%%..%%..%%..%%...%%%%...%%%%%%..%%..%%....%%....%%..%%..%%..%%.
........................................................................
New datasets downloading:   https://cancervar.wglab.org/databases/

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -c config.ini, --config=config.ini
                        The config file of all options. it is for your own
                        configure file.You can edit all the options in the
                        configure and if you use this options,you can ignore
                        all the other options bellow
  -b hg19, --buildver=hg19
                        The genomic build version, it can be hg19 and  hg38
  -i example/ex1.avinput, --input=example/ex1.avinput
                        The input file contains your variants
  --input_type=AVinput  The input file type, it can be  AVinput(Annovar's
                        format),VCF(VCF with single sample),VCF_m(VCF with
                        multiple samples)
  --cancer_type=CANCER  The cancer type, please check the doc for the details
                        of cancer types: Adrenal_Gland Bile_Duct Bladder Blood
                        Bone Bone_Marrow Brain Breast Cancer_all Cervix
                        Colorectal Esophagus Eye Head_and_Neck Inflammatory
                        Intrahepatic Kidney Liver Lung Lymph_Nodes
                        Nervous_System Other Ovary Pancreas Pleura Prostate
                        Skin Soft_Tissue Stomach Testis Thymus Thyroid
                        Uterus),if you are using avinput file, you can can
                        specify the cancer type in the 6th column
  -o example/myanno, --output=example/myanno
                        The prefix of output file which contains the results,
                        the file of results will be as [$$prefix].cancervar

  CancerVar Other Options:
    -t cancervardb, --database_cancervar=cancervardb
                        The  database location/dir for the CancerVar dataset
                        files
    -s your_evidence_file, --evidence_file=your_evidence_file
                        User specified Evidence file for each variant

     How to add your own Evidence for each Variant:
     Prepare your own evidence  file as tab-delimited,the line format:
    The format for upgrad/downgrade of criteria should be like:
    grade_CBPx=2;          3 for Strong; 2 for Moderate; 1 for Supporting)
    Chr Pos_start Pos_end Ref Alt CBP2=1;grade_CBP2=2;CBP9=1

  Annovar Options:
    Caution: check these options from manual of Annovar. The ANNOVAR
    version should be >=  2016-02-01, older verions of ANNOVAR will bring
    problems.

    --table_annovar=./table_annovar.pl
                        The Annovar perl script of table_annovar.pl
    --convert2annovar=./convert2annovar.pl
                        The Annovar perl script of convert2annovar.pl
    --annotate_variation=./annotate_variation.pl
                        The Annovar perl script of annotate_variation.pl
    -d humandb, --database_locat=humandb
                        The  database location/dir for the annotation datasets
    --skip_annovar      Skip the Annovar annotation, this can be true only
                        after you  already got the annovar annotation
                        results

  Examples:
    ./CancerVar.py -c config.ini  # Run the examples in config.ini
    ./CancerVar.py -i your_input  --input_type=VCF  -o your_output

```
test OPAI

```
qli@login1|:~/CancerVar-master> python3.6  OPAI/scripts/opai_predictor.py -h
usage: opai_predictor.py [-h] -i INPUT -v CANCERVAR_PATH [-m METHOD]
                         [-d DEVICE] -c CONFIG -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        the path to input feature
  -v CANCERVAR_PATH, --cancervar_path CANCERVAR_PATH
                        the path to cancervar file
  -m METHOD, --method METHOD
                        use evs features or ensemble features (option: evs,
                        ensemble)
  -d DEVICE, --device DEVICE
                        device used for dl-based predicting (option: cpu,
                        cuda)
  -c CONFIG, --config CONFIG
                        the path to trained model file
  -o OUTPUT, --output OUTPUT
                        the path to output


```


The `CancerVar.py` and OPAI can run when python version >= 3.6 .
If you see above screen output after `python3.6 CancerVar.py` and `python3.6  OPAI/scripts/opai_predictor.py -h` , that mean the CancerVar and OPAI can run on your system without problem.
Otherwise, please check your python version using `env python --version`

Next we need to download the ANNOVAR dataset.

## Download Third-party Program and datasets

Several third-party researchers have provided additional annotation program and datasets that can be used by CancerVar directly. However, users need to agree to specific license terms set forth by the third parties:


* ANNOVAR main package : The latest version of ANNOVAR (2016Feb01) can be downloaded [here](http://www.openbioinformatics.org/annovar/annovar_download_form.php) (registration required).  ANNOVAR is written in Perl and can be run as a standalone application on diverse hardware systems where standard Perl modules are installed.


Now,assume that we have downloaded ANNOVAR package and used `tar xvfz annovar.latest.tar.gz` to unpack the package. You will see that the `bin/` directory contains several Perl programs with .pl suffix. 

Then please copy  ANNOVAR perl files: `annotate_variation.pl` `table_annovar.pl` `convert2annovar.pl` `coding_change.pl` to CancerVar's folder of `CancerVar-master`.

```
qli@sched1|:~/tools/annovar> cp -f annotate_variation.pl table_annovar.pl convert2annovar.pl ~/CancerVar-master

```

Please go to CancerVar's install folder of `CancerVar-master`,check and edit the config.ini using `vim config.ini`, please check these lines in the config.ini . The names and locations in config.ini  should match with your downloaded files:
 

`convert2annovar = ./convert2annovar.pl`
This line is for location of ANNOVAR's `convert2annovar.pl` file

`table_annovar = ./table_annovar.pl`
This line is for location of ANNOVAR's `table_annovar.pl` file

`annotate_variation = ./annotate_variation.pl`
This line is for location of ANNOVAR's `table_annovar.pl` file



## Prepare your input files:
The format of input file can be VCF or AVinput as ANNOVAR input,actually the input file only need information of Chr,Position start,Position end, Reference_allele,Alternative_allele.
```
qli@login1|:~/CancerVar-master> head -4 example/FDA_hg19.av
9 133748283 133748283 C T
14 105246551 105246551 C T
7 140453136 140453136 A T
7 55249071 55249071 C T

```

The input file type can be specified by option of  `--input_type`, there are three types:  AVinput(Annovar's format),VCF(VCF with single sample),VCF_m(VCF with multiple samples)

## Run the Annotation on the  file `example/FDA_hg19.av`.
Please be advice that for the first running, the CancerVar will use the `perl ./annotate_variation.pl` to download the necessary ANNOVAR datasests,it will take some time.(If you also were ANNOVAR user before, you can specify the ANNOVAR's database_location by option `--database_locat`, you also can edit the config.ini, find the line `database_locat = humandb` and replace with your location  , InterVar will check if all database file exist in you provided location). From second running , CancerVar will not  download the same ANNOVAR's datasets again.

```
@cm-a1-n005|:~/CancerVar-master> python CancerVar.py -i example/FDA_hg19.av -o example/FDA
=============================================================================
CancerVar
Interpretation of Pathogenic/Benign for cancer variants (v.1.1)
=============================================================================
........................................................................
..%%%%....%%%%...%%..%%...%%%%...%%%%%%..%%%%%...%%..%%...%%%%...%%%%%..
.%%..%%..%%..%%..%%%.%%..%%..%%..%%......%%..%%..%%..%%..%%..%%..%%..%%.
.%%......%%%%%%..%%.%%%..%%......%%%%....%%%%%...%%..%%..%%%%%%..%%%%%..
.%%..%%..%%..%%..%%..%%..%%..%%..%%......%%..%%...%%%%...%%..%%..%%..%%.
..%%%%...%%..%%..%%..%%...%%%%...%%%%%%..%%..%%....%%....%%..%%..%%..%%.
........................................................................
      New datasets downloading:
  https://cancervar.wglab.org/databases/

Notice: Your command of CancerVar is ['CancerVar.py', '-i', 'example/FDA_hg19.av', '-o', 'example/FDA']
Warning: Your specified evidence file [ None ], the analysis will take your additional evidence.
INFO: The options are {'table_annovar': './table_annovar.pl', 'exclude_snps': 'cancervardb/ext.variants.hg19', 'annotate_variation': './annotate_variation.pl', 'current_version': 'CancerVar_20200119', 'evidence_file': 'None', 'public_dev': 'https://github.com/WGLab/CancerVar/releases', 'otherinfo': 'TRUE', 'database_names': 'refGene esp6500siv2_all 1000g2015aug avsnp147 dbnsfp30a clinvar_20190305 exac03 dbscsnv11 dbnsfp31a_interpro ensGene knownGene cosmic70 icgc21 gnomad_genome', 'mim_pheno': 'cancervardb/mim_pheno.txt', 'cancer_pathway': 'cancervardb/cancers_genes.list_kegg.txt', 'cancers_types': 'cancervardb/cancervar.cancer.types', 'buildver': 'hg19', 'onetranscript': 'FALSE', 'mim2gene': 'cancervardb/mim2gene.txt', 'orpha': 'cancervardb/orpha.txt', 'inputfile_type': 'AVinput', 'knowngenecanonical': 'cancervardb/knownGeneCanonical.txt', 'cancers_genes': 'cancervardb/cancer_census.genes', 'convert2annovar': './convert2annovar.pl', 'database_locat': 'humandb', 'database_cancervar': 'cancervardb', 'lof_genes': 'cancervardb/LOF.genes.exac_me_cancers', 'cancervar_markers': 'cancervardb/cancervar.out.txt', 'outfile': 'example/FDA', 'disorder_cutoff': '0.01', 'mim_orpha': 'cancervardb/mim_orpha.txt', 'inputfile': 'example/FDA_hg19.av'}
Warning: the folder of humandb is already created!
perl ./table_annovar.pl example/FDA_hg19.av humandb -buildver hg19 -remove -out example/FDA -protocol refGene,ensGene,knownGene,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,dbnsfp30a,dbscsnv11,dbnsfp31a_interpro,clinvar_20190305,cosmic91,icgc28,gnomad_genome  -operation  g,g,g,f,f,f,f,f,f,f,f,f,f,f  -nastring . --otherinfo
NOTICE: the --polish argument is set ON automatically (use --nopolish to change this behavior)
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=refGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile example/FDA.refGene -exonsort -nofirstcodondel example/FDA_hg19.av humandb>
NOTICE: Output files are written to example/FDA.refGene.variant_function, example/FDA.refGene.exonic_variant_function
NOTICE: Reading gene annotation from humandb/hg19_refGene.txt ... Done with 52068 transcripts (including 11837 without coding sequence annotation) for 26464 unique genes
NOTICE: Processing next batch with 22 unique variants in 22 input lines
NOTICE: Reading FASTA sequences from humandb/hg19_refGeneMrna.fa ... Done with 21 sequences
WARNING: A total of 356 sequences will be ignored due to lack of correct ORF annotation

NOTICE: Running with system command <coding_change.pl  example/FDA.refGene.exonic_variant_function.orig humandb/hg19_refGene.txt humandb/hg19_refGeneMrna.fa -alltranscript -out example/FDA.refGene.fa -newevf example/FDA.refGene.exonic_variant_function>
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=ensGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype ensGene -outfile example/FDA.ensGene -exonsort -nofirstcodondel example/FDA_hg19.av humandb>
NOTICE: Output files are written to example/FDA.ensGene.variant_function, example/FDA.ensGene.exonic_variant_function
NOTICE: Reading gene annotation from humandb/hg19_ensGene.txt ... Done with 196501 transcripts (including 101155 without coding sequence annotation) for 57905 unique genes
NOTICE: Processing next batch with 22 unique variants in 22 input lines
NOTICE: Reading FASTA sequences from humandb/hg19_ensGeneMrna.fa ... Done with 40 sequences
WARNING: A total of 6780 sequences will be ignored due to lack of correct ORF annotation

NOTICE: Running with system command <coding_change.pl  example/FDA.ensGene.exonic_variant_function.orig humandb/hg19_ensGene.txt humandb/hg19_ensGeneMrna.fa -alltranscript -out example/FDA.ensGene.fa -newevf example/FDA.ensGene.exonic_variant_function>
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=knownGene

......
......
......

NOTICE: Running system command <annotate_variation.pl -filter -dbtype gnomad_genome -buildver hg19 -outfile example/FDA example/FDA_hg19.av humandb -otherinfo>
NOTICE: Output file with variants matching filtering criteria is written to example/FDA.hg19_gnomad_genome_dropped, and output file with other variants is written to example/FDA.hg19_gnomad_genome_filtered
NOTICE: Processing next batch with 22 unique variants in 22 input lines
NOTICE: Database index loaded. Total number of bins is 28127612 and the number of bins to be scanned is 15
NOTICE: Scanning filter database humandb/hg19_gnomad_genome.txt...Done
-----------------------------------------------------------------
NOTICE: Multianno output file is written to example/FDA.hg19_multianno.txt
annovar_outfile is example/FDA.hg19_multianno.txt
Notice: Begin the variants interpretation by CancerVar
Notice: About 22 lines in your variant file!
Notice: About 22 variants has been processed by CancerVar
Notice: The CancerVar is finished, the output file is [ example/FDA.hg19_multianno.txt.cancervar ]
=============================================================================
........................................................................
..%%%%....%%%%...%%..%%...%%%%...%%%%%%..%%%%%...%%..%%...%%%%...%%%%%..
.%%..%%..%%..%%..%%%.%%..%%..%%..%%......%%..%%..%%..%%..%%..%%..%%..%%.
.%%......%%%%%%..%%.%%%..%%......%%%%....%%%%%...%%..%%..%%%%%%..%%%%%..
.%%..%%..%%..%%..%%..%%..%%..%%..%%......%%..%%...%%%%...%%..%%..%%..%%.
..%%%%...%%..%%..%%..%%...%%%%...%%%%%%..%%..%%....%%....%%..%%..%%..%%.
........................................................................
      New datasets downloading:
  https://cancervar.wglab.org/databases/
Thanks for using CancerVar!

Need OPAI(Oncogenic Prioritization by Artificial Intelligence)?
Please check https://github.com/WGLab/CancerVar

Report bugs to leequan@gmail.com;
CancerVar homepage: <https://CancerVar.wglab.org>
=============================================================================

```
Then you can check the result file `example/FDA.hg19_multianno.txt.cancervar`.

```
qli@cm-a1-n005|:~/CancerVar-master> head -2 example/FDA.hg19_multianno.txt.cancervar
#Chr    Start   End     Ref     Alt     Ref.Gene        Func.refGene    ExonicFunc.refGene      Gene.ensGene    avsnp147        AAChange.ensGene       AAChange.refGene        clinvar: Clinvar         CancerVar: CancerVar and Evidence      Freq_ExAC_ALL   Freq_esp6500siv2_all   Freq_1000g2015aug_all   Freq_gnomAD_genome_ALL  CADD_raw        CADD_phred      SIFT_score      GERP++_RS       phastCons20way_mammalian       dbscSNV_ADA_SCORE       dbscSNV_RF_SCORE        Interpro_domain AAChange.knownGene      MetaSVM_score   Freq_gnomAD_genome_POPs        OMIM    Phenotype_MIM   OrphaNumber     Orpha   Pathway Therap_list     Diag_list       Prog_list       Polyphen2_HDIV_score   FATHMM_score    MetaLR_score    MutationAssessor_score  cosmic91        icgc28  Otherinfo
9       133748283       133748283       C       T       ABL1    exonic  nonsynonymous SNV       ENSG00000097007 rs121913459     ENSG00000097007:ENST00000318560:exon6:c.C944T:p.T315I,ENSG00000097007:ENST00000372348:exon6:c.C1001T:p.T334I   ABL1:NM_005157:exon6:c.C944T:p.T315I,ABL1:NM_007313:exon6:c.C1001T:p.T334I     clinvar: Pathogenic/Likely_pathogenic    CancerVar: 8#Tier_II_potential EVS=[1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1]       .       .       .       .       6.838   33      0.002   6.04    0.992   .       .       Protein kinase domain;Protein kinase-like domain;Serine-threonine/tyrosine-protein kinase catalytic domain;Tyrosine-protein kinase, catalytic domain ABL1:uc004bzv.3:exon6:c.C1001T:p.T334I,ABL1:uc004bzw.3:exon6:c.C944T:p.T315I     -0.417  AFR:.,AMR:.,EAS:.,FIN:.,NFE:.,OTH:.,ASJ:.     189980   .       .       .       ~ko04012 ErbB signaling pathway~ko04014 Ras signaling pathway~ko04110 Cell cycle~ko04360 Axon guidance~ko04722 Neurotrophin signaling pathway~ko05130 Pathogenic Escherichia coli infection~ko05131 Shigellosis~ko05200 Pathways in cancer~ko05206 MicroRNAs in cancer~ko05220 Chronic myeloid leukemia~ko05416 Viral myocarditis ~   308,309,155,,300,301,302,303,304,305,306,307,299       .       .       0.999   -0.02   0.314   0.565   ID=COSV59323790;OCCURENCE=66(haematopoietic_and_lymphoid_tissue),2(skin)      .

```
The result is tab-delimited,you can import this file into Excel. The colunm of "CancerVar: CancerVar and Evidence" give the CancerVar interpretation result with all the criteria. such as this variant as Tier_II_potential with score 8:

CancerVar: 8#Tier_II_potential EVS=[1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1]


## OPAI: 
now we can move to Oncogenic Prioritization by Artificial Intelligence.

OPAI firstly call **feature_preprocess.py** to process the features coding from CancerVar and Annovar output, then call **opai_predictor.py** to predict the oncogenicity.

The OPAI scripts are in the **scripts** folder of **“OPAI”**:
- feature_preprocess.py:
   - preprocessing the ANNOVAR data and CancerVar output to generate OPAI input;
- opai_predictor.py:
   - predicting the oncogenicity of a variant


 There are two trained models for prediction in OPAI, located in the folder of **"saves"**:
- Ensemble-based model:
   - both clinical evidence score and 23 pre-computed in silico scores are taken as input of the model;
   - model file: `ensemble.pt`
- Evidence-based model:
    - only clinical evidence score are taken as input of the model, this is useful for case of a lot or even all the missing values in 23 pre-computed in silico scores.
    -  model file: `evs.pt`

 Users can specify the model by using the `-m ensemble ` or `-m evs` option and then following the `-d model_file_location` option.


 - using Ensemble-based model
```
   python3.6 OPAI/scripts/feature_preprocess.py -a example/FDA.hg19_multianno.txt.grl_p -c  example/FDA.hg19_multianno.txt.cancervar -m ensemble -n 5 -d OPAI/saves/nonmissing_db.npy -o example/FDA.hg19_multianno.txt.cancervar.ensemble.csv

   python3.6 OPAI/scripts/opai_predictor.py -i  example/FDA.hg19_multianno.txt.cancervar.ensemble.csv -m ensemble -c OPAI/saves/ensemble.pt -d cpu -v example/FDA.hg19_multianno.txt.cancervar -o example/FDA.hg19_multianno.txt.cancervar.ensemble.pred
```
The predicted oncogenicity are in the (last)column of **"ensemble_score"** in file `example/FDA.hg19_multianno.txt.cancervar.ensemble.pred`.

- using Evidence-based model
```
   python3.6 OPAI/scripts/feature_preprocess.py -a example/FDA.hg19_multianno.txt.grl_p -c  example/FDA.hg19_multianno.txt.cancervar -m evs -n 5 -d OPAI/saves/nonmissing_db.npy -o example/FDA.hg19_multianno.txt.cancervar.evs.csv

   python3.6 OPAI/scripts/opai_predictor.py -i  example/FDA.hg19_multianno.txt.cancervar.evs.csv -m evs -c OPAI/saves/evs.pt -d cpu -v example/FDA.hg19_multianno.txt.cancervar -o example/FDA.hg19_multianno.txt.cancervar.evs.pred

```
The predicted oncogenicity are in the (last)column of **"evs_score"** in file `example/FDA.hg19_multianno.txt.cancervar.evs.pred`.

