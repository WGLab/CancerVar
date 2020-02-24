# CancerVar
Clinical interpretation of somatic mutations in cancer
........................................................................
..%%%%....%%%%...%%..%%...%%%%...%%%%%%..%%%%%...%%..%%...%%%%...%%%%%..
.%%..%%..%%..%%..%%%.%%..%%..%%..%%......%%..%%..%%..%%..%%..%%..%%..%%.
.%%......%%%%%%..%%.%%%..%%......%%%%....%%%%%...%%..%%..%%%%%%..%%%%%..
.%%..%%..%%..%%..%%..%%..%%..%%..%%......%%..%%...%%%%...%%..%%..%%..%%.
..%%%%...%%..%%..%%..%%...%%%%...%%%%%%..%%..%%....%%....%%..%%..%%..%%.
........................................................................

## SYNOPSIS

CancerVar.py [options]

## WHAT DOES IT DO

CanverVar is a python script for cancer variant interpretation of clinical significance.

## PREREQUISITE

1. You need install Python >=2.6.6.
2. You need install [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) version >=  2016-02-01.
3. You need download other files such as mim2gene.txt from [OMIM](http://www.omim.org/downloads).
4. Please use the updated files(should be generated: >= 2016-09), outdated files will bring problems of InterVar.

## OPTIONS

- -h, --help
show this help message and exit

- --version
show program''s version number and exit

- --config=config.ini
Load your config file. The config file contains all options.

if you use this options,you can ignore all the other options bellow.

- -i INPUTFILE, --input=INPUTFILE
input file of  variants for analysis

- --input_type=AVinput
The input file type, it can be  AVinput(Annovar''sformat),VCF

- --cancer_type=CANCER  
The cancer type, please check the help for the details of cancer type


- -o OUTPUTFILE, --output=OUTPUTFILE
prefix the output file (default:output)

- -b BUILDVER, --buildver=BUILDVER
version of reference genome: hg38, hg19(default)

- -t cancervardb, --database_intervar=cancervardb
The database location/dir for the CancerVar dataset files

- --table_annovar=./table_annovar.pl
The Annovar perl script of table_annovar.pl

- --convert2annovar=./convert2annovar.pl
The Annovar perl script of convert2annovar.pl

- --annotate_variation=./annotate_variation.pl
The Annovar perl script of annotate_variation.pl

-  -d humandb, --database_locat=humandb
The database location/dir for the Annovar annotation datasets


## EXAMPLE
```
    ./CancerVar.py -c config.ini  # Run the examples in config.ini
    ./CancerVar.py  -b hg19 -i your_input  --input_type=VCF  -o your_output
```

## HOW DOES IT WORK

CancerVar takes either pre-annotated files, or unannotated input files in VCF format or ANNOVAR input format, where each line corresponds to one genetic variant; CancerVar will call ANNOVAR to generate necessary annotations.
In the output, based on all 10 pieces of evidence, each variant will be assigned as "pathogenic", "likely pathogenic", "uncertain significance", "likely benign/benign" by rules specified in the AMP/ASCO/CAP 2017 guidelines.

## Web server
CancerVar:  [http://cancervar.wglab.org](http://cancervar.wglab.org)

## LICENSE

CancerVar is free for non-commercial use without warranty. Users need to obtain licenses such as OMIM and ANNOVAR by themselves. Please contact the authors for commercial use.

## REFERENCE


Quan Li,Yunyun Zhou and Kai Wang. CancerVar: a web server for improved evidence-based clinical interpretation of cancer somatic mutations and copy number abnormalities (Under Review,2020)

Quan Li and Kai Wang. InterVar: Clinical interpretation of genetic variants by ACMG-AMP 2015 guideline. The American Journal of Human Genetics 100, 1-14, February 2, 2017,[http://dx.doi.org/10.1016/j.ajhg.2017.01.004](http://dx.doi.org/10.1016/j.ajhg.2017.01.004)

[The  AMP/ASCO/CAP 2017 guidelines ](https://www.ncbi.nlm.nih.gov/pubmed/27993330)
Li MM, Datto M, Duncavage EJ, Kulkarni S, Lindeman NI, Roy S, Tsimberidou AM, Vnencak-Jones CL, Wolff DJ, Younes A, Nikiforova MN.
Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists.

[The  ACMG/CGC 2019 guidelines ](https://www.ncbi.nlm.nih.gov/pubmed/31138931)
2.Mikhail FM, et al. Technical laboratory standards for interpretation and reporting of acquired copy-number abnormalities and copy-neutral loss of heterozygosity in neoplastic disorders: a joint consensus recommendation from the American College of Medical Genetics and Genomics (ACMG) and the Cancer Genomics Consortium (CGC). Genet Med. 2019 Sep;21(9):1903-1916. doi: 10.1038/s41436-019-0545-7.

[The ACMG 2015 guide](http://www.ncbi.nlm.nih.gov/pubmed/25741868)
Richards, S. et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genetics in medicine : official journal of the American College of Medical Genetics 17, 405-424 (2015).

