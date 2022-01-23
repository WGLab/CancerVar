# CancerVar & OPAI
Clinical interpretation of Cancer somatic Variants(CancerVar) and Oncogenic Prioritization by Artificial Intelligence(OPAI)

## HOW DOES IT WORK

CancerVar takes either pre-annotated files, or unannotated input files in VCF format or ANNOVAR input format, where each line corresponds to one genetic variant; CancerVar will call ANNOVAR to generate necessary annotations.
In the output, based on all 12 pieces of evidence, each variant will be assigned as "Tier_I_strong", "Tier_II_potential", "Tier_IV_benign" and "Tier_III_Uncertain" by rules specified in the AMP/ASCO/CAP 2017 guidelines.

OPAI takes 12 clinical evidence prediction scores from CancerVar and 23 pre-computed in silico scores predicted by other computational tools as input from ANNOVAR, and predicts oncogenicity scores by a semi-supervised deep-learning model.

CanverVar and OPAI are Python based scripts. The user need to run CancerVar firstly as **step 1** to get clinical evidence-based interpretation results and then run OPAI as **step 2** if they want to get the deep-learning model-based prediction results.

## CancerVar(step 1)

#### SYNOPSIS

CancerVar.py [options]

#### WHAT DOES IT DO

CanverVar is a python script for cancer variant interpretation of clinical significance.

#### PREREQUISITE

1. You need install **Python >=3.6**
2. You need install **[ANNOVAR]**(http://annovar.openbioinformatics.org/en/latest/) version >=  2016-02-01.
3. Most of the datases can be downloaded automatically.
4. Some updated datasets(c**osmic and icgc**) for Annovar:  [https://cancervar.wglab.org/databases/](https://cancervar.wglab.org/databases/) (download and gunzip, put in the Annovar db folder)
5. Please use the updated files, outdated files will bring some problems of running CancerVar.


#### OPTIONS of CancerVar script

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


#### EXAMPLE of CancerVar
```
    python3.6 ./CancerVar.py -c config.ini  # Run the examples in config.ini
    python3.6 ./CancerVar.py  -b hg19 -i your_input  --input_type=VCF  -o your_output
    python3.6 ./CancerVar.py  -b hg19 -i example/FDA_hg19.av -o example/FDA
```

## OPAI(step 2)

After running CancerVar correclty and getting the output files of **"*.cancervar"** and **"*.grl_p"**,we are ready to run Oncogenic Prioritization by Artificial Intelligenc.

### WHAT AND HOW DOES IT DO

OPAI is a python script for Oncogenic Prioritization by Artificial Intelligence after CancerVar.
OPAI firstly call **feature_preprocess.py** to process the features coding from CancerVar and Annovar output, then call **opai_predictor.py** to precit the oncogenicity.

The OPAI scripts are in the **scripts** folder of **“OPAI”**:
- feature_preprocess.py: 
   - preprocessing the ANNOVAR data and CancerVar output to generate OPAI input;
- opai_predictor.py: 
   - predicting the oncogenicity of a variant.



#### PREREQUISITE
OPAI has currently only been tested with **Python 3.6+**, and requires four Python modules to be installed and in path. These are **numpy** https://numpy.org, **pandas** https://pandas.pydata.org , **scikit-learn** https://scikit-learn.org and **pytorch** https://pytorch.org. 

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
#### MODELS
 There are two trained models for prediction in OPAI, located in the folder of **"saves"**:
- Ensemble-based model: 
   - both clinical evidence prediction scores and 23 pre-computed in silico scores are taken as input of the model;
   - model file: `ensemble.pt`
- Evidence-based model: 
    - only clinical evidence prediction scores are taken as input of the model, this is useful for case of a lot or even all the missing values in 23 pre-computed in silico scores.
    -  model file: `evs.pt`

 Users can specify the model by using the `-m ensemble ` or `-m evs` option and then following the `-d model_file_location` option.
 
 #### EXAMPLE of OPAI
 After running of `python3.6 ./CancerVar.py  -b hg19 -i example/FDA_hg19.av -o example/FDA`, check files of `example/FDA.hg19_multianno.txt.grl_p` and `example/FDA.hg19_multianno.txt.cancervar`, see if they are generated correctly.

 Then,

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
 
#### OPTIONS OF OPAI SCRIPTS 
- Feature process using `feature_preprocess.py`
```bash
python3.6  OPAI/scripts/feature_preprocess.py -h
usage: feature_preprocess.py [-h] -a ANNOVAR_PATH -c CANCERVAR_PATH [-m METHOD] [-n MISSING_COUNT] -d DATABASE -o OUTPUT

feature creator from cancervar output

optional arguments:
  -h, --help            show this help message and exit
  -a ANNOVAR_PATH, --annovar_path ANNOVAR_PATH
                        the path to annovar file
  -c CANCERVAR_PATH, --cancervar_path CANCERVAR_PATH
                        the path to cancervar file
  -m METHOD, --method METHOD
                        output evs features or ensemble features (option: evs, ensemble)
  -n MISSING_COUNT, --missing_count MISSING_COUNT
                        variant with more than N missing features will be discarded, (default: 5)
  -d DATABASE, --database DATABASE
                        database for feature normalization
  -o OUTPUT, --output OUTPUT
                        the path to output

```

- Prediction using `opai_predictor.py`
```bash
python3.6 OPAI/scripts/opai_predictor.py -h
usage: opai_predictor.py [-h] -i INPUT -v CANCERVAR_PATH [-m METHOD] [-d DEVICE] -c CONFIG -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        the path to input feature
  -v CANCERVAR_PATH, --cancervar_path CANCERVAR_PATH
                        the path to cancervar file
  -m METHOD, --method METHOD
                        use evs features or ensemble features (option: evs, ensemble)
  -d DEVICE, --device DEVICE
                        device used for dl-based predicting (option: cpu, cuda)
  -c CONFIG, --config CONFIG
                        the path to trained model file
  -o OUTPUT, --output OUTPUT
                        the path to output

```

 
 
## Web server
CancerVar:  [http://cancervar.wglab.org](http://cancervar.wglab.org)
The web server provided pre-compiled 13M mutations annotation results and OPAI scores.

## LICENSE

CancerVar and OPAI is free for non-commercial use without warranty. Users need to obtain licenses such as OMIM and ANNOVAR by themselves. Please contact the authors for commercial use.

## REFERENCE

Quan Li, Zilin Ren, Kajia Cao, Marilyn M. Li, Yunyun Zhou and Kai Wang. CancerVar: an Artificial Intelligence empowered platform for clinical interpretation of somatic mutations in cancer.(Under Review,2022)[BioRxiv](https://doi.org/10.1101/2020.10.06.323162)

Quan Li and Kai Wang. InterVar: Clinical interpretation of genetic variants by ACMG-AMP 2015 guideline. The American Journal of Human Genetics 100, 1-14, February 2, 2017,[http://dx.doi.org/10.1016/j.ajhg.2017.01.004](http://dx.doi.org/10.1016/j.ajhg.2017.01.004)

[The  AMP/ASCO/CAP 2017 guidelines ](https://www.ncbi.nlm.nih.gov/pubmed/27993330)
Li MM, Datto M, Duncavage EJ, Kulkarni S, Lindeman NI, Roy S, Tsimberidou AM, Vnencak-Jones CL, Wolff DJ, Younes A, Nikiforova MN.
Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists.

[The  ACMG/CGC 2019 guidelines ](https://www.ncbi.nlm.nih.gov/pubmed/31138931)
2.Mikhail FM, et al. Technical laboratory standards for interpretation and reporting of acquired copy-number abnormalities and copy-neutral loss of heterozygosity in neoplastic disorders: a joint consensus recommendation from the American College of Medical Genetics and Genomics (ACMG) and the Cancer Genomics Consortium (CGC). Genet Med. 2019 Sep;21(9):1903-1916. doi: 10.1038/s41436-019-0545-7.

[The ACMG 2015 guide](http://www.ncbi.nlm.nih.gov/pubmed/25741868)
Richards, S. et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genetics in medicine : official journal of the American College of Medical Genetics 17, 405-424 (2015).

## Acknowledges

Thanks to all who provided bug reports.

