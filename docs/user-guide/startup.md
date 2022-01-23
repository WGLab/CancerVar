## Quick guide to CancerVar

For beginners, the easiest way to use CancerVar is to annotate a VCF file by ANNOVAR and generate a multianno file, then CancerVar will use this file directly as input  for clinical interpretation.

1. You need install Python >=3.6.
2. You need install [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
3. You need download [CancerVar](https://github.com/WGLab/CancerVar/archive/master.zip) or git clone using URL  https://github.com/WGLab/CancerVar.git
4. Downlad all necessary datasets from ANNOVAR( CancerVar also will download automatically)
5. Prepare your input file, it can be VCF file or AVinput as ANNOVAR input,actually the input file only need information of Chr,Position start,Position end, Reference_allele,Alternative_allele 
6. Run the CancerVar!



### EXAMPLE OF CANCERVAR

    ./CancerVar.py -c config.ini  # Run the examples in config.ini

    ./CancerVar.py  -b hg19 -i your_input  --input_type=VCF  -o your_output


## Quick guide to OPAI

After running CancerVar correclty and getting the output files of "*.cancervar" and "*.grl_p",we are ready to run Oncogenic Prioritization by Artificial Intelligenc.


OPAI has currently only been tested with Python 3.6+, and requires four Python modules to be installed and in path. These are numpy https://numpy.org, pandas https://pandas.pydata.org , scikit-learn https://scikit-learn.org and pytorch https://pytorch.org.


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
 


## Additional resources

A step-by-step protocol on using CancerVar & OPAI  is available at the manual document.

