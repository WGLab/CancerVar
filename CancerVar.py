##################################################################################
# Author: Li Quan (leequan@gmail.com)
# Created Time: 2017-01-13 20:38:23 Friday 
# File Name: CancerVar.py File type: python
# Last Change:.
# Description: python script for Interpretation of Pathogenetic of Cancer Variants
##################################################################################
#!/usr/bin/env python

import string,copy,logging,os,io,re,time,sys,platform,optparse,gzip,csv,glob

prog="CancerVar"

version = """%prog 1.1
Copyright (C) 2019 Wang Genomic Lab
CancerVar is free for non-commercial use without warranty.
Please contact the authors for commercial use.
Written by Quan LI,leequan@gmail.com.
============================================================================
........................................................................
..%%%%....%%%%...%%..%%...%%%%...%%%%%%..%%%%%...%%..%%...%%%%...%%%%%..
.%%..%%..%%..%%..%%%.%%..%%..%%..%%......%%..%%..%%..%%..%%..%%..%%..%%.
.%%......%%%%%%..%%.%%%..%%......%%%%....%%%%%...%%..%%..%%%%%%..%%%%%..
.%%..%%..%%..%%..%%..%%..%%..%%..%%......%%..%%...%%%%...%%..%%..%%..%%.
..%%%%...%%..%%..%%..%%...%%%%...%%%%%%..%%..%%....%%....%%..%%..%%..%%.
........................................................................

"""

usage = """Usage: %prog [OPTION] -i  INPUT -o  OUTPUT ...
       %prog  --config=config.ini ...
"""

description = """=============================================================================
CancerVar
Interpretation of Pathogenic/Benign for cancer variants using python script.
=============================================================================
........................................................................
..%%%%....%%%%...%%..%%...%%%%...%%%%%%..%%%%%...%%..%%...%%%%...%%%%%..
.%%..%%..%%..%%..%%%.%%..%%..%%..%%......%%..%%..%%..%%..%%..%%..%%..%%.
.%%......%%%%%%..%%.%%%..%%......%%%%....%%%%%...%%..%%..%%%%%%..%%%%%..
.%%..%%..%%..%%..%%..%%..%%..%%..%%......%%..%%...%%%%...%%..%%..%%..%%.
..%%%%...%%..%%..%%..%%...%%%%...%%%%%%..%%..%%....%%....%%..%%..%%..%%.
........................................................................

"""
end = """=============================================================================
........................................................................
..%%%%....%%%%...%%..%%...%%%%...%%%%%%..%%%%%...%%..%%...%%%%...%%%%%..
.%%..%%..%%..%%..%%%.%%..%%..%%..%%......%%..%%..%%..%%..%%..%%..%%..%%.
.%%......%%%%%%..%%.%%%..%%......%%%%....%%%%%...%%..%%..%%%%%%..%%%%%..
.%%..%%..%%..%%..%%..%%..%%..%%..%%......%%..%%...%%%%...%%..%%..%%..%%.
..%%%%...%%..%%..%%..%%...%%%%...%%%%%%..%%..%%....%%....%%..%%..%%..%%.
........................................................................
Thanks for using CancerVar!
Report bugs to leequan@gmail.com;
CancerVar homepage: <https://CancerVar.wglab.org>
=============================================================================
"""



if platform.python_version()< '3.0.0' :
    import ConfigParser
else:
    import configparser

paras = {}

def ConfigSectionMap(config,section):
    global paras
    options = config.options(section)
    for option in options:
        try:
            paras[option] = config.get(section, option)
            if paras[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            paras[option] = None
    return

user_evidence_dict={}


class myGzipFile(gzip.GzipFile):
    def __enter__(self):
        if self.fileobj is None:
            raise ValueError("I/O operation on closed GzipFile object")
        return self

    def __exit__(self, *args):
        self.close()


#begin read some important datsets/list firstly;
lof_genes_dict={}
mim2gene_dict={}
mim2gene_dict2={}
exclude_snps_dict={}
mim_pheno_dict={}
mim_orpha_dict={}
orpha_dict={}
knownGeneCanonical_dict={}
knownGeneCanonical_st_dict={}
knownGeneCanonical_ed_dict={}

add_markers_dict={}
cancervar_markers_dict={}
civ_markers_dict={}
cancer_pathway_dict={}
cancers_gene_dict={}
cancers_types_dict={}

cancervar_d=[]

def flip_ACGT(acgt):
    nt="";
    if acgt=="A":
        nt="T"
    if acgt=="T":
        nt="A"
    if acgt=="C":
        nt="G"
    if acgt=="G":
        nt="C"
    if acgt=="N":
        nt="N"
    if acgt=="X":
        nt="X"
    return(nt)

def read_datasets():
#0. read the user specified evidence file
    if os.path.isfile(paras['evidence_file']):
        try:
            fh=open(paras['evidence_file'], "r")
            strs = fh.read()
            for line2 in strs.split('\n'):
                cls2=line2.split('\t')
                if len(cls2)>1:
                    keys=cls2[0]+"_"+cls2[1]+"_"+cls2[2]+"_"+cls2[3]+"_"+cls2[4]
                    keys=re.sub("[Cc][Hh][Rr]","",keys)
                    user_evidence_dict[keys]=cls2[5].upper()
                    print("%s %s\n" %(keys,user_evidence_dict[keys]))
        except IOError:
            print("Error: can\'t read the user specified evidence file %s" % paras['evidence_file'])
        else:
            fh.close()

#1.LOF gene list
    try:
        fh = open(paras['lof_genes'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                lof_genes_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the LOF genes file %s" % paras['lof_genes'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()


#1. cancervar_markers
    global cancervar_d
    try:
        with open(paras['cancervar_markers']) as fh:
            reader = csv.reader(fh, delimiter="\t")
            cancervar_d = list(reader)

    except IOError:
        print("Error: can\'t read the cancervar_markers file %s" % paras['cancervar_markers'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()

        
        print len(cancervar_d)
        for ii in range(0,len(cancervar_d)):
            gene=cancervar_d[ii][0]
            #print ii,gene
            mut_type=cancervar_d[ii][1]
            mut=cancervar_d[ii][2]
            mut_types=mut_type.split(';');  # BIA CNA EXPR FUS MUT
            list_mut_size=len(mut_types)
            mut_alts=mut.split(';');
            for jj in range(0,list_mut_size):
                mut_list=mut_alts[jj].split(':')
                if(len(mut_list)>1):
                    gene1=mut_list[0]
                    mutb=mut_list[1]
                else:
                    gene1=gene
                    mutb=mut_alts[jj]
                    #print ii,gene1,mutb
                if(mut_types[jj]=="FUS"):                    
                    mutt="fus_"+re.sub('__','-',mut_alts[jj])
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    #print key,ii,cancervar_markers_dict[key]
                elif(mut_types[jj]=="EXP" or mut_types[jj]=="EXPR"):                    
                    mutt=gene1+"_"+"expr_"+mutb
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                elif(mut_types[jj]=="CNV" or mut_types[jj]=="CNA"):                    
                    mutt=gene1+"_"+"cna_"+mutb
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                elif(mut_types[jj]=="BIA"):                    
                    mutt=gene1+"_"+"bia_"+mutb
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                elif(mut_types[jj]=="OTH"):                    
                    mutt=gene1+"_"+"oth_"+mutb
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                #start from mutb
                elif(mut_types[jj]=="MUT"):                    
                    if(mutb=="any mutation"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                        #print key,ii,cancervar_markers_dict[key]
                    elif(mutb=="any insertion"): 
                        mutt="insertion_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any deletion"): 
                        mutt="deletion_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any frameshift" or mutb=="frameshift"): 
                        mutt="frameshift_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any indel"): 
                        mutt="indel_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any missense"): 
                        mutt="missense_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="mutant"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="mutation" or mutb=="ALTERATION"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="MUTATION"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="oncogenic mutation"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="Oncogenic Mutations"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                        #print key,ii,cancervar_markers_dict[key]
                    elif(mutb=="any nonsense"): 
                        mutt="nonsense_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)

                    # format as " exon(s) 10, 20, 21 any"
                    elif(  re.findall('exon\(s\) ', mutb, flags=re.IGNORECASE)):
                    #print("exon(s) ");
                        str1,str2=mutb.split('exon(s) ',1);
                        if( not  re.findall(',', mutb, flags=re.IGNORECASE)): # for single position
                            poss,mutc=str2.split(' ',1);
                            #mut0=re.sub('"','',mutb)
                            poss_t=int(poss)
                            muts="exon_"+str(poss_t)+"_"+mutc
                            key=gene+'_'+muts
                            default_s=''
                            cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                            #print key,ii,cancervar_markers_dict[key]
                        else: # for multiple positions
                            list_pos=str2.split(',');
                            list_pos_size=len(list_pos)

                            list_tt=list_pos[list_pos_size-1].split(' ')
                            list_tt_size=len(list_tt)

                            mutc=list_tt[list_tt_size-1]
                            pos_end=str(int(list_tt[list_tt_size-2]))
                            if(mutc=="mutation"): mutc="any"
                            muts="exon_"+pos_end+"_"+mutc
                            key=gene+'_'+muts
                            default_s=''
                            cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                            #print key,ii,add_markers_dict[key]
                            for jj in range(0,list_pos_size-1):
                                pos_be=str(int(list_pos[jj]))
                                if(mutc=="mutation"): mutc="any"
                                muts="exon_"+pos_be+"_"+mutc
                                key=gene+'_'+muts
                                default_s=''
                                cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                                #print key,ii,cancervar_markers_dict[key]
                    # " codon(s) 289, 596, 598 any"   codon(s) 132 any
                    elif(  re.findall('codon\(s\) ', mutb, flags=re.IGNORECASE)):
                        #print("codon(s) ");
                        str1,str2=mutb.split('codon(s) ',1);
                        if( not  re.findall(',', mutb, flags=re.IGNORECASE)): # for single position
                            poss,mutc=str2.split(' ',1);
                            #mut0=re.sub('"','',mutb)
                            poss_t=int(poss)
                            muts="codon_"+str(poss_t)+"_"+mutc
                            key=gene+'_'+muts
                            default_s=''
                            cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                            #print key,ii,cancervar_markers_dict[key]
                        else: # for multiple positions
                            list_pos=str2.split(',');
                            list_pos_size=len(list_pos)
 
                            list_tt=list_pos[list_pos_size-1].split(' ')
                            list_tt_size=len(list_tt)
 
                            mutc=list_tt[list_tt_size-1]
                            pos_end=str(int(list_tt[list_tt_size-2]))
                            muts="codon_"+pos_end+"_"+mutc
                            key=gene+'_'+muts
                            default_s=''
                            cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                            #print key,ii,add_markers_dict[key]
                            for jj in range(0,list_pos_size-1):
                                pos_be=str(int(list_pos[jj]))
                                muts="codon_"+pos_be+"_"+mutc
                                key=gene+'_'+muts
                                default_s=''
                                cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                                #print key,ii,cancervar_markers_dict[key]
                    else:
                        mutt=gene1+"_"+mutb
                        key=mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                        #print key,ii,cancervar_markers_dict[key]
                else:
                    key=gene1+"_"+mutb
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    #print key,ii,cancervar_markers_dict[key]
                                        
                
#        print cancervar_markers_dict
# end deal with cancervar_markers





#4. OMIM mim2gene.txt file
    try:
        fh = open(paras['mim2gene'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2)>1:
                cls0=cls2[4].split(',')
                keys=cls0[0]
                mim2gene_dict[keys]=cls2[0]
                keys1=cls2[3]
                keys=keys1.upper()
                mim2gene_dict2[keys]=cls2[0]
    except IOError:
        print("Error: can\'t read the OMIM  file %s" % paras['mim2gene'])
        print("Error: Please download it from http://www.omim.org/downloads")
        sys.exit()
    else:
        fh.close()


#5. read the user specified SNP list, the variants will pass the frequency check.
    if os.path.isfile(paras['exclude_snps']):
        try:
            fh=open(paras['exclude_snps'], "r")
            strs = fh.read()
            for line2 in strs.split('\n'):
                cls2=line2.split('\t')
                if len(cls2)>1:
                    keys=cls2[0]+"_"+cls2[1]+"_"+cls2[2]+"_"+cls2[3]
                    keys=re.sub("[Cc][Hh][Rr]","",keys)
                    exclude_snps_dict[keys]="1"
        except IOError:
            print("Error: can\'t read the user specified SNP list file %s" % paras['exclude_snps'])
        else:
            fh.close()


#6. knownGeneCanonical exon file  # caution the build ver, now it is hg19
    try:
        fh = open(paras['knowngenecanonical'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split(' ')
            if len(cls2)>1:
                keys=cls2[0]
                knownGeneCanonical_dict[keys]=cls2[1]
                knownGeneCanonical_st_dict[keys]=cls2[2]
                knownGeneCanonical_ed_dict[keys]=cls2[3]
                #print("%s %s" %(keys,knownGeneCanonical_dict[keys]))
    except IOError:
        print("Error: can\'t read the knownGeneCanonical  file %s" % paras['knowngenecanonical'])
        print("Error: Please download it from the source website")
        sys.exit()
    else:
        fh.close()


#7. OMIM mim_pheno.txt file

    try:
        fh = open(paras['mim_pheno'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split(' ')
            if len(cls2)>1:
                keys=cls2[0]
                mim_pheno_dict[keys]=cls2[1]
                #print("%s %s" %(keys,mim_pheno_dict[keys]))
    except IOError:
        print("Error: can\'t read the MIM  file %s" % paras['mim_pheno'])
        print("Error: Please download it from CancerVar source website")
        sys.exit()
    else:
        fh.close()




#8. OMIM mim_orpha.txt file
    try:
        fh = open(paras['mim_orpha'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split(' ')
            if len(cls2)>1:
                keys=cls2[0]
                mim_orpha_dict[keys]=cls2[1]
                #print("%s %s" %(keys,mim_orpha_dict[keys]))
    except IOError:
        print("Error: can\'t read the MIM  file %s" % paras['mim_orpha'])
        print("Error: Please download it from CancerVar source website")
        sys.exit()
    else:
        fh.close()

#9.  orpha.txt file
    try:
        fh = open(paras['orpha'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2)>1:
                keys=cls2[0]
                orpha_dict[keys]=cls2[1]
                #print("%s %s" %(keys,mim_orpha_dict[keys]))
    except IOError:
        print("Error: can\'t read the Orpha  file %s" % paras['orpha'])
        print("Error: Please download it from CancerVar source website")
        sys.exit()
    else:
        fh.close()

#10 cancer_pathway=%(database_cancervar)s/cancer_pathway.list
    global cancer_pathway_dict
    try:
        fh = open(paras['cancer_pathway'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                #cancer_pathway_dict[cls2[0]]='1'
                cancer_pathway_dict[cls2[0]]=cls2[1]
    except IOError:
        print("Error: can\'t read the cancer_pathway genes file %s" % paras['cancer_pathway'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()

#11 cancers_genes=%(database_cancervar)s/cancers_genes.list
    global cancers_gene_dict
    try:
        fh = open(paras['cancers_genes'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('\t')
            if len(cls2[0])>1:
                cancers_gene_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the cancers diseases genes file %s" % paras['cancer_pathway'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()

#12 cancers_types=%(database_cancervar)s/cancervar.cancer.types
    global cancers_types_dict
    try:
        fh = open(paras['cancers_types'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2=line2.split('/')
            if len(cls2)>1:
                cancers_types_dict[cls2[1]]=cls2[0];
                #print("%s %s" %(cls2[1],cls2[0]))
    except IOError:
        print("Error: can\'t read the cancers types file %s" % paras['cancers_types'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()


#end read datasets
    return



def check_downdb():
    path=paras['database_locat']
    path=path.strip()
    path=path.rstrip("\/")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        print("Notice: the folder of %s is created!" % path)
    else:
        print("Warning: the folder of %s is already created!" % path)
    ds=paras['database_names']
    ds.expandtabs(1);
    # database_names = refGene 1000g2014oct esp6500siv2_all avsnp147 ljb26_all clinvar_20150629 exac03 hg19_dbscsnv11 dbnsfp31a_interpro rmsk ensGene

    for dbs in ds.split():
        # os.path.isfile(options.table_annovar)
        file_name=dbs
        if dbs=="1000g2015aug":
            file_name="ALL.sites.2015_08"

        dataset_file=paras['database_locat']+"/"+paras['buildver']+"_"+file_name+".txt"
        if dbs != 'rmsk':
            cmd="perl "+paras['annotate_variation']+" -buildver "+paras['buildver']+" -downdb -webfrom annovar "+file_name+" "+paras['database_locat']
        if dbs == 'rmsk':
            cmd="perl "+paras['annotate_variation']+" -buildver "+paras['buildver']+" -downdb "+file_name+" "+paras['database_locat']
        if  not os.path.isfile(dataset_file):
            if dbs=="1000g2015aug":
                file_name="1000g2015aug"
                dataset_file=paras['database_locat']+"/"+paras['buildver']+"_"+file_name+".txt"
                cmd="perl "+paras['annotate_variation']+" -buildver "+paras['buildver']+" -downdb -webfrom annovar "+file_name+" "+paras['database_locat']
            print("Warning: The Annovar dataset file of %s is not in %s,begin to download this %s ..." %(dbs,paras['database_locat'],dataset_file))
            print("%s" %cmd)
            os.system(cmd)
    return

def check_input():
    inputft= paras['inputfile_type']
    if inputft.lower() == 'avinput' :
        return
    if inputft.lower() == 'vcf':
        if os.path.isfile(paras['convert2annovar']):
        #convert2annovar.pl -format vcf4 variantfile > variant.avinput
            cmd="perl "+paras['convert2annovar']+" -format vcf4 "+ paras['inputfile']+"> "+paras['inputfile']+".avinput"
            print("Warning: Begin to convert your vcf file of %s to AVinput of Annovar ..." % paras['inputfile'])
            print("%s" %cmd)
            os.system(cmd)
        else:
            print("Error: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                    % paras['convert2annovar'])
            sys.exit()
    if inputft.lower() == 'vcf_m':
        if os.path.isfile(paras['convert2annovar']):
        #convert2annovar.pl -format vcf4 variantfile > variant.avinput
            cmd="perl "+paras['convert2annovar']+" -format vcf4 "+ paras['inputfile']+" --allsample   --outfile "+ paras['outfile']
            print("Warning: Begin to convert your vcf file with multiple samples of %s to AVinput of Annovar with All.raw.highqc.vcf.<samplename>.avinput..." % paras['inputfile'])
            print("Warning: Please attention that the sample names in VCF file should  contain letters/numners only, otherwise the converting may be failure!")
            print("%s" %cmd)
            os.system(cmd)
        else:
            print("Error: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                    % paras['convert2annovar'])
            sys.exit()
    return

def check_annovar_result():
# table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp144,ljb26_all,CLNSIG,exac03   -operation  g,f,f,f,f,f,f   -nastring . -csvout
    inputft= paras['inputfile_type']
    annovar_options=" --otherinfo "
    #if re.findall('true',paras['otherinfo'], flags=re.IGNORECASE)  :
    #    annovar_options=annovar_options+"--otherinfo "
    if re.findall('true',paras['onetranscript'], flags=re.IGNORECASE) :
        annovar_options=annovar_options+"--onetranscript "

    if not os.path.isfile(paras['table_annovar']):
        print("Error: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                    % paras['table_annovar'])
        sys.exit()
    if inputft.lower() == 'avinput' :
        #cmd="perl "+paras['table_annovar']+" "+paras['inputfile']+" "+paras['database_locat']+" -buildver "+paras['buildver']+" -remove -out "+ paras['outfile']+" -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,dbnsfp30a,clinvar_20160302,exac03,dbscsnv11,dbnsfp31a_interpro,rmsk,ensGene,knownGene  -operation  g,f,f,f,f,f,f,f,f,r,g,g   -nastring ."+annovar_options
        cmd="perl "+paras['table_annovar']+" "+paras['inputfile']+" "+paras['database_locat']+" -buildver "+paras['buildver']+" -remove -out "+ paras['outfile']+" -protocol refGene,ensGene,knownGene,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,dbnsfp30a,dbscsnv11,dbnsfp31a_interpro,clinvar_20190305,cosmic70,icgc21,gnomad_genome  -operation  g,g,g,f,f,f,f,f,f,f,f,f,f,f  -nastring ."+annovar_options
        print("%s" %cmd)
        os.system(cmd)
    if inputft.lower() == 'vcf' :
        cmd="perl "+paras['table_annovar']+" "+paras['inputfile']+".avinput "+paras['database_locat']+" -buildver "+paras['buildver']+" -remove -out "+ paras['outfile']+" -protocol  refGene,ensGene,knownGene,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,dbnsfp30a,dbscsnv11,dbnsfp31a_interpro,clinvar_20190305,cosmic70,icgc21,gnomad_genome  -operation  g,g,g,f,f,f,f,f,f,f,f,f,f,f   -nastring ."+annovar_options
        print("%s" %cmd)
        os.system(cmd)
    if inputft.lower() == 'vcf_m' :
        for f in glob.iglob(paras['outfile']+"*.avinput"):
            print("INFO: Begin to annotate sample file of %s ...." %(f))
            new_outfile=re.sub(".avinput","",f)
            cmd="perl "+paras['table_annovar']+" "+f+" "+paras['database_locat']+" -buildver "+paras['buildver']+" -remove -out "+ new_outfile +" -protocol  refGene,ensGene,knownGene,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,dbnsfp30a,dbscsnv11,dbnsfp31a_interpro,clinvar_20161128,cosmic70,icgc21  -operation  g,g,g,f,f,f,f,f,f,f,f,f,f    -nastring ."+annovar_options
            print("%s" %cmd)
            os.system(cmd)




    return



def check_genes(anvfile):
#check with multiple genes, so one gene by one gene  to annote
    newoutfile=anvfile+".grl_p"
    try:
        fh = open(anvfile, "r")
        fw = open(newoutfile, "w")
        strs = fh.read()
        sum=0
        otherinf_pos=1
        for line in strs.split('\n'):
            cls=line.split('\t')
            if len(cls)>1:
                if sum==0 and re.findall('true',paras['otherinfo'], flags=re.IGNORECASE) :
                    for ii in range(0,len(cls)):
                        if  re.findall('otherinfo',cls[ii], flags=re.IGNORECASE) :
                            otherinf_pos=ii

                gene_name=cls[6]
                if cls[6] == 'Gene.refGene':
                    gene_name='Gene'
#some with multiple genes, so one gene by one gene  to annote
                sum=sum+1
                for gg in gene_name.split(','):
                    if not re.findall('true',paras['otherinfo'], flags=re.IGNORECASE) :
                        line_out=line+"\t"+gg
                    else:
                        line_out=cls[0]
                        for ii in range(1,len(cls)):
                            if ii != otherinf_pos :
                                line_out=line_out+"\t"+cls[ii]
                            if ii == otherinf_pos :
                                line_out=line_out+"\t"+gg+"\t"+cls[ii]
                    if sum >1: line_out=re.sub("^[Cc][Hh][Rr]","",line_out)
                    fw.write("%s\t\n" % line_out)

    except IOError:
        print("Error: can\'t read/write the annovar output file %s %s" % (anvfile,newoutfile))
        sys.exit()
        return
    else:
        pass
        fh.close()
        fw.close()

    return(sum)



def sum_of_list(list):
    sum=0
    for i in list:
        sum=sum+i
    return(sum)

def classfy(CBP,Allels_flgs,cls):
    BPS=["Pathogenic","Likely_pathogenic","Benign/Likely_benign","Uncertain_significance"]
    PAS_out=-1
    BES_out=-1
    BPS_out=3 # BPS=[3]:Uncertain significance
    CBP_sum=0;
    # CBP[0]:Therapeutic(2100) CBP[1]:Diagno  CBP[2]:Progno CBP[3]:Mutation(1100) CBP[4]:Variant_freq CBP[5]:Potential_germ
    # CBP[6]: Populatio(1100)  CBP[7]:Germline dat(2201) CBP[8]:Somatic dat(2100) CBP[9]:Predict_dama(2201) 
    # CBP[10]:  Path(2100)  CBP[11] : Pubs

    #print("Before up/down grade, the sum of CBP %s, PM %s,PP %s,BS %s,BP %s" %(PS_sum,PM_sum,PP_sum,BS_sum,BP_sum));
    #begin process the user's flexible grade  to get the final interpretation

    if os.path.isfile(paras['evidence_file']):
        keys=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
        keys=re.sub("[Cc][Hh][Rr]","",keys)
        try:
            evds=user_evidence_dict[keys] #PS1=1;PM1=1;BA1=1;PVS1 PP BS BP
            for evd in evds.split(';'):
                evd_t=evd.split('=')
                if(len(evd_t)>1 and  re.findall('grade', evd_t[0], flags=re.IGNORECASE) ):
                    #10  104353782   G   A   PVS1=1;PP1=1;PM3=1;grade_PP1=2;
                    if int(evd_t[1])<=3:
                        #print ("%s %s %s " %(keys,evd_t[1],evd_t[0]))
                        if(evd_t[0].find('CBP')!=-1):
                            t=evd_t[0].find('CBP');
                            tt=evd_t[0];
                            tt3=int(tt[t+3:t+4])
                            if(t<len(evd_t[0])-2 and tt3<=12 ): CBP[tt3-1]=int(evd_t[1])
        except KeyError:
            pass
        else:
            pass

    # end process the user's flexible grade
    if CBP[0]==2 and CBP[3]==1 and CBP[6]==1 and CBP[7]==2 and CBP[8]==2 and  CBP[9]==2 and  CBP[10]==2:
        BPS_out=0
    if CBP[0]==1 and CBP[3]==1 and CBP[6]==1 and CBP[7]==2 and CBP[8]>=1 and  CBP[9]==2 and  CBP[10]>=1:
        BPS_out=1

    if CBP[0]==0 and CBP[3]==0 and CBP[6]==0 and CBP[7]==0 and CBP[8]==0 and  CBP[9]==0 and  CBP[10]==0:
        BPS_out=3
    if CBP[0]==0 and CBP[3]==0 and CBP[6]==0 and CBP[7]>=0 and CBP[8]==0 and  CBP[9]<=1 and  CBP[10]==0:
        BPS_out=2
# EVS=[1, None, None, 1, None, None, 1, 2, 2, 2, 2, None]
#     [1, None, None, 1, None, None, 1, 0, 1, 2, 2, None]
    for i in range(0,11):
        CBP_sum=CBP_sum+CBP[i]
    if(CBP_sum<=3):  BPS_out=2 # benign
    if(CBP_sum>3 and CBP_sum<8):  BPS_out=3 # VUS
    if(CBP_sum>=8 and CBP_sum<=11):  BPS_out=1 # potential path
    if(CBP_sum>=12):  BPS_out=0 # strong path


    #CBPout=BPS[BPS_out]
    CBPout=str(CBP_sum)+"_"+BPS[BPS_out]
    return(CBPout)



def check_Thera(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Therapeutic: 
    2 FDA approved or investigational with strong evidence
    1 FDA approved for different tumor type; investigational therapies with some evidence
    0 Cancer genes: none
    0 None  
    # the function of check_Thera  check_Diagno check_Progno are similar, should combine as one function,
      but in future for  different purpose and criteria, leave three functions.
    '''

    Thera=0;
    global cancervar_d;

    level=0 # ABCD
    cls=line.split('\t')
    clstt=cls[Funcanno_flgs['Otherinfo']].split(';')
    cancer_type="Cancer"
    #if (len(clstt[0])>0):
    #    cancer_type=clstt[0]
    try:
        cancer_type=paras['cancer_type']
    except KeyError:
        if (len(clstt[0])>0): cancer_type=clstt[0]
    else:
        pass
 

    gene_tr=cls[Funcanno_flgs['Gene']]
    func=cls[Funcanno_flgs['Func.refGene']]
    exonfunc=cls[Funcanno_flgs['Func.refGene']]
    # ABL1:NM_005157:exon6:c.C944T:p.T315I,ABL1:NM_007313:exon6:c.C1001T:p.T334I
    # AAChange.refGene  only check variant with  AAchange information
    line_tmp=cls[Funcanno_flgs['Gene']]+" "+cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]+" "+cancer_type
    line_tmp2=cls[Funcanno_flgs['AAChange.refGene']]
    #line_tmp2_sp=line_tmp2.split(',')
    line_tmp2_sp=re.split(';|,',line_tmp2)
    exon_func=cls[Funcanno_flgs['ExonicFunc.refGene']]
    marker_key0=gene_tr+"_"+"mut_any"
    for cls0 in line_tmp2_sp:
        cls0_1=cls0.split(':')

        if(len(cls0)>1 and len(line_tmp2_sp)>0 and len(cls0_1)==5 ):
            cls0_1=cls0.split(':')
            gene=cls0_1[0]
            transc=cls0_1[1]
            exon=re.sub('exon','',cls0_1[2])
            cdna_change=re.sub('c.','',cls0_1[3])
            amino_change=re.sub('p.','',cls0_1[4])
            ltt=len(amino_change)
            codon_num=amino_change[1:ltt-1]
            #print gene, transc,exon,cdna_change,amino_change,cancer_type 
            marker_key=gene+"_"+amino_change
            marker_key0=gene+"_"+"mut_any"
            marker_key1=gene+"_"+"exon_"+exon+"_any"
            marker_key2=gene+"_"+"codon_"+codon_num+"_any"
            marker_key00=""       
            marker_key11=""       
            marker_key22=""       
            if(cls[Funcanno_flgs['ExonicFunc.refGene']]=="nonsynonymous SNV"):
                marker_key00=gene+"_"+"missense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_missense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_missense"
                
            if exon_func.find("frameshift")>=0 and exon_func.find("nonframe")<0 : 
                marker_key00=gene+"_"+"frameshift_any"
                marker_key11=gene+"_"+"exon_"+exon+"_frameshift"
                marker_key22=gene+"_"+"codon_"+codon_num+"_frameshift"

            if exon_func.find("stopgain")>=0 or exon_func.find("stoploss")>=0:  
                marker_key00=gene+"_"+"nonsense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_nonsense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_nonsense"

            if exon_func.find("deletion")>=0 :  
                marker_key00=gene+"_"+"deletion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_deletion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_deletion"
            if exon_func.find("insertion")>=0 :  
                marker_key00=gene+"_"+"insertion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_insertion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_insertion"
 
            #deletion frameshift indel missense nonsense_any
            # nonsynonymous SNV;stopgain/stoploss(nonsense);frameshift insertion/deletion/substitution; nonframeshift insertion/deletion/substitution;
            default_s=""
            add_list=cancervar_markers_dict.get(marker_key,default_s)+cancervar_markers_dict.get(marker_key0,default_s)+cancervar_markers_dict.get(marker_key00,default_s)+cancervar_markers_dict.get(marker_key1,default_s)+cancervar_markers_dict.get(marker_key11,default_s)+cancervar_markers_dict.get(marker_key2,default_s)+cancervar_markers_dict.get(marker_key22,default_s)
            cgi_list=add_list
            #print gene,add_list
            level_AB="0"
            level_CD="0"
            level=0 # ABCD
            for i in cgi_list.split(","):
                #print i
                if(len(i)>0):
                    pos=int(i)
                    if(cancervar_d[pos][9]=="Therapeutic"):
                        t_level=cancervar_d[pos][10]
                        #print gene,pos,t_level,cancervar_d[pos][10],cancervar_d[pos][9]
                        if(t_level=='A'): level=2;
                        if(t_level=='B'): level=2;
                        if(t_level=='C'): level=1;
                        if(t_level=='D'): level=1;
                        if( not  re.findall(cancer_type, cancervar_d[pos][8], flags=re.IGNORECASE)   ):
                            if(t_level=='A' or t_level=='B'):  # did not find the specific type, decrease level
                                level=1;
 
            #print gene,level,"Therapeutic"
    Thera=level 

       
    return(Thera)
    #print line_tmp

def check_Diagno(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Diagnostic: 
    In PG or reported evidence with consensus
    not in PG but with convincing published data
    Cancer genes: none
    None
    ''' 
    Diagno=0
    global cancervar_d;
    level=0 # ABCD
    cls=line.split('\t')
    clstt=cls[Funcanno_flgs['Otherinfo']].split(';')
    cancer_type="Cancer"
    #if (len(clstt[0])>0):
    #    cancer_type=clstt[0]
    try:
        cancer_type=paras['cancer_type']
    except KeyError:
        if (len(clstt[0])>0): cancer_type=clstt[0]
    else:
        pass
 

    gene_tr=cls[Funcanno_flgs['Gene']]
    func=cls[Funcanno_flgs['Func.refGene']]
    exonfunc=cls[Funcanno_flgs['Func.refGene']]
    # ABL1:NM_005157:exon6:c.C944T:p.T315I,ABL1:NM_007313:exon6:c.C1001T:p.T334I
    # AAChange.refGene  only check missense with AAchange
    line_tmp=cls[Funcanno_flgs['Gene']]+" "+cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]+" "+cancer_type
    line_tmp2=cls[Funcanno_flgs['AAChange.refGene']]
    #line_tmp2_sp=line_tmp2.split(',')
    line_tmp2_sp=re.split(';|,',line_tmp2)
    exon_func=cls[Funcanno_flgs['ExonicFunc.refGene']]
    marker_key0=gene_tr+"_"+"mut_any"
    for cls0 in line_tmp2_sp:
        cls0_1=cls0.split(':')
        if(len(cls0)>1 and len(line_tmp2_sp)>0 and len(cls0_1)==5 ):
            cls0_1=cls0.split(':')
            gene=cls0_1[0]
            transc=cls0_1[1]
            exon=re.sub('exon','',cls0_1[2])
            cdna_change=re.sub('c.','',cls0_1[3])
            amino_change=re.sub('p.','',cls0_1[4])
            ltt=len(amino_change)
            codon_num=amino_change[1:ltt-1]
            #print gene, transc,exon,cdna_change,amino_change,cancer_type 
            marker_key=gene+"_"+amino_change
            marker_key0=gene+"_"+"mut_any"
            marker_key1=gene+"_"+"exon_"+exon+"_any"
            marker_key2=gene+"_"+"codon_"+codon_num+"_any"
            marker_key00=""       
            marker_key11=""       
            marker_key22=""       
            if(cls[Funcanno_flgs['ExonicFunc.refGene']]=="nonsynonymous SNV"):
                marker_key00=gene+"_"+"missense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_missense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_missense"
                
            if exon_func.find("frameshift")>=0 and exon_func.find("nonframe")<0 : 
                marker_key00=gene+"_"+"frameshift_any"
                marker_key11=gene+"_"+"exon_"+exon+"_frameshift"
                marker_key22=gene+"_"+"codon_"+codon_num+"_frameshift"

            if exon_func.find("stopgain")>=0 or exon_func.find("stoploss")>=0:  
                marker_key00=gene+"_"+"nonsense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_nonsense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_nonsense"

            if exon_func.find("deletion")>=0 :  
                marker_key00=gene+"_"+"deletion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_deletion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_deletion"
            if exon_func.find("insertion")>=0 :  
                marker_key00=gene+"_"+"insertion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_insertion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_insertion"
 
            #deletion frameshift indel missense nonsense_any
            # nonsynonymous SNV;stopgain/stoploss(nonsense);frameshift insertion/deletion/substitution; nonframeshift insertion/deletion/substitution;
            default_s=""
            add_list=cancervar_markers_dict.get(marker_key,default_s)+cancervar_markers_dict.get(marker_key0,default_s)+cancervar_markers_dict.get(marker_key00,default_s)+cancervar_markers_dict.get(marker_key1,default_s)+cancervar_markers_dict.get(marker_key11,default_s)+cancervar_markers_dict.get(marker_key2,default_s)+cancervar_markers_dict.get(marker_key22,default_s)
            cgi_list=add_list
            #print gene,add_list
            level_AB="0"
            level_CD="0"
            level=0 # ABCD
            for i in cgi_list.split(","):
                #print i
                if(len(i)>0):
                    pos=int(i)
                    if(cancervar_d[pos][9]=="Diagnostic"):
                        t_level=cancervar_d[pos][10]
                        #print gene,pos,t_level,cancervar_d[pos][10],cancervar_d[pos][9]
                        if(t_level=='A'): level=2;
                        if(t_level=='B'): level=2;
                        if(t_level=='C'): level=1;
                        if(t_level=='D'): level=1;
                        if( not  re.findall(cancer_type, cancervar_d[pos][8], flags=re.IGNORECASE)   ):
                            if(t_level=='A' or t_level=='B'):  # did not find the specific type, decrease level
                                level=1;
 
            #print gene,level
            

    Diagno=level;

    return(Diagno)

def check_Progno(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Prognostic: 
    In PG or reported evidence with consensus
    not in PG but with convincing published data
    Cancer genes: none
    None
    ''' 
    Progno=0
    global cancervar_d;
    level=0 # ABCD
    cls=line.split('\t')
    clstt=cls[Funcanno_flgs['Otherinfo']].split(';')
    cancer_type="Cancer"
    #if (len(clstt[0])>0):
    #    cancer_type=clstt[0]
    try:
        cancer_type=paras['cancer_type']
    except KeyError:
        if (len(clstt[0])>0): cancer_type=clstt[0]
    else:
        pass
 

    gene_tr=cls[Funcanno_flgs['Gene']]
    func=cls[Funcanno_flgs['Func.refGene']]
    exonfunc=cls[Funcanno_flgs['Func.refGene']]
    # ABL1:NM_005157:exon6:c.C944T:p.T315I,ABL1:NM_007313:exon6:c.C1001T:p.T334I
    # AAChange.refGene  only check missense with AAchange
    line_tmp=cls[Funcanno_flgs['Gene']]+" "+cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]+" "+cancer_type
    line_tmp2=cls[Funcanno_flgs['AAChange.refGene']]
    #line_tmp2_sp=line_tmp2.split(',')
    line_tmp2_sp=re.split(';|,',line_tmp2)
    exon_func=cls[Funcanno_flgs['ExonicFunc.refGene']]
    marker_key0=gene_tr+"_"+"mut_any"
    for cls0 in line_tmp2_sp:
        cls0_1=cls0.split(':')
        if(len(cls0)>1 and len(line_tmp2_sp)>0 and len(cls0_1)==5 ):
            cls0_1=cls0.split(':')
            gene=cls0_1[0]
            transc=cls0_1[1]
            exon=re.sub('exon','',cls0_1[2])
            cdna_change=re.sub('c.','',cls0_1[3])
            amino_change=re.sub('p.','',cls0_1[4])
            ltt=len(amino_change)
            codon_num=amino_change[1:ltt-1]
            #print gene, transc,exon,cdna_change,amino_change,cancer_type 
            marker_key=gene+"_"+amino_change
            marker_key0=gene+"_"+"mut_any"
            marker_key1=gene+"_"+"exon_"+exon+"_any"
            marker_key2=gene+"_"+"codon_"+codon_num+"_any"
            marker_key00=""       
            marker_key11=""       
            marker_key22=""       
            if(cls[Funcanno_flgs['ExonicFunc.refGene']]=="nonsynonymous SNV"):
                marker_key00=gene+"_"+"missense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_missense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_missense"
                
            if exon_func.find("frameshift")>=0 and exon_func.find("nonframe")<0 : 
                marker_key00=gene+"_"+"frameshift_any"
                marker_key11=gene+"_"+"exon_"+exon+"_frameshift"
                marker_key22=gene+"_"+"codon_"+codon_num+"_frameshift"

            if exon_func.find("stopgain")>=0 or exon_func.find("stoploss")>=0:  
                marker_key00=gene+"_"+"nonsense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_nonsense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_nonsense"

            if exon_func.find("deletion")>=0 :  
                marker_key00=gene+"_"+"deletion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_deletion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_deletion"
            if exon_func.find("insertion")>=0 :  
                marker_key00=gene+"_"+"insertion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_insertion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_insertion"
 
            #deletion frameshift indel missense nonsense_any
            # nonsynonymous SNV;stopgain/stoploss(nonsense);frameshift insertion/deletion/substitution; nonframeshift insertion/deletion/substitution;
            default_s=""
            add_list=cancervar_markers_dict.get(marker_key,default_s)+cancervar_markers_dict.get(marker_key0,default_s)+cancervar_markers_dict.get(marker_key00,default_s)+cancervar_markers_dict.get(marker_key1,default_s)+cancervar_markers_dict.get(marker_key11,default_s)+cancervar_markers_dict.get(marker_key2,default_s)+cancervar_markers_dict.get(marker_key22,default_s)
            cgi_list=add_list
            #print gene,add_list
            level_AB="0"
            level_CD="0"
            level=0 # ABCD
            for i in cgi_list.split(","):
                #print i
                if(len(i)>0):
                    pos=int(i)
                    if(cancervar_d[pos][9]=="Prognostic"):
                        t_level=cancervar_d[pos][10]
                        #print gene,pos,t_level,cancervar_d[pos][10],cancervar_d[pos][9]
                        if(t_level=='A'): level=2;
                        if(t_level=='B'): level=2;
                        if(t_level=='C'): level=1;
                        if(t_level=='D'): level=1;
                        if( not  re.findall(cancer_type, cancervar_d[pos][8], flags=re.IGNORECASE)   ):
                            if(t_level=='A' or t_level=='B'):  # did not find the specific type, decrease level
                                level=1;
 
            #print gene,level



    Progno=level

    return(Progno)

def check_Mut(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Mutation type:
    1 Activating, LOF (missense, nonsense, indel, splicing), CNVs, fusions
    1 Activating, LOF (missense, nonsense, indel, splicing), CNVs, fusions
    0 Functionally unknown; mostly missense, in-frame indels; less commonly,other types
    0 Functionally benign or unknown; mostly missense; less commonly, other types
    '''

    cls=line.split('\t')
    funcs_tmp=["nonsynonymous","missense","nonsense","frameshift","splic","stopgain","stoplost","CNV","fusion"]
    funcs_tmp2="nonframe"
    line_tmp=cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]
    Mut=0 # 0 for Tire3/4; 1 for Tire 1/2
    VS_t1=0
    VS_t2=0
    VS_t3=0
    dbscSNV_cutoff=0.6    #either score(ada and rf) >0.6 as splicealtering
    # Funcanno_flgs={'Func.refGene':0,'ExonicFunc.refGene':0
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 and line_tmp.find(funcs_tmp2)<0 :
            VS_t1=1
            break
    try:
        if lof_genes_dict[ cls[Funcanno_flgs['Gene']] ] == '1' :
            VS_t2=1
    except KeyError:
        VS_t2=0
    else:
        pass
    # begin check the site in  affect the splicing
    try:
        if float(cls[Funcanno_flgs['dbscSNV_RF_SCORE']])>dbscSNV_cutoff or float(cls[Funcanno_flgs['dbscSNV_ADA_SCORE']])>dbscSNV_cutoff:
            VS_t3=1
    except ValueError:
        pass
    else:
        pass
    if VS_t1 !=0 and VS_t2 != 0 :
        Mut=1
    if VS_t3 !=0 and VS_t2 != 0:
        Mut=1
    #print "mut=",Mut
    return (Mut)



def check_VF(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Variant frequencies
    1 Mostly mosaic
    1 Mostly mosaic
    0 Mosaic or nonmosaic
    0 Mostly nonmosaic (VAF, approximately 50% or 100%)
    '''
    VF=0;
    return(VF)
def check_PotG(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Potential germline
    Mostly nonmosaic (VAF approximately 50% or 100%)
    Mostly nonmosaic (VAF approximately 50% or 100%)
    Mostly nonmosaic (VAF approximately 50% or 100%)
    Mostly nonmosaic (VAF, approximately 50% or 100%)
    '''
    PotG=0
    return (PotG)

def check_PopD(line,Freqs_flgs,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Population database: ESP, dbSNP, 1000Genome, ExAC, gnomad
    1  Absent or extremely low MAF
    1  Absent or extremely low MAF
    1  Absent or extremely low MAF
    0  MAF>0.5% in the general population; or high MAF in some ethnic populations
    
    '''
    MAF_cutoff=0.005
    PopD=0;
    cls=line.split('\t')
    Freqs_3pops={'1000g2015aug_all':0,'esp6500siv2_all':0,'ExAC_ALL':0,'gnomAD_genome_ALL':0}
    # Freqs_flgs
 
    tt=1;
    for key in Freqs_flgs.keys():
        if(cls[Freqs_flgs[key]]!='.'): # absent in all  controls
            tt=tt*0;
    if tt==1:
        PopD=1

    for key in Freqs_3pops.keys():
        try:
            #if (cls[Freqs_flgs[key]]!='.' and float(cls[Freqs_flgs[key]])>0.01): PopD=0 #  MAF>1%
            if (cls[Freqs_flgs[key]]!='.' and float(cls[Freqs_flgs[key]])<MAF_cutoff): PopD=0  #  extremely low MAF

        except ValueError:
            pass
        else:
            pass

    #print "PopD=",PopD,cls[Freqs_flgs['1000g2015aug_all']],cls[Freqs_flgs['esp6500siv2_all']],cls[Freqs_flgs['ExAC_ALL']]
    return (PopD)

def check_GermD(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Germline database: HGMD, ClinVar
    2 May or may not be present
    2 May or may not be present
    0 Absent or downgraded from pathogenic to VUS
    1 Absent or present but downgraded to benign/likely benign
    '''
    GermD=0
    cls=line.split('\t')

    line_tmp2=cls[Funcanno_flgs['CLNSIG']]
    if line_tmp2.find("enign")<0 and line_tmp2.find("athogenic")>=0:
        GermD=2
    if line_tmp2.find("ikely benign")>=0 or line_tmp2.find("enign")>=0:
        GermD=0
    if line_tmp2.find("enign")>=0 and line_tmp2.find("athogenic")>=0:
        GermD=1
    if line_tmp2.find("ncertain significance")>=0 :
        GermD=0
    #print "GermD=",GermD,cls[Funcanno_flgs['CLNSIG']]
    return(GermD)


def check_SomD(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Somatic database: COSMIC, My Cancer Genome, TCGA
    2 Most likely present
    1 Likely present
    0 Absent or present without association to specific tumors (potential germline VUS); present but in very few cases
    0 Absent or present without association to specific tumors (potential rare germline polymorphism)
    ''' # cosmic70    ID=COSM12560;OCCURENCE=60(haematopoietic_and_lymphoid_tissue)
    #ICGC_Id ICGC_Occurrence MU31370893  COCA-CN|1|187|0.00535,PRAD-CA|1|124|0.00806,SKCA-BR|1|66|0.01515,MELA-AU|2|183|0.01093
    SomD=0;
    cls=line.split('\t')

    if cls[Funcanno_flgs['cosmic70']]!="." or cls[Funcanno_flgs['ICGC_Id']]!=".":
        SomD=1 
    if cls[Funcanno_flgs['cosmic70']]=="." and cls[Funcanno_flgs['ICGC_Id']]==".":
        SomD=0
    if cls[Funcanno_flgs['cosmic70']]!="." and cls[Funcanno_flgs['ICGC_Id']]!=".":
        SomD=2 
    return(SomD)
    


def check_PreP(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Predictive software: SIFT, PolyPhen2, MutTaster, CADD, MetaSVM
   2 Mostly damaging; information to be used for reference only
   2 Mostly damaging; information to be used for reference only
   0 Variable; information to be used for reference only
   1 Mostly benign; information to be used for reference only
    '''
    # MetaSVM SIFT Polyphen2_HDIV MetaLR FATHMM  MutationTaster
    sift_cutoff=0.05 #SIFT_score,SIFT_pred, The smaller the score the more likely the SNP has damaging effect
    metasvm_cutoff=0 # greater scores indicating more likely deleterious effects
    cutoff_conserv=2 # for GERP++_RS

    dam=0;
    var=0;
    ben=0;
    PreP=0;

    cls=line.split('\t')
    try:
        if float(cls[Funcanno_flgs['MetaSVM_score']]) <  metasvm_cutoff:
            ben=ben+1
        else:
            dam=dam+1
    except ValueError:  
        var=var+1
    else:
        pass

    try:
        if float(cls[Funcanno_flgs['SIFT_score']]) >= sift_cutoff:
            ben=ben+1
        else:
            dam=dam+1
    except ValueError:  # the sift absent means many:  synonymous indel  stop, but synonymous also is no impact
        var=var+1
    else:
        pass


    if cls[Funcanno_flgs['Polyphen2_HDIV_pred']] == "P" or cls[Funcanno_flgs['Polyphen2_HDIV_pred']] == "D":
        dam=dam+1
    if cls[Funcanno_flgs['Polyphen2_HDIV_pred']] == "B" :
        ben=ben+1
    if cls[Funcanno_flgs['Polyphen2_HDIV_pred']] == "." :
        var=var+1

    if cls[Funcanno_flgs['MetaLR_pred']] == "D":
        dam=dam+1
    if cls[Funcanno_flgs['MetaLR_pred']] == "T" :
        ben=ben+1
    if cls[Funcanno_flgs['MetaLR_pred']] == "." :
        var=var+1


    if cls[Funcanno_flgs['FATHMM_pred']] == "D":
        dam=dam+1
    if cls[Funcanno_flgs['FATHMM_pred']] == "T" :
        ben=ben+1
    if cls[Funcanno_flgs['FATHMM_pred']] == "." :
        var=var+1

    if cls[Funcanno_flgs['MutationTaster_pred']] == "A" or cls[Funcanno_flgs['MutationTaster_pred']] == "D":
        dam=dam+1
    if cls[Funcanno_flgs['MutationTaster_pred']] == "P":
        ben=ben+1
    if cls[Funcanno_flgs['MutationTaster_pred']] == "." :
        var=var+1
    
    if cls[Funcanno_flgs['GERP++_RS']] == ".": 
        var=var+1
    else:
        if float(cls[Funcanno_flgs['GERP++_RS']])>= cutoff_conserv:
            dam=dam+1
        else:
            ben=ben+1
        




    if dam >4: PreP=2;
    if ben >3: PreP=0;
    if var >3: Prep=1;
    if dam==ben: PreP=1;

    #print "PreP=",PreP,dam,ben,var
    return(PreP)

def check_Path(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Pathway involvement
    2 Disease-associated pathways
    1 Involve disease-associated pathways or pathogenic pathways
    0 May or may not involve disease-associated pathways
    0 May or may not involve disease-associated pathways
    '''
    Path=0;
#   cancer_pathway_dic{}
#   cancers_gene_dict{}

    cls=line.split('\t')
    try:
        if cancer_pathway_dict[ cls[Funcanno_flgs['Gene']] ] != '' :
            Path=1
    except KeyError:
        pass
    else:
        pass
    try:
        if cancers_gene_dict[ cls[Funcanno_flgs['Gene']] ] == '1' :
            Path=1
    except KeyError:
        pass
    else:
        pass
    #print "Path=",Path
    return(Path)
    



def check_Pubs(line,Funcanno_flgs,Allels_flgs,lof_genes_dict):
    '''Publications: functional study, population study, other
    Therapeutic/Diagnostic/Prognostic: reported evidence with consensus
    Therapeutic: evidence of using FDA-approved therapies for different tumor types; phase 2 or 3 clinical trials for investigational therapies; Diagnostic/Prognostic: multiple lines of reported evidence without consensus
    None or no convincing evidence to determine clinical/biological significance
    Reported evidence supportive of benign/likely benign; or none
    '''
    Pubs=0;

    level=0
    global cancervar_d;
    cls=line.split('\t')
    clstt=cls[Funcanno_flgs['Otherinfo']].split(';')
    cancer_type="Cancer"
    #if (len(clstt[0])>0):
    #    cancer_type=clstt[0]
    try:
        cancer_type=paras['cancer_type']
    except KeyError:
        if (len(clstt[0])>0): cancer_type=clstt[0]
    else:
        pass
 

    gene_tr=cls[Funcanno_flgs['Gene']]
    func=cls[Funcanno_flgs['Func.refGene']]
    exonfunc=cls[Funcanno_flgs['Func.refGene']]
    # ABL1:NM_005157:exon6:c.C944T:p.T315I,ABL1:NM_007313:exon6:c.C1001T:p.T334I
    # AAChange.refGene  only check missense with AAchange
    line_tmp=cls[Funcanno_flgs['Gene']]+" "+cls[Funcanno_flgs['Func.refGene']]+" "+cls[Funcanno_flgs['ExonicFunc.refGene']]+" "+cancer_type
    line_tmp2=cls[Funcanno_flgs['AAChange.refGene']]
    #line_tmp2_sp=line_tmp2.split(',')
    line_tmp2_sp=re.split(';|,',line_tmp2)
    exon_func=cls[Funcanno_flgs['ExonicFunc.refGene']]
    marker_key0=gene_tr+"_"+"mut_any"
    for cls0 in line_tmp2_sp:
        cls0_1=cls0.split(':')
        if(len(cls0)>1 and len(line_tmp2_sp)>0 and len(cls0_1)==5 ):
            cls0_1=cls0.split(':')
            gene=cls0_1[0]
            transc=cls0_1[1]
            exon=re.sub('exon','',cls0_1[2])
            cdna_change=re.sub('c.','',cls0_1[3])
            amino_change=re.sub('p.','',cls0_1[4])
            ltt=len(amino_change)
            codon_num=amino_change[1:ltt-1]
            #print gene, transc,exon,cdna_change,amino_change,cancer_type 
            marker_key=gene+"_"+amino_change
            marker_key0=gene+"_"+"mut_any"
            marker_key1=gene+"_"+"exon_"+exon+"_any"
            marker_key2=gene+"_"+"codon_"+codon_num+"_any"
            marker_key00=""       
            marker_key11=""       
            marker_key22=""       
            if(cls[Funcanno_flgs['ExonicFunc.refGene']]=="nonsynonymous SNV"):
                marker_key00=gene+"_"+"missense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_missense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_missense"
                
            if exon_func.find("frameshift")>=0 and exon_func.find("nonframe")<0 : 
                marker_key00=gene+"_"+"frameshift_any"
                marker_key11=gene+"_"+"exon_"+exon+"_frameshift"
                marker_key22=gene+"_"+"codon_"+codon_num+"_frameshift"

            if exon_func.find("stopgain")>=0 or exon_func.find("stoploss")>=0:  
                marker_key00=gene+"_"+"nonsense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_nonsense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_nonsense"

            if exon_func.find("deletion")>=0 :  
                marker_key00=gene+"_"+"deletion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_deletion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_deletion"
            if exon_func.find("insertion")>=0 :  
                marker_key00=gene+"_"+"insertion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_insertion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_insertion"
 
            #deletion frameshift indel missense nonsense_any
            # nonsynonymous SNV;stopgain/stoploss(nonsense);frameshift insertion/deletion/substitution; nonframeshift insertion/deletion/substitution;
            default_s=""
            add_list=cancervar_markers_dict.get(marker_key,default_s)+cancervar_markers_dict.get(marker_key0,default_s)+cancervar_markers_dict.get(marker_key00,default_s)+cancervar_markers_dict.get(marker_key1,default_s)+cancervar_markers_dict.get(marker_key11,default_s)+cancervar_markers_dict.get(marker_key2,default_s)+cancervar_markers_dict.get(marker_key22,default_s)
            cgi_list=add_list
            #print gene,add_list
            level=0 # ABCD
            for i in cgi_list.split(","):
                #print i
                if(len(i)>0):
                    pos=int(i)
                    if(cancervar_d[pos][7]!=";" and cancervar_d[pos][7]!=""):
                        level=2
                        if( not  re.findall(cancer_type, cancervar_d[pos][8], flags=re.IGNORECASE)   ):
                            level=1;  # did not find the specific type, decrease level
 
            #print gene,level

    Pubs=level
    return(Pubs)




def assign(BP,line,Freqs_flgs,Funcanno_flgs,Allels_flgs):

    CBP=[0,0,0,0,0,0,0,0,0,0,0,0]

    Therapeutic=check_Thera(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[0]=Therapeutic

    Diagnosis=check_Diagno(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[1]=Diagnosis

    Prognosis=check_Progno(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[2]=Prognosis

    Mutation_type=check_Mut(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[3]=Mutation_type

    Variant_freq=check_VF(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[4]=Variant_freq

    Potential_germ=check_PotG(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[5]=Potential_germ

    Population_data=check_PopD(line,Freqs_flgs,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[6]=Population_data

    Germline_data=check_GermD(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[7]=Germline_data

    Somatic_data=check_SomD(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[8]=Somatic_data
    
    Predict_pathoge=check_PreP(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[9]=Predict_pathoge


    Pathway_invol=check_Path(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[10]=Pathway_invol

    Publications=check_Pubs(line,Funcanno_flgs,Allels_flgs,lof_genes_dict)
    CBP[11]=Publications

    #print CBP

    cls=line.split('\t')

    #begin process the user's evidence file
    if os.path.isfile(paras['evidence_file']):
        keys=cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
        keys=re.sub("[Cc][Hh][Rr]","",keys)
        try:
            evds=user_evidence_dict[keys] #CBP1=1;CBP2=1;Cancer=LUAD;
            for evd in evds.split(';'):
                evd_t=evd.split('=')
                if(len(evd_t)>1 and (not re.findall('grade', evd_t[0], flags=re.IGNORECASE)) ):
                    if int(evd_t[1])<=1:
                        #print ("%s %s %s " %(keys,evd_t[1],evd_t[0]))
                        if(evd_t[0].find('CBP')!=-1):
                            t=evd_t[0].find('CBP');
                            tt=evd_t[0];
                            tt3=int(tt[t+3:t+4])
                            if(t<len(evd_t[0])-2 and tt3<=12 ): CBP[tt3-1]=int(evd_t[1])
        except KeyError:
            pass
        else:
            pass

    # end process the user's evidence file





    cls=line.split('\t')
    if len(cls)>1:#esp6500siv2_all 1000g2015aug_all ExAC_ALL
        BP_out=classfy(CBP,Allels_flgs,cls)
        line_t="%s EVS=%s" %(BP_out,CBP)

        #print("%s " % BP_out)
        BP_out=line_t
        pass
    #BP=BP_out
    return(BP_out)


def search_key_index(line,dict):
    cls=line.split('\t')
    for key in dict.keys():
        for i in range(1,len(cls)):
            ii=i-1
            if key==cls[ii]:
                dict[key]=ii
                break
    return

def my_inter_var_can(annovar_outfile):
    newoutfile=annovar_outfile+".grl_p"
    newoutfile2=annovar_outfile+".cancervar"

    Freqs_flgs={'1000g2015aug_all':0,'esp6500siv2_all':0,'ExAC_ALL':0,'ExAC_AFR':0,'ExAC_AMR':0,'ExAC_EAS':0,'ExAC_FIN':0,'ExAC_NFE':0,'ExAC_OTH':0,'ExAC_SAS':0,'gnomAD_genome_ALL':0,'gnomAD_genome_AFR':0,'gnomAD_genome_AMR':0,'gnomAD_genome_EAS':0,'gnomAD_genome_FIN':0,'gnomAD_genome_NFE':0,'gnomAD_genome_OTH':0,'gnomAD_genome_ASJ':0}
    Funcanno_flgs={'Func.refGene':0,'ExonicFunc.refGene':0,'AAChange.refGene':0,'Gene':0,'Gene damage prediction (all disease-causing genes)':0,'CLNDBN':0,'CLNACC':0,'CLNDSDB':0,'dbscSNV_ADA_SCORE':0,'dbscSNV_RF_SCORE':0,'GERP++_RS':0,'LoFtool_percentile':0,'Interpro_domain':0,'rmsk':0,'SIFT_score':0,'phastCons20way_mammalian':0,'Gene.ensGene':0,'CLNSIG':0,'CADD_raw':0,'CADD_phred':0,'avsnp147':0,'AAChange.ensGene':0,'AAChange.knownGene':0,'MetaSVM_score':0,'cosmic70':0,'ICGC_Id':0,'ICGC_Occurrence':0,'Otherinfo':0,'Polyphen2_HDIV_pred':0,'MetaLR_pred':0,'MutationTaster_pred':0,'FATHMM_pred':0,'Otherinfo':0}
    Allels_flgs={'Chr':0,'Start':0,'End':0,'Ref':0,'Alt':0}
# ExAC_ALL esp6500siv2_all   1000g2015aug_all  SIFT_score    CADD_raw    CADD_phred  GERP++_RS   phastCons20way_mammalian  dbscSNV_ADA_SCORE   dbscSNV_RF_SCORE   Interpro_domain

    try:
        fh=open(newoutfile, "r")
        fw=open(newoutfile2, "w")
        strs=fh.read()
        line_sum=0;
        print("Notice: Begin the variants interpretation by CancerVar ")
        if re.findall('true',paras['otherinfo'], flags=re.IGNORECASE)  :
            fw.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \t CancerVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Chr","Start","End","Ref","Alt","Ref.Gene","Func.refGene","ExonicFunc.refGene", "Gene.ensGene","avsnp147","AAChange.ensGene","AAChange.refGene","Clinvar","CancerVar and Evidence","Freq_ExAC_ALL", "Freq_esp6500siv2_all","Freq_1000g2015aug_all", "Freq_gnomAD_genome_ALL","CADD_raw","CADD_phred","SIFT_score","GERP++_RS","phastCons20way_mammalian","dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "Interpro_domain","AAChange.knownGene","MetaSVM_score","Freq_ExAC_POPs","OMIM","Phenotype_MIM","OrphaNumber","Orpha","Otherinfo"  ))
        else:
            fw.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \t CancerVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Chr","Start","End","Ref","Alt","Ref.Gene","Func.refGene","ExonicFunc.refGene", "Gene.ensGene","avsnp147","AAChange.ensGene","AAChange.refGene","Clinvar","CancerVar and Evidence","Freq_ExAC_ALL", "Freq_esp6500siv2_all","Freq_1000g2015aug_all", "Freq_gnomAD_genome_ALL","CADD_raw","CADD_phred","SIFT_score","GERP++_RS","phastCons20way_mammalian","dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "Interpro_domain","AAChange.knownGene","MetaSVM_score","Freq_ExAC_POPs","OMIM","Phenotype_MIM","OrphaNumber","Orpha"  ))
        for line in strs.split('\n'):
            BP="UNK" # the inter of pathogenetic/benign
            clinvar_bp="UNK"
            cls=line.split('\t')
            if len(cls)<2: break
            if line_sum==0:
                search_key_index(line,Freqs_flgs)
                search_key_index(line,Funcanno_flgs)
                search_key_index(line,Allels_flgs)

            else:
                #begin check the BP status from clinvar
                line_tmp2=cls[Funcanno_flgs['CLNSIG']]
                if line_tmp2 != '.':
                    cls3=line_tmp2.split(';')
                    clinvar_bp=cls3[0]

                cancervar_bp=assign(BP,line,Freqs_flgs,Funcanno_flgs,Allels_flgs)
                Freq_ExAC_POPs="AFR:"+cls[Freqs_flgs['ExAC_AFR']]+",AMR:"+cls[Freqs_flgs['ExAC_AMR']]+",EAS:"+cls[Freqs_flgs['ExAC_EAS']]+",FIN:"+cls[Freqs_flgs['ExAC_FIN']]+",NFE:"+cls[Freqs_flgs['ExAC_NFE']]+",OTH:"+cls[Freqs_flgs['ExAC_OTH']]+",SAS:"+cls[Freqs_flgs['ExAC_SAS']]
                OMIM="."
                mim2=mim2gene_dict2.get(cls[Funcanno_flgs['Gene']],".")
                mim1=mim2gene_dict.get(cls[Funcanno_flgs['Gene.ensGene']],".")
                if(mim1 !="."):
                   OMIM=mim1
                if(mim2 !="."):
                   OMIM=mim2
                Pheno_MIM=mim_pheno_dict.get(OMIM,".")
                orpha="";
                orpha_details="";
                # .;442835;;306;;.;
                for ort2 in Pheno_MIM.split(';'):
                    ort3=mim_orpha_dict.get(ort2,".")
                    if(ort3 !="."):
                        orpha=ort3+orpha
                for ort4 in orpha.split(';'):
                    if len(ort4)>0:
                         orpha_details=orpha_details+orpha_dict.get(ort4,".")+"~"
                if(orpha ==""):
                    orpha="."
                if(orpha_details ==""):
                    orpha_details="."


                if re.findall('true',paras['otherinfo'], flags=re.IGNORECASE)  :
                    fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \t CancerVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cls[Allels_flgs['Chr']],cls[Allels_flgs['Start']],cls[Allels_flgs['End']],cls[Allels_flgs['Ref']],cls[Allels_flgs['Alt']],cls[Funcanno_flgs['Gene']],cls[Funcanno_flgs['Func.refGene']],cls[Funcanno_flgs['ExonicFunc.refGene']], cls[Funcanno_flgs['Gene.ensGene']],cls[Funcanno_flgs['avsnp147']],cls[Funcanno_flgs['AAChange.ensGene']],cls[Funcanno_flgs['AAChange.refGene']],clinvar_bp,cancervar_bp,cls[Freqs_flgs['ExAC_ALL']], cls[Freqs_flgs['esp6500siv2_all']], cls[Freqs_flgs['1000g2015aug_all']],cls[Freqs_flgs['gnomAD_genome_ALL']], cls[Funcanno_flgs['CADD_raw']],cls[Funcanno_flgs['CADD_phred']],cls[Funcanno_flgs['SIFT_score']],  cls[Funcanno_flgs['GERP++_RS']],cls[Funcanno_flgs['phastCons20way_mammalian']], cls[Funcanno_flgs['dbscSNV_ADA_SCORE']], cls[Funcanno_flgs['dbscSNV_RF_SCORE']], cls[Funcanno_flgs['Interpro_domain']],cls[Funcanno_flgs['AAChange.knownGene']],cls[Funcanno_flgs['MetaSVM_score']],Freq_ExAC_POPs,OMIM,Pheno_MIM,orpha,orpha_details,cls[Funcanno_flgs['Otherinfo']]     ))
                else:
                    fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tclinvar: %s \t CancerVar: %s \t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cls[Allels_flgs['Chr']],cls[Allels_flgs['Start']],cls[Allels_flgs['End']],cls[Allels_flgs['Ref']],cls[Allels_flgs['Alt']],cls[Funcanno_flgs['Gene']],cls[Funcanno_flgs['Func.refGene']],cls[Funcanno_flgs['ExonicFunc.refGene']], cls[Funcanno_flgs['Gene.ensGene']],cls[Funcanno_flgs['avsnp147']],cls[Funcanno_flgs['AAChange.ensGene']],cls[Funcanno_flgs['AAChange.refGene']],clinvar_bp,cancervar_bp,cls[Freqs_flgs['ExAC_ALL']], cls[Freqs_flgs['esp6500siv2_all']], cls[Freqs_flgs['1000g2015aug_all']], cls[Freqs_flgs['gnomAD_genome_ALL']],cls[Funcanno_flgs['CADD_raw']],cls[Funcanno_flgs['CADD_phred']],cls[Funcanno_flgs['SIFT_score']],  cls[Funcanno_flgs['GERP++_RS']],cls[Funcanno_flgs['phastCons20way_mammalian']], cls[Funcanno_flgs['dbscSNV_ADA_SCORE']], cls[Funcanno_flgs['dbscSNV_RF_SCORE']], cls[Funcanno_flgs['Interpro_domain']],cls[Funcanno_flgs['AAChange.knownGene']],cls[Funcanno_flgs['MetaSVM_score']],Freq_ExAC_POPs,OMIM,Pheno_MIM,orpha,orpha_details   ))
                #print("%s\t%s %s" % (line,clinvar_bp,cancervar_bp))

            line_sum=line_sum+1

    except IOError:
        print("Error: can\'t readi/write the annovar output files %s" % (newoutfile,newoutfile2))
        sys.exit()
        return
    else:
        fh.close()
        fw.close()
    return(line_sum)


def main():


    if platform.python_version()< '3.0.0'  :
        config=ConfigParser.ConfigParser()
    else:
        config=configparser.ConfigParser()





    parser = optparse.OptionParser(usage=usage, version=version, description=description)


    parser.add_option("-?", action="help", help=optparse.SUPPRESS_HELP, dest="help")
    parser.add_option("-v", action="version", help=optparse.SUPPRESS_HELP, dest="version")

    parser.add_option("-c", "--config", dest="config", action="store",
                  help="The config file of all options. it is for your own configure file.You can edit all the options in the configure and if you use this options,you can ignore all the other options bellow", metavar="config.ini")

    parser.add_option("-b", "--buildver", dest="buildver", action="store",
                  help="The genomic build version, it can be hg19 and will support GRCh37 hg18 GRCh38 later", metavar="hg19")


    parser.add_option("-i", "--input", dest="input", action="store",
                  help="The input file contains your variants", metavar="example/ex1.avinput")

    parser.add_option("--input_type", dest="input_type", action="store",
                  help="The input file type, it can be  AVinput(Annovar's format),VCF(VCF with single sample),VCF_m(VCF with multiple samples)", metavar="AVinput")

    parser.add_option("--cancer_type", dest="cancer_type", action="store",
                  help="The cancer type, please check the doc for the details of cancer types: Adrenal_Gland Bile_Duct Bladder Blood Bone Bone_Marrow Brain Breast Cancer Cancervar_type Cervix Colorectal Esophagus Eye Head_and_Neck Inflammatory Intrahepatic Kidney Liver Lung Lymph_Nodes Nervous_System Other Ovary Pancreas Pleura Prostate Skin Soft_Tissue Stomach Testis Thymus Thyroid Uterus)", metavar="CANCER")

    parser.add_option("-o", "--output", dest="output", action="store",
                  help="The prefix of output file which contains the results, the file of results will be as [$$prefix].cancervar ", metavar="example/myanno")


    group = optparse.OptionGroup(parser, "CancerVar Other Options")
    group.add_option("-t", "--database_cancervar", dest="database_cancervar", action="store",
            help="The  database location/dir for the CancerVar dataset files", metavar="cancervardb")
    group.add_option("-s", "--evidence_file", dest="evidence_file", action="store",
            help="User specified Evidence file for each variant", metavar="your_evidence_file")
    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "   How to add your own Evidence for each Variant",
    """ Prepare your own evidence  file as tab-delimited,the line format:
         The format for upgrad/downgrade of criteria should be like: grade_CBPx=2;
         3 for Strong; 2 for Moderate; 1 for Supporting)
            Chr Pos_start Pos_end Ref Alt CBP2=1;grade_CBP2=2;CBP9=1
                                """)
    parser.add_option_group(group)



    group = optparse.OptionGroup(parser, "Annovar Options",
                                "Caution: check these options from manual of Annovar. The ANNOVAR version should be >=  2016-02-01, older verions of ANNOVAR will bring problems.")
    group.add_option("--table_annovar", action="store", help="The Annovar perl script of table_annovar.pl",metavar="./table_annovar.pl",dest="table_annovar")
    group.add_option("--convert2annovar", action="store", help="The Annovar perl script of convert2annovar.pl",metavar="./convert2annovar.pl",dest="convert2annovar")
    group.add_option("--annotate_variation", action="store", help="The Annovar perl script of annotate_variation.pl",metavar="./annotate_variation.pl",dest="annotate_variation")
    group.add_option("-d", "--database_locat", dest="database_locat", action="store",
            help="The  database location/dir for the annotation datasets", metavar="humandb")
    group.add_option("--skip_annovar", action="store_true", help="Skip the Annovar annotation, this can be true only after you  already got the annovar's annotation results",dest="skip_annovar")

    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "Examples",
    """./CancerVar.py -c config.ini  # Run the examples in config.ini
    ./CancerVar.py -i your_input  --input_type=VCF  -o your_output
                                """)
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    #(options,args) = parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    print("%s" %description)
    print("%s" %version)
    print("Notice: Your command of CancerVar is %s" % sys.argv[:])




    if os.path.isfile("config.ini"):
        config.read("config.ini")
        sections = config.sections()
        for section in sections:
            ConfigSectionMap(config,section)
    else:
        print("Error: The default configure file of [ config.ini ] is not exit! Please redownload the CancerVar.")
        sys.exit()

#begin to process user's options:
    if options.config != None:
        if os.path.isfile(options.config):
            config.read(options.config)
            sections = config.sections()
            for section in sections:
                ConfigSectionMap(config,section)
        else:
            print("Error: The config file [ %s ] is not here,please check the path of your config file." % options.config)
            sys.exit()

    if options.buildver != None:
        paras['buildver']=options.buildver
    if options.database_locat != None:
        paras['database_locat']=options.database_locat
    if options.input != None:
        paras['inputfile']=options.input
    if options.input_type != None:
        paras['inputfile_type']=options.input_type
    if options.output != None:
        paras['outfile']=options.output
    if options.evidence_file != None:
        paras['evidence_file']=options.evidence_file
        print("Warning: You provided your own evidence file [ %s ] for the CancerVar." % options.evidence_file)
    if options.cancer_type != None:
        paras['cancer_type']=options.cancer_type
    if options.database_cancervar != None:
        paras['database_cancervar']=options.database_cancervar
        paras['lof_genes'] = paras['database_cancervar']+'/PVS1.LOF.genes'
        paras['mim2gene'] =paras['database_cancervar']+'/mim2gene.txt'
        paras['mim_pheno'] = paras['database_cancervar']+'/mim_pheno.txt'
        paras['mim_orpha'] = paras['database_cancervar']+'/mim_orpha.txt'
        paras['orpha'] = paras['database_cancervar']+'/orpha.txt'
        paras['knowngenecanonical'] = paras['database_cancervar']+'/knownGeneCanonical.txt'
        paras['exclude_snps'] = paras['database_cancervar']+'/ext.variants'
        paras['cancer_pathway'] =paras['database_cancervar']+'/cancer_pathway.list'
        paras['cancers_genes'] =paras['database_cancervar']+'/cancers_genes.list'
        paras['cancers_types']=paras['database_cancervar']+'/cancervar.cancer.types'


    #paras['ps1_aa'] = paras['ps1_aa']+'.'+paras['buildver']
    #paras['ps4_snps'] = paras['ps4_snps']+'.'+paras['buildver']
    #paras['bs2_snps'] = paras['bs2_snps']+'.'+paras['buildver']
    paras['exclude_snps'] = paras['exclude_snps']+'.'+paras['buildver']

    if options.table_annovar != None:
        if os.path.isfile(options.table_annovar):
            paras['table_annovar']=options.table_annovar
        else:
            print("Error: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                    % options.table_annovar)
            sys.exit()
    if options.convert2annovar != None:
        if os.path.isfile(options.convert2annovar):
            paras['convert2annovar']=options.convert2annovar
        else:
            print("Error: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                    % options.convert2annovar)
            sys.exit()
    if options.annotate_variation != None:
        if os.path.isfile(options.annotate_variation):
            paras['annotate_variation']=options.annotate_variation
        else:
            print("Error: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                    % options.annotate_variation)
            sys.exit()


    if not os.path.isfile(paras['inputfile']):
        print("Error: Your input file [ %s ] is not here,please check the path of your input file." % paras['inputfile'])
        sys.exit()
    if  not os.path.isfile(paras['evidence_file']) and paras['evidence_file']!="None":
        print("Warning: Your specified evidence file [ %s ] is not here,please check the path of your evidence file." % paras['evidence_file'])
        print("         Your analysis will begin without your specified evidence.")
    else: 
        print("Warning: Your specified evidence file [ %s ], the analysis will take your additional evidence." % paras['evidence_file'])





    print ("INFO: The options are %s " % paras)
    try:
        print("Warning: Your specified the cancer type: %s in command, the cancer type in otherinfo column will be replaced for all your variants!!!" % paras['cancer_type'])
    except KeyError:
        pass
    else:
        pass

    check_downdb()
    check_input()
    if  options.skip_annovar != True   :
        check_annovar_result() #  to obtain myanno.hg19_multianno.csv
    else:
         print ("Warning: You activated the option of --skip_annovar, the Annovar will not run!")
         print ("Warning: The InterVar will interpret the variants based on your old annotation information!")

    #check_annovar_result() #  to obtain myanno.hg19_multianno.csv
    read_datasets()
    inputft= paras['inputfile_type']
    some_file_fail=0
    out_annf=0;
    for annovar_outfile  in glob.iglob(paras['outfile']+"*."+paras['buildver']+"_multianno.txt"):
        print("annovar_outfile is %s" % annovar_outfile)
        sum1=check_genes(annovar_outfile)
        sum2=my_inter_var_can(annovar_outfile)
        out_annf=out_annf+1;

        outfile=annovar_outfile+".cancervar"
        if os.path.isfile(outfile):
            print ("Notice: About %d lines in your variant file! " % (sum1-1))
            print ("Notice: About %d variants has been processed by CancerVar" % (sum2-1))
            if inputft.lower() != 'vcf_m' :
                print ("Notice: The CancerVar is finished, the output file is [ %s.cancervar ]" % annovar_outfile)
        else:
            some_file_fail=some_file_fail+1
            print ("Warning: The CancerVar seems not run correctly, please check your inputs and options in configure file")

    if inputft.lower() == 'vcf_m' :
        print ("Notice: The CancerVar for VCF with multiple samples is finished, the output file is as [ %s.<samplename>.cancervar ]" % annovar_outfile)
        sum_sample=1;
        for f in glob.iglob(paras['outfile']+"*."+paras['buildver']+"_multianno.txt.cancervar"):
            print ("Notice: The CancerVar for VCF with multiple samples is finished, The %d sample output file is [ %s]" %(sum_sample,f))
            sum_sample=sum_sample+1;
        if some_file_fail>=1:
            print ("Warning: The CancerVar seems not run correctly for your %d samples in the VCF, please check your inputs and options in configure file" %  some_file_fail )
    if out_annf==0:
         print ("Warning: The InterVar seems not run correctly, please check your inputs , options and configure file!")
         print ("ERROR: The InterVar did not find the annotation result file from ANNOVAR!")
         print ("ERROR: The name of annotation result file should be like %s*.%s__multianno.txt" % (paras['outfile'],paras['buildver']))
    print("%s" %end)






if __name__ == "__main__":
    main()

