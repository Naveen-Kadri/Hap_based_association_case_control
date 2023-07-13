'''
Input files :
  genotypes in phased .vcf
  phenotypes in a .txt file with two columns : id & phenotype (0,1)
Test :
Fisher's exact test
Naveen.Kadri@usys.ethz.ch Nov 29 2021
'''

import gzip
from collections import defaultdict
from scipy.stats import fisher_exact as fe
        
##FILES
vcf_file="/cluster/work/pausch/arnav/imputed_pigs/1.vcf.gz"
out_file = "test.txt"
pheno_file="/cluster/work/pausch/arnav/HIS/gwas/haplotype_gwas/pheno.txt"


##PARAMETERS
inheritance="additive"
window_size=10
overlap=3
freq_thresh=0.10

cases, controls =dict (), dict()
phenotyped=dict()
hap1, hap2=list (), list ()
pos=list ()
nvar = -1
nwin = 0
ngenopheno=0
ncases=0
ncontrols=0


def FishersTest ():
    hap_freq =dict ()
    genotypes = defaultdict (dict)
    for myid, myhap1, myhap2 in zip (ids, hap1, hap2):
        if myid in phenotyped:
            hap_freq [myhap1] = hap_freq.get (myhap1,0)+1
            hap_freq [myhap2] = hap_freq.get (myhap2,0)+1
            genotypes [myhap1] [myid] = genotypes [myhap1].get (myid,0) +1
            genotypes [myhap2] [myid] = genotypes [myhap2].get (myid,0) +1
    hap_freq = {myhap : hap_freq [myhap]/ nhap for myhap in hap_freq}
    ##check for frequencies
    #print (sum (list (hap_freq.values ())))
    for myhap, myfreq in hap_freq.items ():
        if myfreq < freq_thresh or myfreq > (1-freq_thresh):
            continue
        incase = 0
        incontrol=0
        for myid in genotypes [myhap]:
            mygenotypes = genotypes [myhap]


            if inheritance == "additive":
                if myid in cases and  mygenotypes [myid] >0:
                    incase+=1
                if myid in controls and mygenotypes [myid] >0:
                    incontrol+=1
            else:
                if myid in cases and  mygenotypes [myid] == 2:
                    incase+=1
                if myid in controls and mygenotypes [myid] !=2:
                    incontrol+=1
                    
        mytable=[[incase,ncases-incase], [incontrol,ncontrols-incontrol]]
        myor, myp=fe(mytable, alternative='two-sided')
        myp = f"{myp:.2e}"
        mytable = [sel for fel in mytable for sel in fel ]
        tw = "\t".join ( [str (el)  for el in [pos[0], pos[-1],myhap, round(myfreq,2)] +  mytable + [round(myor,2), myp ]] )
        print (nwin, pos[0], pos[1])
        out.write (f"{tw}\n")
                
out = open (out_file, "w")
header = "\t".join ("startpos endpos hap hapfrq case+ case- control+ control-  OR pvalue".split (" "))
out.write (f"{header}\n")

print ('reading the phenotype file')
with open (pheno_file) as inf:
    for line in inf:
        myid, mypheno=line.rstrip().split()
        phenotyped [myid] =1
        if mypheno=="1":
            cases[myid]=1
        else:
            controls [myid]=1
print (f"Number of phenotyped : {len (phenotyped)}")            
print (f"Number of cases in the phenotype file : {len (cases)  }")
print (f"Number of controls in the phenotype file : {len (controls)  }")


with gzip.open (vcf_file, "rt") as inf:
    for line in inf:
        if line [0:2] != "##":
            spl=line.rstrip ().split ()
            if line [0:6] == "#CHROM":
                ids = spl [9:]
                pids=[myid for myid in ids if myid in phenotyped]
                print (len (pids))
                nani = len (ids)
                ##how many cases among samples with genotypes
                for myid in ids:
                    if myid in cases:
                        ncases+=1
                    elif myid in controls:
                        ncontrols+=1
                ngenopheno = ncases + ncontrols
                nhap = ngenopheno *2
                if ngenopheno ==0:
                    print ("No overlap between geno and pheno files")
                    exit ()
                else:
                    print ("--------------------------------------------------------\n")
                    print (f"Number of samples with geno and pheno : {ngenopheno}")
                    print (f"Number of cases, controls among geno + pheno : {ncases}, {ncontrols}")
                    print ("\n")
                    print ("--------------------------------------------------------\n")
            else:
                if "," not in spl[3] and "," not in spl[4]:       
                    gts = spl [9:]
                    nvar+=1
                    pos.append (int(spl [1]))
                    if nvar==window_size:
                        nwin +=1
                        #print (pos )
                        FishersTest()
                        nvar = overlap  #already you have n=overlap snps !
                        if overlap ==0:
                            hap1, hap2, pos = [], [], []
                        else:
                            hap1 = [myhap [-overlap:]  for myhap in hap1]
                            hap2 = [myhap [-overlap:]  for myhap in hap2]
                            pos  = pos [-overlap:]

                    for i, gt in enumerate (gts):
                        if nvar ==0:
                            hap1.append(gt [0])
                            hap2.append(gt [2])
                        else:
                            hap1 [i]= hap1 [i] + gt [0] 
                            hap2 [i]= hap2 [i] + gt [2]
                            
##For the last window
FishersTest()
out.close ()
