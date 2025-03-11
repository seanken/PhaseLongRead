import pandas 
import sys

##Loads VCF
def LoadVCF(vcf):
    dat=pandas.read_csv(vcf,sep="\t",comment="#",header=None)
    cols=["col"+str(i) for i in range(0,10)]
    dat.columns=cols
    return(dat)


def PhaseGeno(call1,call2):
    if call1=="0/0" and call2=="1/1":
        return("0/1");
    if call2=="0/0" and call1=="1/1":
        return("1/0");
    return("./.")

def CombineVCF(vcf1,vcf2):
    dat1=LoadVCF(vcf1)
    dat2=LoadVCF(vcf2)
    dat3=dat1.copy()
    geno1=[i.split(":")[0] for i in dat1["col9"]]
    geno2=[i.split(":")[0] for i in dat2["col9"]]
    combineGeno=[PhaseGeno(geno1[i],geno2[i]) for i in range(0,dat1.shape[0])]
    dat3["col9"]=combineGeno
    dat3=dat3[dat3["col9"]!="./."]
    return(dat3)

if __name__=="__main__":
    args=sys.argv
    vcf1=args[1]
    vcf2=args[2]
    vcf_out=args[3]
    dat=CombineVCF(vcf1,vcf2)
    dat.to_csv(vcf_out,header=False,sep="\t",index=False)