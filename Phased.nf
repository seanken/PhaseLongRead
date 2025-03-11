#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.snpVCF //the vcf with SNP information
params.svVCF //the vcf from sniffles
params.bam //the bam of long reads mapped to the genome
params.vcf_col //column with sample information
params.vcf_col_sv=params.vcf_col //column with sample information in the sv vcf if different than the name for the SNP one
phaseReadReadsScript="$projectDir/scripts/AlleleMASSeq.jar" //the script that phases reads
combineVCFScript="$projectDir/scripts/combineVCF.py" //combines the allele specific VCFs into one joint, phased vcf
params.outdir="output"

workflow{
    newvcf=PrepVCF(params.vcf_col,params.snpVCF)
    PhasedReads(params.bam,newvcf,phaseReadReadsScript)
    PrepSV_VCF(params.vcf_col_sv,params.svVCF)
    CallAlleleOneSV(PhasedReads.out.bam1,PrepSV_VCF.out)
    CallAlleleTwoSV(PhasedReads.out.bam2,PrepSV_VCF.out)
    CombineVCF(CallAlleleOneSV.out,CallAlleleTwoSV.out,combineVCFScript)
}



process PrepVCF
{

    input:
    env vcf_col 
    path "input.vcf.gz" 

    output:
    path "new.vcf" 


    '''
    tabix -p vcf input.vcf.gz
    bcftools view -H -O v -s $vcf_col input.vcf.gz | grep -v "0|0" | grep -v "1|1" | grep -v "\\.|\\."  > new.vcf
    '''


}


process PrepSV_VCF
{

    input:
    env vcf_col 
    path "input.vcf" 

    output:
    path "new.sv.vcf" 


    '''
    bcftools sort input.vcf > input.sort.vcf
    bgzip -c input.sort.vcf > input.sort.vcf.gz
    tabix -p vcf input.sort.vcf.gz
    bcftools view -H -O v -s $vcf_col input.sort.vcf.gz | grep -v "0/0" | grep -v "1/1" | grep -v "\\./\\."  > new.sv.vcf
    '''


}

process PhasedReads
{
    input:
    path "input.bam"
    path "snps.vcf"
    path "phaseReads.jar"

    output:
    path "allele1.bam", emit: bam1
    path "allele2.bam", emit: bam2
    path "failed.bam", emit: bam3

    '''
    java -jar phaseReads.jar input.bam snps.vcf allele1.bam allele2.bam failed.bam
    '''
}


process CallAlleleOneSV
{
    input:
    path "allele1.bam"
    path "sv.vcf"

    output:
    path "allele1_sv.vcf"

    '''
    samtools view -b --expr '[AL]>3 || [AM]>3' allele1.bam > allele1.top.bam
    samtools index allele1.top.bam
    sniffles --input allele1.top.bam --genotype-vcf sv.vcf --vcf allele1_sv.vcf
    '''
}


process CallAlleleTwoSV
{
    input:
    path "allele2.bam"
    path "sv.vcf"

    output:
    path "allele2_sv.vcf"

    '''
    samtools view -b --expr '[AL]>3 || [AM]>3' allele2.bam > allele2.top.bam
    samtools index allele2.top.bam
    sniffles --input allele2.top.bam --genotype-vcf sv.vcf --vcf allele2_sv.vcf
    '''
}

process CombineVCF
{
    publishDir "${params.outdir}/PhasedVCF", mode: 'copy'

    input:
    path "allele1.vcf"
    path "allele2.vcf"
    path "combineVCF.py"

    output:
    path "comb.vcf"

    '''
    python combineVCF.py allele1.vcf allele2.vcf comb.vcf
    '''
}