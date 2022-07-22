import os
import numpy
from multiprocessing import Pool
class anao_process(obeject):
    """
     :param path_bed: bea files
     :param promoter_files: the area of promoter and enhancer data
     :param outdir: the directory of the ld socre and annotation files
     :param type_d: the file type[enhancer or promoter-presufix of generate files
     :param snp_ano: filter the SNP files
     :param exp_files: the gene region for annotion
     """
    def __init__(self,path_bed,promoter_files,outdir,type_d,snp_ano,exp_files):
        self.path_bed=path_bed
        self.promoter_files=promoter_files
        self.outdir=outdir
        self.type_d=type_d
        self.snp_ano=snp_ano
        self.exp_files=exp_files
    def make_ano(self):
        for j in range(1,23):
            tem_file='1000G.EUR.QC.'+str(j)+'.bim'
            annofile=os.path.join(self.path_bed,tem_file)
            os.system('python /home1/zhaoxz/ldsc/ldsc-master/make_annot.py --bed-file '+
                      self.promoter_files+' --bimfile '+annofile+' --annot '+self.outdir+'/'+self.type_d+'.'+str(j)+
                      '.annot.gz')
    def partion_ldsc_cal(self,annotfile,snps,anoto_type,outdir,file_prex):
        for j in range(1,23):
            tem_file = '1000G.EUR.QC.' + str(j)
            tem_anno=file_prex+'.'+str(j)+'.annot.gz'
            annot_d=os.path.join(self.outdir,tem_anno)
            bed_using=os.path.join(self.path_bed,tem_file)
            os.system('python /home1/zhaoxz/ldsc/ldsc-master/ldsc.py --l2'+
            ' --bfile '+bed_using+
            ' --ld-wind-cm 1'+
            ' --annot ' +annot_d+
            ' --thin-annot '+
            '--out ' +outdir+anoto_type+'.'+str(j)+
            ' --print-snps ' +self.snps)

if __name__ == "__main__":
    path_bed='/share/inspurStorage/home1/zhaoxz/ldsc/baseline/1000G_EUR_Phase3_plink'
    snp_ano='/share/inspurStorage/home1/zhaoxz/ldsc/list.txt'
    exp_files='/home1/zhaoxz/ldsc/baseline/annot/enhancer/'
    outdir='/home1/zhaoxz/ldsc/baseline/annot/allannot/'
    #annotfile='/home1/zhaoxz/ldsc/baseline/annot/promoter/'
    enhance_file=os.listdir(exp_files)
    kernal=20
    pool = Pool(processes = kernal)
    for i in enhance_file:
        path_file=i.split('.')[0]
        file_prefix=os.path.join(outdir,path_file)
        path_all=os.path.join(exp_files,i)
        #print(path_all)
        #print(file_prefix)
        make_ano(path_bed,path_all,outdir,path_file)
        partion_ldsc_cal(path_bed,outdir,snp_ano,path_file,outdir,file_prefix)




