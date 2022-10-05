import os
import numpy
from multiprocessing import Pool
path_bed='/share/inspurStorage/home1/zhaoxz/ldsc/baseline/1000G_EUR_Phase3_plink'
snps='/share/inspurStorage/home1/zhaoxz/ldsc/list.txt'
outdir='/home1/zhaoxz/ldsc/baseline/annot/allanote/allannote/'
def make_ano(promoter_index):
    path_file=promoter_index.split('.')[0]
   # file_prefix=os.path.join(outdir,path_file)
    promoter_files=os.path.join(exp_files,promoter_index)
    for j in range(1,23):
        tem_file='1000G.EUR.QC.'+str(j)+'.bim'
        annofile=os.path.join(path_bed,tem_file)
        os.system('python /home1/zhaoxz/ldsc/ldsc-master/make_annot.py --bed-file '+
                  promoter_files+' --bimfile '+annofile+' --annot '+outdir+'/'+path_file+'.'+str(j)+
                  '.annot.gz')
def partion_ldsc_cal(annotfile):
    annotfile=annotfile.split('.')[0]
    file_prex=os.path.join(outdir,annotfile)
    for j in range(1,23):
        tem_file = '1000G.EUR.QC.' + str(j)
        tem_anno=file_prex+'.'+str(j)+'.annot.gz'
        annot_d=os.path.join(outdir,tem_anno)
        bed_using=os.path.join(path_bed,tem_file)
        os.system('python /home1/zhaoxz/ldsc/ldsc-master/ldsc.py --l2'+
		' --bfile '+bed_using+
		' --ld-wind-cm 1'+
		' --annot ' +annot_d+
		' --thin-annot '+
		'--out ' +outdir+annotfile+'.'+str(j)+
		' --print-snps ' +snps)
def combined_fc(data):
    make_ano(data)
    partion_ldsc_cal(data)
if __name__ == "__main__":
    exp_files='EPI_folder'
    enhance_file=os.listdir(exp_files)
    kernal=20
    pool = Pool(processes = kernal)
    for i in range(len(enhance_file)):
        sub=enhance_file[i]
        #combined_fc(sub)
        pool.apply_async(combined_fc,(sub,))
    pool.close()
    pool.join()



