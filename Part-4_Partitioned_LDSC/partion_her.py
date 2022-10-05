import os
from multiprocessing import Pool
gwas_file = '17_disorders_and_traits_folder/'
outdir = 'Output_dir'
annot_path = os.listdir('/share/inspurStorage/home1/zhaoxz/ldsc/baseline/annot/EPI_annot/')
def get_h2_res(phone_type):
    phone_name=phone_type.split('.')[0]
    phone_type = os.path.join(gwas_file, phone_type)
    for i in annot_path:
        path_data=i.split('.')[0]
        path_data_en=path_data+'.'
        out_put_file_en=os.path.join(outdir,phone_name+'_'+path_data)
        #out_put_file_pro_en=os.path.join(outdir,phone_name+'_'+annotfile_promote_enhancer)
        os.system('python /share/inspurStorage/home1/zhaoxz/ldsc/ldsc-master/ldsc.py --h2 '+
         phone_type+' --w-ld-chr /home1/zhaoxz/ldsc/baseline/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. '+
        '--ref-ld-chr/share/inspurStorage/home1/zhaoxz/ldsc/baseline/annot/EPI_annot/'+
        path_data_en+',/home1/zhaoxz/ldsc/baseline/baselineLD_v2.1/baselineLD.'+
         ' --overlap-annot --frqfile-chr /home1/zhaoxz/ldsc/baseline/1000G_Phase3_frq/1000G.EUR.QC.'+
        ' --out ' +out_put_file_en  +' --print-coefficients')
if __name__ == '__main__':
    kernal=20
    pool = Pool(processes = kernal)
    gwas_file='17_disorders_and_traits_folder/'
    allfile=os.listdir(gwas_file)
    n=len(allfile)
    num=0
    for i in range(num,n):
        sub = allfile[i]
        pool.apply_async(get_h2_res,(sub,))
pool.close()
pool.join()