import os
from multiprocessing import Pool
def get_h2_res(phone_type):
  
    phone_name=phone_type.split('.')[0]
    print(phone_name)
    gwas_file='/home1/zhaoxz/zhao/brain_volume/volume_sumstats/'
    outdir='/home1/zhaoxz/ldsc/imgae_single_cell'
    phone_type=os.path.join(gwas_file,phone_type)
    annot_path=os.listdir('/share/inspurStorage/home1/zhaoxz/ldsc/baseline/annot/cell_type/annote/link_annot/')
    for i in annot_path:
        path_data=i.split('.')[0]
        path_data_en=path_data+'.'
        out_put_file_en=os.path.join(outdir,phone_name+'_'+path_data)
        #out_put_file_pro_en=os.path.join(outdir,phone_name+'_'+annotfile_promote_enhancer)
        os.system('python /share/inspurStorage/home1/zhaoxz/ldsc/ldsc-master/ldsc.py --h2 '+
         phone_type+' --w-ld-chr /home1/zhaoxz/ldsc/baseline/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. '+
        '--ref-ld-chr /share/inspurStorage/home1/zhaoxz/ldsc/baseline/annot/cell_type/annote/link_annot/'+
        path_data_en+',/home1/zhaoxz/ldsc/baseline/baselineLD_v2.1/baselineLD.'+
         ' --overlap-annot --frqfile-chr /home1/zhaoxz/ldsc/baseline/1000G_Phase3_frq/1000G.EUR.QC.'+
        ' --out ' +out_put_file_en  +' --print-coefficients')
        # os.system('python /share/inspurStorage/home1/zhaoxz/ldsc/ldsc-master/ldsc.py --h2 '+
        # phone_type+' --w-ld-chr /home1/zhaoxz/ldsc/baseline/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. '+
        # '--ref-ld-chr /home1/zhaoxz/ldsc/baseline/annot/allannot/'+
        # annotfile_promote+',/home1/zhaoxz/ldsc/baseline/baselineLD_v2.1/baselineLD.'+
        # ' --overlap-annot --frqfile-chr /home1/zhaoxz/ldsc/baseline/1000G_Phase3_frq/1000G.EUR.QC.'+
        # ' --out ' +out_put_file_pro  +' --print-coefficients')
        # os.system('python /share/inspurStorage/home1/zhaoxz/ldsc/ldsc-master/ldsc.py --h2 '+
        # phone_type+' --w-ld-chr /home1/zhaoxz/ldsc/baseline/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. '+
        # '--ref-ld-chr /home1/zhaoxz/ldsc/baseline/annot/allannot/'+
        # annotfile_promote_enhancer+',/home1/zhaoxz/ldsc/baseline/baselineLD_v2.1/baselineLD.'+
        # ' --overlap-annot --frqfile-chr /home1/zhaoxz/ldsc/baseline/1000G_Phase3_frq/1000G.EUR.QC.'+
        # ' --out ' +out_put_file_pro_en  +' --print-coefficients')
if __name__ == '__main__':
    kernal=20
    pool = Pool(processes = kernal)
    gwas_file='/home1/zhaoxz/zhao/brain_volume/volume_sumstats/'
    allfile=os.listdir(gwas_file)
    #allfile='MDD.sumstats.gz'
    n=len(allfile)
    num=0
   # get_h2_res(allfile)
    #get_h2_res(allfile[0])
    for i in range(num,n):
        sub = allfile[i]
        pool.apply_async(get_h2_res,(sub,))
pool.close()
pool.join()