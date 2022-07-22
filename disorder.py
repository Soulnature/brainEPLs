##ldsc_data
import os  
from multiprocessing import Pool
# ASD='/home1/zhaoxz/ldsc/pgc_disorders/ASD.txt'
# ADHD='/home1/zhaoxz/ldsc/pgc_disorders/ADHD.txt'
# BIP='/home1/zhaoxz/ldsc/pgc_disorders/BIP.txt'
# MDD='/home1/zhaoxz/ldsc/pgc_disorders/MDD.txt'
# SCZ='/home1/zhaoxz/ldsc/pgc_disorders/SCZ.txt'
# PTSD='/home1/zhaoxz/ldsc/pgc_disorders/PTSD.txt'
disorder_root='/home1/zhaoxz/ldsc/pgc_disorders_9.13/'
disorder_name=['ASD','ALCDEP','BIP','MDD','SCZ','PTSD',
'AUD','OCD','ALCDEP','EPILEPSY','Parkinson','TS','AD_JANSEN','AD_KUNKLE','MS',
'ANXIETY','ASB','INSOMNIA','INTELLIGENCE','RISK_behavior','Neuroticism']
# disorder_name=['MDD','AUD','AD_JANSEN',
# 'ANXIETY','INTELLIGENCE','Neuroticism']
disorder_num=[46351,46568, 51710,500199,79845,
174659,141932,6518,46568,44889,482730,14307,455258,
63926,50000,17310,16400,386533,
269867,466571,380506]
def disease_annote(i):
        tem_files=disorder_name[i]+'.txt'
        disorder_path=os.path.join(disorder_root,tem_files)
        if disorder_name[i]=='Parkinson' or disorder_name[i]=='MS' or disorder_name[i]=='ANXIETY' or disorder_name[i]=='INSOMNIA' or disorder_name[i]=='INTELLIGENCE' or disorder_name[i]=='Neuroticism' :      
            os.system('python /share/inspurStorage/home1/zhaoxz/ldsc/ldsc-master/munge_sumstats.py --sumstats '+disorder_path +' --merge-alleles /home1/zhaoxz/cluster_gwas/eur_w_ld_chr/w_hm3.snplist --a1-inc  --out '+disorder_name[i] )
        else:       
            os.system('python /share/inspurStorage/home1/zhaoxz/ldsc/ldsc-master/munge_sumstats.py --sumstats '+disorder_path +'  --N '+str(disorder_num[i])+' --merge-alleles /home1/zhaoxz/cluster_gwas/eur_w_ld_chr/w_hm3.snplist --a1-inc --out '+disorder_name[i])
kernal=10
num=len(disorder_num)
pool = Pool(processes = kernal)
for i in range(num):
        pool.apply_async(disease_annote,(i,))
pool.close()
pool.join()