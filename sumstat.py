##ldsc_data
import os  
root_dir='/share/inspurStorage/home1/zhaoxz/ldsc/batch2/'
statldsc_dir='/share/inspurStorage/home1/zhaoxz/ldsc/batch2/'
allfile=os.listdir(root_dir)
for i in allfile:
    temfile=os.path.join(root_dir,i)
    str_first=i.split('.')[0]
    outfile=os.path.join(statldsc_dir,str_first)
    os.system('python /share/inspurStorage/home1/zhaoxz/ldsc/ldsc-master/munge_sumstats.py --p pval\(-log10\) --sumstats '+temfile +' --N 33000 --out '+outfile)
#python ./ldsc-master/munge_sumstats.py --p pval --sumstats ./pho_data/0438.txt --N 40000 --out 0438
