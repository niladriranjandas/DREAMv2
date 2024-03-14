import sys
import findViolations as fv

def getViolsForProtein(fileslist,protname):
    ''''

    '''
    
    c=0

    list_fileslist = fileslist.split(',')
    for files in list_fileslist:
           c=c+1
           protname_i = '%s_%d'%(protname,c)
           print('\n%s'%(files))
           #len_pseudo_viol_strong,ps_ret_mean,ps_ret_std,len_pseudo_viol_medium,pm_ret_mean,pm_ret_std,len_pseudo_viol_weak,pw_ret_mean,pw_ret_std,len_nonpseudo_viol_strong,nps_ret_mean,nps_ret_std,len_nonpseudo_viol_medium,npm_ret_mean,npm_ret_std,len_nonpseudo_viol_weak,npw_ret_mean,npw_ret_std = fv.parseViolFile_v2(files,protname_i)
           len_pseudo_viol_strong,ps_ret_mean,ps_ret_std,len_pseudo_viol_medium,pm_ret_mean,pm_ret_std,len_pseudo_viol_weak,pw_ret_mean,pw_ret_std,len_nonpseudo_viol_strong,nps_ret_mean,nps_ret_std,len_nonpseudo_viol_medium,npm_ret_mean,npm_ret_std,len_nonpseudo_viol_weak,npw_ret_mean,npw_ret_std = fv.parseViolFile(files,protname_i)

           print('Pseudo: strong: (%f,%f,%f) medium: (%f,%f,%f) weak: (%f,%f,%f)'%(len_pseudo_viol_strong,ps_ret_mean,ps_ret_std,len_pseudo_viol_medium,pm_ret_mean,pm_ret_std,len_pseudo_viol_weak,pw_ret_mean,pw_ret_std))
           print('Non-Pseudo: strong: (%f,%f,%f) medium: (%f,%f,%f) weak: (%f,%f,%f)'%(len_nonpseudo_viol_strong,nps_ret_mean,nps_ret_std,len_nonpseudo_viol_medium,npm_ret_mean,npm_ret_std,len_nonpseudo_viol_weak,npw_ret_mean,npw_ret_std))



def main():
    getViolsForProtein(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()

