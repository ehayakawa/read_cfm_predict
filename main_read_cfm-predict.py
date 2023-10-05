import read_cfm_predict_a1



filename = "AAEVYOVXGOFMJO-UHFFFAOYSA-N.log"
o = read_cfm_predict_a1.read_cfm_predict_a1(filename)

print o.name
print o.dic_energy_vs_list_mz_int_fragid[2]
print o.dic_fragid_vs_list_fraginfo


#######################
### multiple files in folder
######################

"""
import glob

for filename in glob.glob('cfm_out\*.log'):


    o = read_cfm_predict_a1.read_cfm_predict_a1(filename)

    print "name" , o.name
    print o.dic_energy_vs_list_mz_int_fragid[3]
    print o.dic_fragid_vs_list_fraginfo
"""