import read_cfm_predict_a1

######################################################################
# specify file path of single cfm--predict output and ge as object
######################################################################

import read_cfm_predict_a1

filename = "E:/Dropbox/Dropbox/Programming_drpbx/python/mass_spec_related/files/DBXBTMSZEOQQDU-VKHMYHEASA-N.log"
o = read_cfm_predict_a1.read_cfm_predict_a1(filename)

print o.name
print o.list_mz_int
print o.dic_energy_vs_list_mz_int_fragid[3]
print o.dic_fragid_vs_list_fraginfo



######################################################################
# specify file path of single cfm--predict output and get as spec_uni_a2 object.
# filename will be used as name ob the object
# refer to spec_uni_a2 soure for the details
######################################################################
print "\nas spec_uni_a2 object"
s = read_cfm_predict_a1.get_as_spec_uni_a2(filename)
""":type: spec_uni_a2.spectrum_uni_class"""
print s.name
print s.peak_list_mz_int_rel
print s.peak_list_mz_int_abs

print s.peak_list_mz_annoid

print s.dic_annoid_struct


######################################################################
# specify the folder path containing multi[le cfm-predict files and get as spec_uni_a2 object.
# filename will be used as name ob the object
# refer to spec_uni_a2 soure for the details
######################################################################
import glob

print "\n\njmultiple fiiles in folder"

cfm_predict_output_path =  "E:/Dropbox/Dropbox/Programming_drpbx/python/mass_spec_related/files/*.log"

for filename in glob.glob(cfm_predict_output_path):

    s = read_cfm_predict_a1.get_as_spec_uni_a2(filename)
    """:type: spec_uni_a2.spectrum_uni_class"""
    print s.name
    print s.inchi_key
    print s.peak_list_mz_int_rel
    print s.peak_list_mz_int_abs

    print s.peak_list_mz_annoid

    print s.dic_annoid_struct

