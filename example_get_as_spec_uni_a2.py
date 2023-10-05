print "\nas spec_uni_a2 object"
import read_cfm_predict_a1

filename = "E:/Dropbox/Dropbox/Programming_drpbx/python/mass_spec_related/files/DBXBTMSZEOQQDU-VKHMYHEASA-N.log"


s = read_cfm_predict_a1.get_as_spec_uni_a2(filename)
""":type: spec_uni_a2.spectrum_uni_class"""
print s.name
print s.peak_list_mz_int_rel
print s.peak_list_mz_int_abs

print s.peak_list_mz_annoid

print s.dic_annoid_struct