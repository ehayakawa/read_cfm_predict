import read_cfm_predict_a1

filename = "AAEVYOVXGOFMJO-UHFFFAOYSA-Nb.log"
o = read_cfm_predict_a1.read_cfm_predict_a1(filename, score_th = 1)

print o.name
print o.list_mz_int
print o.dic_energy_vs_list_mz_int_fragid_fragscore[3]
print o.dic_fragid_vs_list_fraginfo


