import read_cfm_predict_a1

from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit.Chem import rdMolDescriptors






#  path to cfm-predict output files
#cfm_predict_output_path =  'cfm_out\*.log'
#cfm_predict_output_path =  'cfm_predict_aracyc\*.log'
cfm_predict_output_path =  'output_webkis_metab_ce_cfm\*.log'

# filepath to inchi list
# you need this to retrieve inchi string
#filepath_inchi_list = "cfmpredict_input_aracyc.txt"
filepath_inchi_list = "selected_inchi_webkis201810.txt"


# pathway list
filepath_pathway_list = "arabi_compound_types_list.txt"



########################
# reading pathway list
# this pathway list is created by read_metacy_compoun_pathway, when input for cfm-predict is made.


f = open(filepath_pathway_list, "r")

dic_name_vs_list_pathway = {}

for each_line in f:
    #print each_line

    list_line = each_line.strip().split("\t")

    list_pathway = []
    name =""
    if len(list_line) > 1 :
        name = list_line[0]
        list_pathway = list_line[1].split(",")
        dic_name_vs_list_pathway[name] = list_pathway




################################



#############################
# reading inchi string file ( the input file for cfm-id)

print "inchi list"



fi = open(filepath_inchi_list, "r")

dic_name_vs_inchi = {}

for each_line in fi:

    list_line = each_line.strip().split(" ")

    name =""
    if len(list_line) > 1 :
        name = list_line[0]
        inchi = list_line[1]
        dic_name_vs_inchi[name] = inchi




#######################
### multiple files in folder
######################
import sys
import glob

import json


import collections as cl

# initialize dictionary with collections to keep order of elements.
dic_compound =cl.OrderedDict()





list_entry_json = []

for filename in glob.glob(cfm_predict_output_path):


    o_cfmp = read_cfm_predict_a1.read_cfm_predict_a1(filename)

    print "name" , o_cfmp.name
    print o_cfmp.dic_energy_vs_list_mz_int_fragid[2]
    print o_cfmp.dic_fragid_vs_list_fraginfo

    dic_one_entry = {}
    list_as_compound_member_value = []
    dic_for_compound_content = {}

    inchi_string =""


    ####################
    ####################
    # compound+++++++++++++++

    # compound : [     dic_as_compound_member_value , dic_as_compound_member_value ,dic_as_compound_member_value    ]
    # but normally elemnt of list is 1

    dic_as_compound_member_value = {}

    #################
    # inchi

    if not o_cfmp.name in dic_name_vs_inchi.keys():
        print "Inchi data not found in  -name vs inchi- file"
        sys.exit()


    dic_as_compound_member_value["inchi"] = dic_name_vs_inchi[o_cfmp.name]

    ###################
    # Metadata

    list_as_metadata_member_value = []


    ## total exact mass+++++++++++++++++++++

    my_mol = Chem.MolFromInchi( str(dic_name_vs_inchi[o_cfmp.name]))

    dic_total_exact_mass={}

    dic_total_exact_mass["name"] = "total exact mass"
    dic_total_exact_mass["value"] =  rdMolDescriptors.CalcExactMolWt(my_mol)
    list_as_metadata_member_value.append(dic_total_exact_mass)

    ## finish metadata ++++++++++++

    dic_as_compound_member_value["metaData"] = list_as_metadata_member_value





    ################
    # name
    list_as_names_member_value = []
    dic_name = {}
    dic_name["name"] = o_cfmp.name
    list_as_names_member_value.append(dic_name)

    dic_as_compound_member_value["names"] = list_as_names_member_value

    ################
    # Classification
    # "compound":  [  {  classification : XXX , }   ]

    list_as_classification_member_value = []

    dic_classification = {}
    dic_classification["name"] = "kingdom"
    dic_classification["value"] = "notavailableforcfmpredict"
    list_as_classification_member_value.append(dic_classification)



    #++++++++++++++++
    #!!!!!!!!!!!!!!
    #  here we add pathway info to classification temporarily....
    # this is because spectral_network_generator use classification for making multilayer network

    #  pathway


    if o_cfmp.name in dic_name_vs_list_pathway.keys():
        for pth in dic_name_vs_list_pathway[o_cfmp.name] :

            dic_classification = {}
            dic_classification["name"] = "pathway"
            dic_classification["value"] = pth
            list_as_classification_member_value.append(dic_classification)


    dic_as_compound_member_value["classification"] = list_as_classification_member_value





    ########################
    # .....  compound member value finished
    ########################
    # in case you have multiple compound entry in one spectrum entry you have to modify here.
    dic_one_entry["compound"] = [dic_as_compound_member_value]


    ####################
    # spectrum
    #####################


    # get spectrum of particular energy (0,1,2)
    list_mz_int_fragid = o_cfmp.dic_energy_vs_list_mz_int_fragid[3]

    str_for_pk_list = ""
    pk=""
    for n in range( 0 , len(list_mz_int_fragid) ):

        pk = str(   list_mz_int_fragid[n][0]) + ":"  +  str(    list_mz_int_fragid[n][1]   )
        print pk

        str_for_pk_list = str_for_pk_list + pk

        # avoid adding extra " " to the end of line
        if not n == len( list_mz_int_fragid) -1 :
            str_for_pk_list = str_for_pk_list +  " "

    print str_for_pk_list
    dic_one_entry["spectrum"] = str_for_pk_list

    #dic_one_entry[
    #    "spectrum"] = "243:3.252998 254:6.595885 255:21.797459 270:7.047031 271:37.950510 272:3.855091 298:8.320754 299:30.438765 300:100.000000 301:25.906955"



    ########################
    # spectrum with annotation
    ###################################



    # get spectrum of particular energy (0,1,2)
    list_mz_int_fragid = o_cfmp.dic_energy_vs_list_mz_int_fragid[3]

    str_for_pk_list_w_anno_id = ""
    pk=""
    for n in range( 0 , len(list_mz_int_fragid) ):

        pk = str(   list_mz_int_fragid[n][0]) + ":"  +  str(    list_mz_int_fragid[n][1]   ) + ":"  +   ','.join(str(x) for x in list_mz_int_fragid[n][2])
        print "PPKK* " , pk

        str_for_pk_list_w_anno_id = str_for_pk_list_w_anno_id + pk

        # avoid adding extra " " to the end of line
        if not n == len( list_mz_int_fragid) -1 :
            str_for_pk_list_w_anno_id = str_for_pk_list_w_anno_id +  " "

    print "str_for_pk_list_w_anno_id:  " , str_for_pk_list_w_anno_id
    dic_one_entry["spectrum_w_anno_id"] = str_for_pk_list_w_anno_id


    ################################
    # for peak annotation
    ##############################
    count = 0
    str_for_id_vs_annotate = ""
    for key, value in o_cfmp.dic_fragid_vs_list_fraginfo.iteritems():

        str_for_id_vs_annotate = str_for_id_vs_annotate  + str(key) + ":" + str(value[1])

        # avoid adding extra " " to the end of line
        if not count == len( o_cfmp.dic_fragid_vs_list_fraginfo) -1 :
            str_for_id_vs_annotate = str_for_id_vs_annotate +  " "

        count = count + 1

    print "str_for_id_vs_annotate" , str_for_id_vs_annotate
    dic_one_entry["peak_annotation"] = str_for_id_vs_annotate
    #dic_one_entry[
    #    "spectrum"] = "243:3.252998 254:6.595885 255:21.797459 270:7.047031 271:37.950510 272:3.855091 298:8.320754 299:30.438765 300:100.000000 301:25.906955"





    ####################
    ####################
    # id
    dic_one_entry["id"] = o_cfmp.name


    #####################
    ####################
    # metaData

    list_as_metadata_member_value = []

    #######################
    # precursor m/z
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!! now precursor m/z is simply exact mass + H. its for POSITIVE MODE!!!!!!!!!!!!!!!!!!!!!!!!

    dic_precursor_mz = {}
    dic_precursor_mz["name"] = "precursor m/z"
    dic_precursor_mz["value"] = rdMolDescriptors.CalcExactMolWt(my_mol)  + 1.00728

    list_as_metadata_member_value.append(dic_precursor_mz)


    #######################
    #  instrument
    str_instrument = "Orbitrap"
    dic_instrument ={}
    dic_instrument["name"] = str_instrument
    list_as_metadata_member_value.append(dic_instrument)


    #####################
    # ionization

    str_ionization = "ESI"
    dic_ionization ={}
    dic_ionization["name"] = str_ionization
    list_as_metadata_member_value.append(dic_ionization)

    ######################
    # adduct
    #
    str_adduct_type =   "M+H"
    dic_adduct ={}
    dic_adduct["name"] = str_adduct_type
    list_as_metadata_member_value.append(dic_adduct)

    ####
    # finishe meatadata (of spec entry, not compound)

    dic_one_entry["metaData"] = list_as_metadata_member_value

    list_entry_json.append(dic_one_entry)



print "len" , len(list_entry_json)
json_str = json.dumps(list_entry_json,indent=2)
print "----------\n"
print json_str

print "\n"

for j in list_entry_json:
    print j

json_file_path = "json_made.json"
fj = open(json_file_path , "w")
fj.write(json_str)