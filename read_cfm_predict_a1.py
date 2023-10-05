__author__ = 'esk'

# last modified 20190104
#    error handling when no file was found

import re
import copy
import spectrum_uni_a2

import ntpath
import sys

from pathlib2 import Path


##########
#  NOTE  idx 3 of dic_energy_vs_list_mz_int_fragid is merged spectrum.






class Cfm_predict_res_a1:
    def __init__(self):
        self.name = ""
        self.list_mz_int = []
        ##   energy0: [100,50]
        self.dic_energy_vs_list_mz_int_fragid_fragscore = {}
        ###   id :  [   mass , SMILES]
        self.dic_fragid_vs_list_fraginfo = {}
        self.note = ""


debug = 0


#  read_cfm_predict_a1
#       core function. receive filepath (single filename)of cfm-predict

# returns object:  print o.name
#  *  dic_energy_vs_list_mz_int_fragid   :
#        dictonary where key is id of energy state (normally 0,1,2). 3 is all combined version.
#       value is list of  peak info consisting ( mz, intensity, list of frag id)
#  *   dic_fragid_vs_list_fraginfo :
#         dictionary where key is frag id and value is [ mz,  frag structure].



def read_cfm_predict_a1(filepath, score_th = 0 , flag_exit_if_file_not_found=0 , debug = 0 ):
    # read and process only if file actually exists
    my_file = Path(filepath)

    predict_o = Cfm_predict_res_a1()

    if my_file.is_file():

        # opening file
        f = open(filepath, "r")

        if debug: print "filepathh", filepath
        # check file path
        if debug: print "splitted", filepath.split("\\")

        filename = ""

        # check if the specifile file(path) is a path incl. folder  or just filename

        # if file path is actually path to particular folder.
        if (len(filepath.split("\\")) > 1):
            if debug: print "rel"
            filename = filepath.split("\\")[len(filepath.split("\\")) - 1]

        # or just filename
        else:
            filename = filepath

        if debug: print "filename:", filename
        cmpdname = filename.split(".")[0]

        if debug: print cmpdname

        pattern = r"\."

        re_pat = re.compile(pattern)

        list_mz_int_list_fragid_fragscore_energy0 = []
        list_mz_int_list_fragid_fragscore_energy1 = []
        list_mz_int_list_fragid_fragscore_energy2 = []

        dic_fragid_vs_list_fraginfo = {}

        count = 0
        energy = 0
        for each_line in f:
            # print each_line


            if each_line == "energy0\n":
                energy = 0
            if each_line == "energy1\n":
                energy = 1

            if each_line == "energy2\n":
                energy = 2

            list_contents = each_line.split(" ")
            # print list_contents[0]

            # checking if the first element is float "xx.xx"
            matchOB_a = re_pat.search(list_contents[0])

            list_fragid = []
            list_fragscore = []


            # if the current line is for fragment peak content..
            reg_score = re.compile("\\(([0-9.\s]*)\\)")
            if matchOB_a:

                # getting score    [xx xx]
                match = reg_score.search(each_line)
                if match != None:
                    if debug: print match.group(1)
                    list_fragscore =  [ float(x)  for x in match.group(1).split(" ") ]
                    if debug: print l_score

                each_line_blktcontent_removed = re.sub(r"\(.*\)", "", each_line)
                list_contents = each_line_blktcontent_removed.rstrip().split(" ")
                count = 0

                for con in list_contents:
                    matchOB_b = re_pat.search(con)
                    if count > 1:
                        list_fragid.append(int(con))
                    count = count + 1

                if debug: print "double"
                if debug: print "mz", list_contents[0], "int", list_contents[1]

                mz = float(list_contents[0])
                intensity = float(list_contents[1])
                if debug: print list_fragid

                if energy == 0:
                    list_mz_int_list_fragid_fragscore_energy0.append([mz, intensity, list_fragid,list_fragscore])
                if energy == 1:
                    list_mz_int_list_fragid_fragscore_energy1.append([mz, intensity, list_fragid,list_fragscore])
                if energy == 2:
                    list_mz_int_list_fragid_fragscore_energy2.append([mz, intensity, list_fragid,list_fragscore])

            # if the current line is for fragment struct information

            # 231 257.1899918 C#CC#CC#CC=CCCCCCCCCCC[OH2+]
            if (not each_line.startswith("e")) and (not matchOB_a) and (len(each_line) > 10):

                if debug: print "now", each_line
                if debug: print "FRAG"
                if debug: print list_contents[0], list_contents[1], list_contents[2]

                dic_fragid_vs_list_fraginfo[int(list_contents[0])] = [float(list_contents[1]), list_contents[2].strip()]

        if debug: print list_mz_int_list_fragid_fragscore_energy0
        if debug: print list_mz_int_list_fragid_fragscore_energy1
        if debug: print list_mz_int_list_fragid_fragscore_energy2

        # create merged spectrum
        list_mz_int_list_fragid_fragscore_COMBINED = []

        list_l = [ list_mz_int_list_fragid_fragscore_energy0, list_mz_int_list_fragid_fragscore_energy1, list_mz_int_list_fragid_fragscore_energy2]

        count_l = -1
        for li in list_l:
            count_l = count_l + 1
            print "=============now" , count_l , "th list"
            for pk_i in li:

                flag_match = 0

                #for pk_c in list_mz_int_list_fragid_fragscore_COMBINED:
                for i in range(0,  len(list_mz_int_list_fragid_fragscore_COMBINED)):

                    pk_c = list_mz_int_list_fragid_fragscore_COMBINED[i]
                    # if mz match
                    if pk_c[0] == pk_i[0]:
                        flag_match = 1
                        print "\nMATCHED  pk_c mz " , pk_c[0] , pk_i[0]
                        # accumulate intensity
                        pk_c[1] = pk_c[1] + pk_i[1]

                        #pk_c[2] = pk_c[2] + pk_i[2]

                        l_fragid_c = pk_c[2]
                        l_fragscore_c =   pk_c[3]
                        # idx 2 is frag id
                        pk_c[2] = pk_c[2] + pk_i[2]
                        # idx 3 is fragscore
                        pk_c[3] = pk_c[3] + pk_i[3]

                # if this is new mz, just append
                if flag_match == 0:
                    print "new peak" , pk_i
                    list_mz_int_list_fragid_fragscore_COMBINED.append(pk_i)


        list_mz_int_list_fragid_fragscore_COMBINED.sort(key=lambda e: e[0])

        print "list_mz_int_list_fragid_fragscore_COMBINED"
        print list_mz_int_list_fragid_fragscore_COMBINED

        ###########################################
        # making  fragid (and fragscore) non redundatnt and sync
        list_mz_int_list_fragid_fragscore_COMBINED_mod = []

        for pk in list_mz_int_list_fragid_fragscore_COMBINED :

            l_fragid = pk[2]
            l_fragscore = pk[3]
            l_fragid_nr = []
            l_fragscore_nr = []
            for i in range( 0, len(l_fragid)) :

                # if the frag id is new one
                if l_fragid[i] not in l_fragid_nr :
                    l_fragid_nr.append( l_fragid[i] )
                    l_fragscore_nr.append(l_fragscore[i])

                # if the fragid is no new( present in l_nr), but the score is higher, then update with new score
                if ( l_fragid[i] in l_fragid_nr  ) :
                    # get index of l_nr where same frag id can be found
                    idx_of_l_nr_hit_id = l_fragid_nr.index(l_fragid[i])
                    if (    l_fragscore[i]    >    l_fragscore_nr[idx_of_l_nr_hit_id]  )   :
                        l_fragscore_nr[idx_of_l_nr_hit_id] =  l_fragscore[i]

            list_mz_int_list_fragid_fragscore_COMBINED_mod.append(  [pk[0], pk[1], l_fragid_nr, l_fragscore_nr])


        list_all_list_mz_int_list_fragid_fragscore = [ list_mz_int_list_fragid_fragscore_energy0, list_mz_int_list_fragid_fragscore_energy1, list_mz_int_list_fragid_fragscore_energy2, list_mz_int_list_fragid_fragscore_COMBINED_mod ]




        ########
        # thresholding with frag score

        list_all_list_mz_int_list_fragid_fragscore_th = []

        for list_mz_int_list_fragid_fragscore in list_all_list_mz_int_list_fragid_fragscore:
            list_mz_int_list_fragid_fragscore_th =[]

            for pk in list_mz_int_list_fragid_fragscore :
                l_fragid = pk[2]
                l_fragscore = pk[3]
                l_fragid_th = []
                l_fragscore_th = []

                f_any_fragid_pass_th = 0
                for i in range( 0, len(l_fragid)) :
                    if l_fragscore[i] > score_th :
                        f_any_fragid_pass_th = 1
                        l_fragid_th.append( l_fragid[i]  )
                        l_fragscore_th.append(l_fragscore[i])

                if f_any_fragid_pass_th == 1 :
                    list_mz_int_list_fragid_fragscore_th.append( [pk[0], pk[1], l_fragid_th, l_fragscore_th] )

            list_all_list_mz_int_list_fragid_fragscore_th.append(list_mz_int_list_fragid_fragscore_th)





        predict_o.name = cmpdname

        # making dictionary that hold list of peaks in different energy
        count = -1
        for l_pk in list_all_list_mz_int_list_fragid_fragscore_th :
            count = count + 1
            predict_o.dic_energy_vs_list_mz_int_fragid_fragscore[count] = l_pk


        predict_o.dic_fragid_vs_list_fraginfo = dic_fragid_vs_list_fraginfo

    if not my_file.is_file():
        predict_o.note = "NO_CORRESPONDING_FILE_FOUND"

        if flag_exit_if_file_not_found == 1:
            print "ERROR from read_cfm_predict_a1 read_cfm_predict_a1 function.  "
            print "it seems NO corresponding file was found"
            print "this module expect file path like ....  E:/Dropbox/Dropbox/Programming_drpbx/python/mass_spec_related/files/DBXBTMSZEOQQDU-VKHMYHEASA-N.log "
            sys.exit()
    return predict_o


# Note !!! if there is no corresponding file, it returns empty spec_uni_a2
#
#  flag_return_empty_if_no_file_found :  1 : it returns empty spec object.  2: strict. ERROR and exit if no corresponding file found.
def get_as_spec_uni_a2(filepath, id_energy=3 , score_th = 0 , flag_return_empty_if_no_file_found = 1):
    o = read_cfm_predict_a1(filepath , score_th = score_th )

    peak_list_mz_int_abs = []
    peak_list_mz_annoid = []
    peak_list_mz_fragscore = []
    dic_annoid_struct = {}

    # create spec_uni_a2 object
    su = spectrum_uni_a2.spectrum_uni_class()

    ###
    # note if there is no file, dic_energy_vs_list_mz_int_fragid will have 0 length and results in error in the following part
    # This flag is intended to avoid that/

    flag_peaklist_present = 0

    # if the cfm predict spec object has valid dictionary,
    if id_energy in o.dic_energy_vs_list_mz_int_fragid_fragscore:
        flag_peaklist_present = 1


    """
    if flag_peaklist_present == 1:

        for l in o.dic_energy_vs_list_mz_int_fragid_fragscore[id_energy]:
            peak_list_mz_int_abs.append([l[0], l[1]])
            peak_list_mz_annoid.append([l[0], l[2]])
            peak_list_mz_fragscore.append([l[0], l[3]])

        for k, v in o.dic_fragid_vs_list_fraginfo.iteritems():
            dic_annoid_struct[k] = v[1]
    """
    if flag_peaklist_present == 1:

        for l in o.dic_energy_vs_list_mz_int_fragid_fragscore[id_energy]:



            peak_list_mz_int_abs.append([l[0], l[1]])
            peak_list_mz_annoid.append([l[0], l[2]])
            peak_list_mz_fragscore.append([l[0], l[3]])

        for k, v in o.dic_fragid_vs_list_fraginfo.iteritems():
            dic_annoid_struct[k] = v[1]


    su.peak_list_mz_int_abs = peak_list_mz_int_abs
    spectrum_uni_a2.set_peak_list_rel_int(su, max_value=100)

    su.peak_list_mz_annoid = peak_list_mz_annoid
    su.dic_annoid_struct = dic_annoid_struct

    # trying to give name to spec uni object based on filename
    filename = ntpath.basename(filepath)

    reg_filename = re.compile("(.*)\.,*")
    match = reg_filename.search(filename)

    ##
    ## MAtch is None if there is no corresponding cfm-pred file.
    if match == None  and flag_return_empty_if_no_file_found == 0 :
        print "ERROR from read_cfm_predict_a1 get_as_spec_uni_a2 function.  "
        print "it seems something wrong with file path or file name "
        print "this module expect file path like ....  E:/Dropbox/Dropbox/Programming_drpbx/python/mass_spec_related/files/DBXBTMSZEOQQDU-VKHMYHEASA-N.log "
        sys.exit()

    if match != None:
        su.name = match.group(1)
        su.inchi_key = match.group(1)

    # if_correspoding file is not found,
    if match == None:
        su.name = "THEO_SPEC_FILE_NOT_FOUND"
        su.inchi_key = "EMPTY_SINCE_FILE_NOT_FOUND"

    return su


"""
def read_specific_file_from_folder(foldername, filename_specified, flag_ignore_ext=1):

    flag_match = 0
    count_match = 0
    filepath_to_return  = ""
    for filepath in glob.glob(foldername + "\*"):
        flag_match = 0
        # print filepath , " "
        # print ntpath.basename(filepath)

        filename = ntpath.basename(filepath)
        filename_wo_ext = str(ntpath.basename(filepath).split(".")[0].strip("\n"))

        if flag_ignore_ext == 1:
            if filename_wo_ext == filename_specified:
                flag_match = 1

        if flag_ignore_ext == 0:
            if filename == filename_specified:
                flag_match = 1

        if flag_match :
            print "found" , filename_specified
            count_match = count_match + 1
            filepath_to_return = filepath
"""