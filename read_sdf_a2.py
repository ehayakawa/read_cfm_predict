from rdkit import Chem
from rdkit.Chem import AllChem,Draw,Descriptors
from rdkit import rdBase
import math

# check http://cheminformist.itmol.com/rdkit/file/tag


'''
class class_candidate_mol:

    def __init__(self):
            self.pubchem_cid    =    0
            self.smiles    =    ()
            self.exact_mol_weight = 0
'''

class class_cmpd_from_sdf:
    def __init__(self):
            self.name = ""
            self.db_id    =    0
            self.smiles    =    ""
            self.exact_mol_weight = 0
            self.inchi = ""
            self.inchi_key = ""
            self.list_CAS_registory_numbers = []
            self.link = ""
            self.rdkit_mol = Chem

#####
#  this is for Non lite file.
#   I could not solve the problem witherror MolToInchi.

def read_sdf_a2_ChEBI(sdf_filename):
    suppl = Chem.SDMolSupplier(sdf_filename)
    list_obj_candidate_mol = []

    count = 0
    count_valid = 0
    for mol in suppl:
        if count%100 ==0 :print count
        if mol is not None :

            obj_cmpd_from_sdf = class_cmpd_from_sdf()

            if "ChEBI Name" in list(mol.GetPropNames()):
                obj_cmpd_from_sdf.name = mol.GetProp("ChEBI Name")

            #print mol.GetProp("ChEBI Name")


            ###
            ### Note if thre is description about structures and properties, prioritize them.

            # ID.  Chebi provide id differently when downloading batch or separately ?????

            if "ChEBI ID" in list(mol.GetPropNames()):
                obj_cmpd_from_sdf.db_id = str(mol.GetProp("ChEBI ID"))
            if "ID" in list(mol.GetPropNames()):
                obj_cmpd_from_sdf.db_id = str(mol.GetProp("ID"))

            # inchi and related identifier
            if "InChI" in list(mol.GetPropNames()):
                obj_cmpd_from_sdf.inchi = str(mol.GetProp("InChI"))
            if "INCHIKEY" in list(mol.GetPropNames()):
                obj_cmpd_from_sdf.inchi_key = str(mol.GetProp("INCHIKEY"))
            if "InChIKey" in list(mol.GetPropNames()):
                obj_cmpd_from_sdf.inchi_key = str(mol.GetProp("InChIKey"))
            if "SMILES" in list(mol.GetPropNames()):
                obj_cmpd_from_sdf.smiles = str(mol.GetProp("SMILES"))

            if "CAS Registry Numbers" in list(mol.GetPropNames()):
                obj_cmpd_from_sdf.list_CAS_registory_numbers = str(mol.GetProp("CAS Registry Numbers"))


            ####
            # now it is disabled since some entries have more than one exact mol weight
            # see CHEBI:3385
            """
            if "Monoisotopic Mass" in list(mol.GetPropNames()):
                obj_cmpd_from_sdf.exact_mol_weight = float(mol.GetProp("Monoisotopic Mass"))
            """

            """
            # trying to get info from mol. but later rewritten. See below
            if len(obj_cmpd_from_sdf.smiles) < 1:
                obj_cmpd_from_sdf.smiles = Chem.MolToSmiles(mol)

            if len(obj_cmpd_from_sdf.inchi) < 1  and mol is not None :
                obj_cmpd_from_sdf.inchi = Chem.MolToInchi(mol)
            """

            """
            if len(obj_cmpd_from_sdf.inchi) < 1:
                if mol is None : print "none"
                obj_cmpd_from_sdf.inchi = Chem.MolToInchi(mol)
            """

            # try to create mol from inchi
            obj_cmpd_from_sdf.rdkit_mol = Chem.MolFromInchi(obj_cmpd_from_sdf.inchi)

            if obj_cmpd_from_sdf.rdkit_mol == None :
                obj_cmpd_from_sdf.rdkit_mol = Chem.MolFromSmiles(obj_cmpd_from_sdf.smiles)

            if obj_cmpd_from_sdf.exact_mol_weight == 0:
                if obj_cmpd_from_sdf.rdkit_mol is not None :
                    obj_cmpd_from_sdf.exact_mol_weight =  Descriptors.ExactMolWt(obj_cmpd_from_sdf.rdkit_mol)






            # Note !!!!!!!!!!!!
            # Chebi sdf itself lacks stereo info in SDF, its better to retain Inchi Directly from Inchi string

            # in this version, only entries that give valid rdkit mol will be kept
            if obj_cmpd_from_sdf.rdkit_mol is not None:
                list_obj_candidate_mol.append(obj_cmpd_from_sdf)
                count_valid = count_valid  + 1
        count = count + 1

    print "TOTAL" , count
    print "TOTAL VALID" , count_valid
    return list_obj_candidate_mol


def read_sdf_a2_chemspider(sdf_filename):
    suppl = Chem.SDMolSupplier(sdf_filename)

    list_obj_candidate_mol = []

    for mol in suppl:


        print "name",mol.GetProp("_Name")
        print "< name"
        obj_cmpd_from_sdf = class_cmpd_from_sdf()
        #Draw.MolToFile(mol,'cdk2_mol1.png')
        my_mass  = Descriptors.ExactMolWt(mol)
        obj_cmpd_from_sdf.name = mol.GetProp("_Name")
        obj_cmpd_from_sdf.smiles =  Chem.MolToSmiles(mol)
        obj_cmpd_from_sdf.exact_mol_weight = float(Descriptors.ExactMolWt(mol))
        obj_cmpd_from_sdf.db_id = int(mol.GetProp("CSID"))
        obj_cmpd_from_sdf.inchi = str(mol.GetProp("InChI"))
        obj_cmpd_from_sdf.link = str(mol.GetProp("CSURL"))

        list_obj_candidate_mol.append(    obj_cmpd_from_sdf )

    print "-----------------"

    # example with aspirin
    '''
    aspirin = suppl[0]
    for key in aspirin.GetPropNames():
        value = aspirin.GetProp(key)
        print ">",key
        print value
    '''

    # accessing pubchemcid ---------------------------------------


    return list_obj_candidate_mol



def get_candidate_in_db_by_mass(query_mass, ppm_tol, list_mol):

    list_candidate_matched =[]
    for obj in list_mol:
        if obj.exact_mol_weight > 0 and   (math.fabs(((query_mass - obj.exact_mol_weight)/obj.exact_mol_weight)*1000000)<ppm_tol):
            list_candidate_matched.append(obj)
    return list_candidate_matched






def print_hhh(input):
    print "hhh"
    return input + 1




