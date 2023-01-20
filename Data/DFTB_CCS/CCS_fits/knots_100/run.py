#!/Users/tjams20/opt/anaconda3/bin/python

import os

# os.system("rm *.db")
# os.system("../scripts_DFTB/gen_data_list_DFTB.py")
# os.system("../scripts_DFTB/DB2STRUCTURES.py")
os.system("ccs_fit CCS_input.json")
os.system("../../scripts_DFTB/gen_spl.sh")
os.system("../../scripts_DFTB/tilts_test.py")
os.system("../../scripts_DFTB/NEB_static.py")
