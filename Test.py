from time import clock
import re
import numpy as np

import Util
import Logic
import LogicPrep
import Valid

############### start to set env ################
WORK_DIR = "D:/000_WORK/KimYoungGwang/20200519/WORK_DIR/"
CHR_PATH = "D:/000_WORK/FAST_REF/human/hg19/chromosomes/"

MUT_FILE = "Mutation_summary"


############### end setting env ################

def main():
    util = Util.Utils()

    mut_sum_dict = util.read_txt_dvd_by_tab(WORK_DIR + MUT_FILE)

    print(mut_sum_dict)


start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
main()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))


