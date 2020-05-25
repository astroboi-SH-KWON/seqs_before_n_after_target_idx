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
WINDOW_SIZE = 1
MAX_SEQ_LEN = 9

INITIAL_MAIN = [CHR_PATH, MAX_SEQ_LEN, WINDOW_SIZE]
############### end setting env ################

def main_YG():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    mut_sum_dict = util.read_txt_dvd_by_tab(WORK_DIR + MUT_FILE)
    sorted_mut_dict = logic_prep.sort_dict(mut_sum_dict)

    result_dict = logic.get_seqs_bfr_aft_trgt_idx(sorted_mut_dict, INITIAL_MAIN)

    util.make_excel(WORK_DIR + "analyze_hg19", sorted_mut_dict, result_dict)

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
main_YG()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))