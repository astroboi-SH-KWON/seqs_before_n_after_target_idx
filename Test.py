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
WINDOW_SIZE = 2
MAX_SEQ_LEN = 9 #+ WINDOW_SIZE


############### end setting env ################

def main_YG():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    mut_sum_dict = util.read_txt_dvd_by_tab(WORK_DIR + MUT_FILE)
    sorted_mut_dict = logic_prep.sort_dict(mut_sum_dict)
    # print(sorted_mut_dict)

    result_dict = {}
    for key, val_list in sorted_mut_dict.items():
        idx_val = 0
        len_val = len(val_list)
        result_dict.update({key: {}})

        after_strad_dict = {}
        with open(CHR_PATH + key + ".fa", "r") as f:
            header = f.readline()  # header ignored : >chr19
            print("header : " + header)

            idx = 1
            full_cover_p_str = ""
            full_cover_m_str = ""

            while True:
                c = f.read(1)
                if c == "":
                    break
                if "\n" in c:
                    continue
                elif "\r" in c:
                    continue

                tmp_list = []
                full_cover_p_str = full_cover_p_str + c.upper()
                full_cover_m_str = full_cover_m_str + logic.get_complementary(c.upper())
                if len(full_cover_p_str) > MAX_SEQ_LEN + 2:
                    tmp_list.append(full_cover_p_str[-MAX_SEQ_LEN - 1: -1])
                    tmp_list.append(full_cover_p_str[-MAX_SEQ_LEN - 2: -2])
                    tmp_list.append(full_cover_p_str[-MAX_SEQ_LEN - 3: -3])
                    tmp_list.append(full_cover_m_str[-MAX_SEQ_LEN - 1: -1])
                    tmp_list.append(full_cover_m_str[-MAX_SEQ_LEN - 2: -2])
                    tmp_list.append(full_cover_m_str[-MAX_SEQ_LEN - 3: -3])

                if idx in after_strad_dict:
                    for result_dict_key in after_strad_dict[idx]:
                        result_dict[key][result_dict_key].append(tmp_list[result_dict[key][result_dict_key][1]])

                if len_val != idx_val:
                    if val_list[idx_val][0] == idx:
                        tmp_val = val_list[idx_val]
                        Cellline = tmp_val[1]
                        Genome_Change = tmp_val[2]
                        tmp_key = Cellline + "^" + Genome_Change # HCC827^g.chr1:915845C>A
                        Genomesequence = tmp_val[4]
                        Codon_Change = tmp_val[7] # Codon_Change : c.(223-225)Ggt>Tgt
                        Codon_Change_arr = Codon_Change.split(">")

                        # if Genomesequence in ["Missense_Mutation", "Silent", "Splice_Site"]:
                        if len(Codon_Change_arr) > 1:
                            mut_seq = Codon_Change_arr[1] # mut_seq : Tgt
                            if len(mut_seq) == 3:
                                mut_char = Genome_Change.split(">")[1] # A
                                idx_mut_char = mut_seq.find(mut_char)
                                if idx_mut_char == -1:
                                    idx_mut_char = mut_seq.find(logic.get_complementary(mut_char))
                                    result_dict[key].update({tmp_key: ["-", idx_mut_char + 3, tmp_list[idx_mut_char + 3]]})
                                    if idx_mut_char + MAX_SEQ_LEN in after_strad_dict:
                                        after_strad_dict[idx + 3 + MAX_SEQ_LEN].append(tmp_key)
                                    else:
                                        after_strad_dict.update({idx + 3 + MAX_SEQ_LEN: [tmp_key]})
                                else:
                                    result_dict[key].update({tmp_key: ["+", idx_mut_char, tmp_list[idx_mut_char]]})
                                    if idx_mut_char + MAX_SEQ_LEN in after_strad_dict:
                                        after_strad_dict[idx + 3 + MAX_SEQ_LEN].append(tmp_key)
                                    else:
                                        after_strad_dict.update({idx + 3 + MAX_SEQ_LEN: [tmp_key]})
                                # print(after_strad_dict)



                            """
                            {'chr1': 
                                {'HCC827^g.chr1:915845C>A': ['CCCCTTGGC', 'CCCTTGGCA', 'CCTTGGCAC', 'GGGGAACCG', 'GGGAACCGT', 'GGAACCGTG']
                                , 'COLO320^g.chr1:9324611G>A': ['ACAGCAGCT', 'CAGCAGCTT', 'AGCAGCTTC', 'TGTCGTCGA', 'GTCGTCGAA', 'TCGTCGAAG']
                                }
                            }
                            """
                        #     result_dict[key].update({tmp_key: tmp_list})
                        #


                        idx_val += 1

                idx = idx + 1

def test():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    mut_sum_dict = util.read_txt_dvd_by_tab(WORK_DIR + MUT_FILE)
    for key, val_list in mut_sum_dict.items():
        print("")
        print("chr : " + key + ", len : " + str(len(val_list)))
        tmp_dict = {}
        for val in val_list:
            if val[0] in tmp_dict:
                print(str(val[0]) + " is already")
            else:
                tmp_dict.update({val[0]: 1})
        print("chr : " + key + ", len : " + str(len(tmp_dict)))



start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# test()
main_YG()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))


