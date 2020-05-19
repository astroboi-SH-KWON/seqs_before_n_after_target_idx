from os import listdir
from os.path import isfile, join
import pandas as pd
import openpyxl
from time import clock
import random
import math

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    def read_txt_dvd_by_tab(self, path):
        result_dict = {}
        with open(path + self.ext_txt, "r") as f:
            head_list = f.readline().replace("\n", "").split("\t")
            idx = 1
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == '':
                    break
                val_list = tmp_line.split("\t")
                # result_dict.update({idx: {
                #     head_list[0]: val_list[0]
                #     , head_list[1]: val_list[1]
                #     , head_list[2]: val_list[2]
                #     , head_list[3]: val_list[3]
                #     , head_list[4]: val_list[4]
                #     , head_list[5]: val_list[5]
                #     , head_list[6]: val_list[6]
                # }})
                gene_chng_list = val_list[1].replace("g.","").split(":")
                if gene_chng_list[0] in result_dict:
                    tmp_list = [gene_chng_list[1].split("_")[0]]
                    result_dict[gene_chng_list[0]].append(tmp_list + val_list)
                else:
                    tmp_list = [gene_chng_list[1]]
                    result_dict.update({gene_chng_list[0]: [tmp_list + val_list]})
                idx += 1

        return result_dict

