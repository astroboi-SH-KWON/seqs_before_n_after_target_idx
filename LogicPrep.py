import re

import Logic

class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"

    def sort_dict(self, mut_sum_dict):
        result_dict = {}
        for key, vals in mut_sum_dict.items():
            for tmp_list in sorted(vals, key=lambda tmp_list: tmp_list[0]):
                if key in result_dict:
                    result_dict[key].append(tmp_list)
                else:
                    result_dict.update({key: [tmp_list]})
        return result_dict
