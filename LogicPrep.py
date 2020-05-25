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
            result_dict.update({key: {}})
            for tmp_list in sorted(vals, key=lambda tmp_list: tmp_list[0]):
                seq_key = tmp_list[0]
                data_list = tmp_list[1:]
                if seq_key in result_dict[key]:
                    result_dict[key][seq_key].append(data_list)
                else:
                    result_dict[key].update({seq_key: [data_list]})
        return result_dict
