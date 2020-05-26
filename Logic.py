import re

import Util
import LogicPrep

class Logics:
    def __init__(self):
        self.tmp = ""

    def get_complementary(self, c):
        if c == 'C':
            return "G"
        elif c == 'A':
            return "T"
        elif c == 'T':
            return "A"
        elif c == 'G':
            return "C"
        elif c == 'N':
            return "N"
        else:
            print("get_complementary ERROR .... char is [" + c + "]")
            exit()

    def get_complementary_string(self, trgt_str):
        rtrn_str = ""
        for tmp_char in trgt_str:
            rtrn_str += self.get_complementary(tmp_char)
        return rtrn_str

    def check_p_m_strand(self, str_genome, str_change):
        strnd = "+"
        if str_change == "":
            return ""
        elif str_genome[-1:] != str_change[-1:]:
            return "-"
        return strnd

    def get_p_m_str_list(self, full_cover_p_str, full_cover_m_str, max_len, win_len):
        tmp_list = []
        if len(full_cover_p_str) > max_len + win_len:
            tmp_list.append(full_cover_p_str[-max_len - 1: -1])
            tmp_list.append(full_cover_p_str[-max_len - 2: -2])
            tmp_list.append(full_cover_p_str[-max_len - 3: -3])
            tmp_list.append(full_cover_m_str[-max_len - 1: -1])
            tmp_list.append(full_cover_m_str[-max_len - 2: -2])
            tmp_list.append(full_cover_m_str[-max_len - 3: -3])
        return tmp_list, full_cover_p_str[-max_len - win_len - 1:], full_cover_m_str[-max_len - win_len - 1:], full_cover_p_str[1:win_len + 1]

    def get_p_m_str_list_for_n3(self, full_cover_p_str, full_cover_m_str, max_len, win_len):
        tmp_list = []
        if len(full_cover_p_str) > max_len + win_len:
            tmp_list.append(full_cover_p_str[-max_len - 1: -1])
            tmp_list.append(full_cover_p_str[-max_len - 2: -1])
            tmp_list.append(full_cover_p_str[-max_len - 3: -1])
            tmp_list.append(full_cover_m_str[-max_len - 1: -1])
            tmp_list.append(full_cover_m_str[-max_len - 2: -1])
            tmp_list.append(full_cover_m_str[-max_len - 3: -1])
        return tmp_list

    def add_in_frm_forwrd_seq(self, idx, chr_key, after_strd_dict, result_dict, p_m_str_list):
        if idx in after_strd_dict:
            # print(p_m_str_list)
            for result_dict_key in after_strd_dict[idx]:
                # print("result_dict_key : " + result_dict_key)
                # print("p_m_str_list' idx : " + str(result_dict[chr_key][result_dict_key][1]))
                result_dict[chr_key][result_dict_key].append(
                    p_m_str_list[result_dict[chr_key][result_dict_key][1]])
        return result_dict

    def add_out_frm_forwrd_seq(self, idx, chr_key, after_strd_dict, result_dict, p_m_str_list, n_len_str):
        if idx in after_strd_dict:
            for result_dict_key in after_strd_dict[idx]:
                pm_idx = result_dict[chr_key][result_dict_key][1]
                if result_dict[chr_key][result_dict_key][0] == "-":
                    n_len_str = self.get_complementary_string(n_len_str)
                result_dict[chr_key][result_dict_key].append(
                    n_len_str + p_m_str_list[pm_idx])
        return result_dict

    def get_frwrd_str_bckwrd_idx(self, mut_seq, p_m_str_n3_list, forward):
        m_strand_adder = 0
        if forward == "-":
            m_strand_adder = 3
        frwrd_idx = 0
        forward = p_m_str_n3_list[frwrd_idx + m_strand_adder]
        bckwrd_idx = 0
        is_first_low = False
        is_second_low = False
        if mut_seq[0].islower():
            forward = p_m_str_n3_list[frwrd_idx + m_strand_adder]
            frwrd_idx += 1
            is_first_low = True

        if mut_seq[1].islower():
            is_second_low = True
            if is_first_low:
                frwrd_idx += 1
                forward = p_m_str_n3_list[frwrd_idx + m_strand_adder]
            else:
                bckwrd_idx += 1

        if mut_seq[2].islower():
            if not is_first_low:
                bckwrd_idx += 1
            elif not is_second_low:
                bckwrd_idx += 1

        return forward, frwrd_idx, bckwrd_idx

    def get_seqs_bfr_aft_trgt_idx(self, data_dict, init):
        chr_path = init[0]
        max_len = init[1]
        win_len = init[2]
        result_dict = {}

        util = Util.Utils()
        for chr_key, val_list in data_dict.items():
            result_dict.update({chr_key: {}})

            n3_after_strd_dict = {}
            n_len_after_strd_dict = {}
            n_mx_len_after_strd_dict = {}
            with open(chr_path + chr_key + ".fa", "r") as f:
                header = f.readline()  # header ignored : >chr19
                print("header : " + header)

                idx = 1
                full_cover_p_str = ""
                full_cover_m_str = ""

                while True:
                    c = f.read(1)
                    if c == "":
                        print(str(idx))
                        break
                    if "\n" in c:
                        continue
                    elif "\r" in c:
                        continue

                    full_cover_p_str = full_cover_p_str + c.upper()
                    full_cover_m_str = full_cover_m_str + self.get_complementary(c.upper())
                    p_m_str_n3_list = self.get_p_m_str_list_for_n3(full_cover_p_str, full_cover_m_str, max_len, win_len)
                    p_m_str_list, full_cover_p_str, full_cover_m_str, n_len_str = self.get_p_m_str_list(full_cover_p_str, full_cover_m_str, max_len, win_len)

                    # add forward seq into result_dict
                    result_dict = self.add_in_frm_forwrd_seq(idx, chr_key, n3_after_strd_dict, result_dict, p_m_str_n3_list)
                    result_dict = self.add_in_frm_forwrd_seq(idx, chr_key, n_mx_len_after_strd_dict, result_dict, p_m_str_list)
                    result_dict = self.add_out_frm_forwrd_seq(idx, chr_key, n_len_after_strd_dict, result_dict, p_m_str_list, n_len_str)

                    if idx in val_list:
                        for tmp_val in val_list[idx]:
                            Cellline = tmp_val[0]
                            Genome_Change = tmp_val[1]
                            tmp_key = Cellline + "^" + Genome_Change  # HCC827^g.chr1:915845C>A
                            Genomesequence = tmp_val[3]
                            cDNA_change = tmp_val[4]
                            Codon_Change = tmp_val[6]  # Codon_Change : c.(223-225)Ggt>Tgt
                            Codon_Change_arr = Codon_Change.split(">")  # ['c.(223-225)Ggt', 'Tgt']

                            # check +/- strand
                            strnd = self.check_p_m_strand(Genome_Change, cDNA_change)

                            if len(Codon_Change_arr) > 1:
                                mut_seq = Codon_Change_arr[1]  # mut_seq : Tgt
                                mut_char_arr = Genome_Change.split(">")
                                if len(mut_seq) == 3 and len(mut_char_arr) > 1:
                                    # TODO del
                                    mut_char = mut_char_arr[1]  # A
                                    # idx_mut_char = mut_seq.find(mut_char)

                                    if strnd == "-":
                                        forward, frwrd_idx, bckwrd_idx = self.get_frwrd_str_bckwrd_idx(mut_seq, p_m_str_n3_list, strnd)

                                        result_dict[chr_key].update({tmp_key: [strnd, bckwrd_idx + 3, forward]})

                                        if idx + 3 + max_len - frwrd_idx in n3_after_strd_dict:
                                            n3_after_strd_dict[idx + 3 + max_len - frwrd_idx].append(tmp_key)
                                        else:
                                            n3_after_strd_dict.update({idx + 3 + max_len - frwrd_idx: [tmp_key]})
                                    else:
                                        forward, frwrd_idx, bckwrd_idx = self.get_frwrd_str_bckwrd_idx(mut_seq, p_m_str_n3_list, strnd)

                                        result_dict[chr_key].update({tmp_key: [strnd, bckwrd_idx, forward]})

                                        if idx + 3 + max_len - frwrd_idx in n3_after_strd_dict:
                                            n3_after_strd_dict[idx + 3 + max_len - frwrd_idx].append(tmp_key)
                                        else:
                                            n3_after_strd_dict.update({idx + 3 + max_len - frwrd_idx: [tmp_key]})
                                else:
                                    # mut_seq's len is not 3
                                    end_seq = int(re.findall(r'\d+', Genome_Change[Genome_Change.index(":"):])[1])
                                    # if Genomesequence has 'In_Frame'
                                    if 'In_Frame' in Genomesequence:
                                        idx_p_m_str_list = 0
                                        if strnd == "-":
                                            idx_p_m_str_list = 3
                                        forward_str = p_m_str_list[idx_p_m_str_list]
                                        adder = 1
                                        if 'ins'.upper() in Genome_Change.upper():
                                            adder = 0
                                            forward_str = full_cover_p_str[-max_len:]
                                            if strnd == '-':
                                                forward_str = full_cover_m_str[-max_len:]

                                        result_dict[chr_key].update({tmp_key: [strnd, idx_p_m_str_list, forward_str]})
                                        if end_seq + max_len + adder in n_mx_len_after_strd_dict:
                                            n_mx_len_after_strd_dict[end_seq + max_len + adder].append(tmp_key)
                                        else:
                                            n_mx_len_after_strd_dict.update({end_seq + max_len + adder: [tmp_key]})

                                    else:
                                        idx_p_m_str_list = 0
                                        if strnd == "-":
                                            idx_p_m_str_list = 3

                                        forward_str = n_len_str + p_m_str_list[idx_p_m_str_list]
                                        adder = 1
                                        if 'ins'.upper() in Genome_Change.upper():
                                            adder = 0
                                            forward_str = full_cover_p_str[-max_len + win_len:]
                                            if strnd == '-':
                                                forward_str = full_cover_m_str[-max_len + win_len:]

                                        result_dict[chr_key].update(
                                            {tmp_key: [strnd, idx_p_m_str_list, forward_str]})
                                        if end_seq + max_len + win_len + adder in n_len_after_strd_dict:
                                            n_len_after_strd_dict[end_seq + max_len + win_len + adder].append(tmp_key)
                                        else:
                                            n_len_after_strd_dict.update({end_seq + max_len + win_len + adder: [tmp_key]})
                            elif "_" in Genome_Change:
                                # deli with '_'
                                end_seq = int(re.findall(r'\d+', Genome_Change[Genome_Change.index(":"):])[1])
                                # if Genomesequence has 'In_Frame'
                                if 'In_Frame' in Genomesequence:
                                    idx_p_m_str_list = 0
                                    if strnd == "-":
                                        idx_p_m_str_list = 3
                                    forward_str = p_m_str_list[idx_p_m_str_list]
                                    adder = 1
                                    if 'ins'.upper() in Genome_Change.upper():
                                        adder = 0
                                        forward_str = full_cover_p_str[-max_len:]
                                        if strnd == '-':
                                            forward_str = full_cover_m_str[-max_len:]

                                    result_dict[chr_key].update(
                                        {tmp_key: [strnd, idx_p_m_str_list, forward_str]})
                                    if end_seq + max_len + adder in n_mx_len_after_strd_dict:
                                        n_mx_len_after_strd_dict[end_seq + max_len + adder].append(tmp_key)
                                    else:
                                        n_mx_len_after_strd_dict.update({end_seq + max_len + adder: [tmp_key]})

                                else:
                                    idx_p_m_str_list = 0
                                    if strnd == "-":
                                        idx_p_m_str_list = 3

                                    forward_str = n_len_str + p_m_str_list[idx_p_m_str_list]
                                    adder = 1
                                    if 'ins'.upper() in Genome_Change.upper():
                                        adder = 0
                                        forward_str = full_cover_p_str[- max_len - win_len:]
                                        if strnd == "-":
                                            forward_str = full_cover_m_str[- max_len - win_len:]

                                    result_dict[chr_key].update(
                                        {tmp_key: [strnd, idx_p_m_str_list, forward_str]})
                                    if end_seq + max_len + win_len + adder in n_len_after_strd_dict:
                                        n_len_after_strd_dict[end_seq + max_len + win_len + adder].append(tmp_key)
                                    else:
                                        n_len_after_strd_dict.update({end_seq + max_len + win_len + adder: [tmp_key]})

                            else:
                                # out of logic
                                end_seq = idx  # + 1
                                idx_p_m_str_list = 0
                                if strnd == "-":
                                    idx_p_m_str_list = 3

                                forward_str = n_len_str + p_m_str_list[idx_p_m_str_list]
                                adder = 1
                                if 'ins'.upper() in Genome_Change.upper():
                                    adder = 0
                                    forward_str = full_cover_p_str[- max_len - win_len:]
                                    if strnd == "-":
                                        forward_str = full_cover_m_str[- max_len - win_len:]

                                result_dict[chr_key].update(
                                    {tmp_key: [strnd, idx_p_m_str_list, forward_str]})
                                if end_seq + max_len + win_len + adder in n_len_after_strd_dict:
                                    n_len_after_strd_dict[end_seq + max_len + win_len + adder].append(tmp_key)
                                else:
                                    n_len_after_strd_dict.update({end_seq + max_len + win_len + adder: [tmp_key]})

                    idx = idx + 1

        return result_dict





