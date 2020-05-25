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
    # TODO tmp_dict 지우기
    def get_seqs_bfr_aft_trgt_idx(self, data_dict, init):
        chr_path = init[0]
        max_len = init[1]
        win_len = init[2]
        result_dict = {}
        # TODO 정리
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
                    p_m_str_list, full_cover_p_str, full_cover_m_str, n_len_str = self.get_p_m_str_list(full_cover_p_str, full_cover_m_str, max_len, win_len)

                    if idx in n3_after_strd_dict:
                        for result_dict_key in n3_after_strd_dict[idx]:
                            result_dict[chr_key][result_dict_key].append(
                                p_m_str_list[result_dict[chr_key][result_dict_key][1]])
                    if idx in n_mx_len_after_strd_dict:
                        for result_dict_key in n_mx_len_after_strd_dict[idx]:
                            result_dict[chr_key][result_dict_key].append(
                                p_m_str_list[result_dict[chr_key][result_dict_key][1]])
                    if idx in n_len_after_strd_dict:
                        for result_dict_key in n_len_after_strd_dict[idx]:
                            pm_idx = result_dict[chr_key][result_dict_key][1]
                            if result_dict[chr_key][result_dict_key][0] == "-":
                                n_len_str = self.get_complementary_string(n_len_str)
                            result_dict[chr_key][result_dict_key].append(
                                n_len_str + p_m_str_list[pm_idx])

                    if idx in val_list:
                        for tmp_val in val_list[idx]:
                            Cellline = tmp_val[0]
                            Genome_Change = tmp_val[1]
                            tmp_key = Cellline + "^" + Genome_Change  # HCC827^g.chr1:915845C>A
                            Genomesequence = tmp_val[3]
                            cDNA_change = tmp_val[4]
                            Codon_Change = tmp_val[6]  # Codon_Change : c.(223-225)Ggt>Tgt
                            Codon_Change_arr = Codon_Change.split(">")  # ['c.(223-225)Ggt', 'Tgt']

                            # TODO check +/- strand
                            strnd = "+"
                            if cDNA_change == "":
                                strnd = ""
                            elif Genome_Change[-1:] != cDNA_change[-1:]:
                                strnd = "-"

                            # if Genomesequence in ["Missense_Mutation", "Silent", "Splice_Site"]:
                            if len(Codon_Change_arr) > 1:
                                mut_seq = Codon_Change_arr[1]  # mut_seq : Tgt
                                mut_char_arr = Genome_Change.split(">")
                                if len(mut_seq) == 3 and len(mut_char_arr) > 1:
                                    mut_char = mut_char_arr[1]  # A
                                    idx_mut_char = mut_seq.find(mut_char)

                                    if idx_mut_char == -1:
                                        # if mut_char is string not char
                                        tmp_mut = self.get_complementary_string(mut_char)
                                        idx_mut_char = mut_seq.find(tmp_mut)

                                        result_dict[chr_key].update(
                                            {tmp_key: [strnd, idx_mut_char + 3, p_m_str_list[idx_mut_char + 3]]})
                                        if idx + 3 + max_len in n3_after_strd_dict:
                                            # TODO debug
                                            # if tmp_key in tmp_dict:
                                            #     continue
                                            #     print("126")
                                            #     print(tmp_dict[tmp_key])
                                            n3_after_strd_dict[idx + 3 + max_len].append(tmp_key)
                                        else:
                                            # TODO debug
                                            # if tmp_key in tmp_dict:
                                            #     continue
                                            #     print("132")
                                            #     print(tmp_dict[tmp_key])
                                            n3_after_strd_dict.update({idx + 3 + max_len: [tmp_key]})
                                    else:
                                        result_dict[chr_key].update({tmp_key: [strnd, idx_mut_char, p_m_str_list[idx_mut_char]]})
                                        if idx + 3 + max_len in n3_after_strd_dict:
                                            # TODO debug
                                            # if tmp_key in tmp_dict:
                                            #     continue
                                            #     print("141")
                                            #     print(tmp_dict[tmp_key])
                                            n3_after_strd_dict[idx + 3 + max_len].append(tmp_key)
                                        else:
                                            # TODO debug
                                            # if tmp_key in tmp_dict:
                                            #     continue
                                            #     print("143")
                                            #     print(tmp_dict[tmp_key])
                                            n3_after_strd_dict.update({idx + 3 + max_len: [tmp_key]})
                                    # print(n3_after_strd_dict)
                                else:
                                    # TODO mut_seq's len is not 3
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
                                            # TODO debug
                                            # if tmp_key in tmp_dict:
                                            #     continue
                                            #     print("158")
                                            #     print(tmp_dict[tmp_key])
                                            # print(tmp_val)
                                            n_mx_len_after_strd_dict[end_seq + max_len + adder].append(tmp_key)
                                        else:
                                            # TODO debug
                                            # if tmp_key in tmp_dict:
                                            #     continue
                                            #     print("164")
                                            #     print(tmp_dict[tmp_key])
                                            # print(tmp_val)
                                            n_mx_len_after_strd_dict.update({end_seq + max_len + adder: [tmp_key]})

                                    else:
                                        # TODO outrage!!!!!!!!!!!!
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
                                            # TODO debug
                                            # if tmp_key in tmp_dict:
                                            #     continue
                                            #     print("177")
                                            #     print(tmp_dict[tmp_key])
                                            n_len_after_strd_dict[end_seq + max_len + win_len + adder].append(tmp_key)
                                        else:
                                            # TODO debug
                                            # if tmp_key in tmp_dict:
                                            #     continue
                                            #     print("181")
                                            #     print(tmp_dict[tmp_key])
                                            n_len_after_strd_dict.update({end_seq + max_len + win_len + adder: [tmp_key]})
                            elif "_" in Genome_Change:
                                # TODO deli with '_'
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
                                        # TODO debug
                                        # if tmp_key in tmp_dict:
                                        #     continue
                                        #     print("198")
                                        #     print(tmp_dict[tmp_key])
                                        # print(tmp_val)
                                        n_mx_len_after_strd_dict[end_seq + max_len + adder].append(tmp_key)
                                    else:
                                        # TODO debug
                                        # if tmp_key in tmp_dict:
                                        #     continue
                                        #     print("202")
                                        #     print(tmp_dict[tmp_key])
                                        # print(tmp_val)
                                        n_mx_len_after_strd_dict.update({end_seq + max_len + adder: [tmp_key]})

                                else:
                                    # TODO outrage!!!!!!!!!!!!
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
                                        # TODO debug
                                        # if tmp_key in tmp_dict:
                                        #     continue
                                        #     print("217")
                                        #     print(tmp_dict[tmp_key])
                                        n_len_after_strd_dict[end_seq + max_len + win_len + adder].append(tmp_key)
                                    else:
                                        # TODO debug
                                        # if tmp_key in tmp_dict:
                                        #     continue
                                        #     print("221")
                                        #     print(tmp_dict[tmp_key])
                                        n_len_after_strd_dict.update({end_seq + max_len + win_len + adder: [tmp_key]})

                            else:
                                # TODO out of logic
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
                                    # TODO debug
                                    # if tmp_key in tmp_dict:
                                    #     continue
                                    #     print("230_238")
                                    #     print(tmp_dict[tmp_key])
                                    n_len_after_strd_dict[end_seq + max_len + win_len + adder].append(tmp_key)
                                else:
                                    # TODO debug
                                    # if tmp_key in tmp_dict:
                                    #     continue
                                    #     print("230_238")
                                    #     print(tmp_dict[tmp_key])
                                    n_len_after_strd_dict.update({end_seq + max_len + win_len + adder: [tmp_key]})

                    idx = idx + 1

        return result_dict

