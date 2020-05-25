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
    # print(sorted_mut_dict)
    # sorted_mut_dict = {'chr11': {614355: [['HCC827', 'g.chr11:614355G>C', 'IRF7', 'Silent', 'c.537C>G', 'p.L179L', 'c.(535-537)ctC>ctG']], 4936194: [['COLO320', 'g.chr11:4936194C>T', 'OR51G2', 'Missense_Mutation', 'c.700G>A', 'p.A234T', 'c.(700-702)Gcc>Acc']], 4944864: [['COLO320', 'g.chr11:4944864G>A', 'OR51G1', 'Nonsense_Mutation', 'c.706C>T', 'p.R236*', 'c.(706-708)Cga>Tga']], 5269657: [['COLO320', 'g.chr11:5269657C>T', 'HBG1', 'Missense_Mutation', 'c.376G>A', 'p.E126K', 'c.(376-378)Gag>Aag']], 6048332: [['NCIH522', 'g.chr11:6048332A>T', 'OR56A1', 'Silent', 'c.603T>A', 'p.L201L', 'c.(601-603)ctT>ctA']], 6578097: [['NCIH522', 'g.chr11:6578097C>T', 'DNHD1', 'Silent', 'c.7572C>T', 'p.V2524V', 'c.(7570-7572)gtC>gtT']], 6585903: [['NCIH522', 'g.chr11:6585903G>A', 'DNHD1', 'Missense_Mutation', 'c.10625G>A', 'p.R3542Q', 'c.(10624-10626)cGa>cAa']], 6644232: [['COLO320', 'g.chr11:6644232C>T', 'DCHS1', 'Missense_Mutation', 'c.8675G>A', 'p.R2892Q', 'c.(8674-8676)cGg>cAg']], 7981578: [['HCC827', 'g.chr11:7981578C>T', 'NLRP10', 'Silent', 'c.1581G>A', 'p.S527S', 'c.(1579-1581)tcG>tcA']], 12023940: [['NCIH522', 'g.chr11:12023940G>C', 'DKK3', 'Missense_Mutation', 'c.258C>G', 'p.N86K', 'c.(256-258)aaC>aaG']], 14825503: [['COLO320', 'g.chr11:14825503C>G', 'PDE3B', 'Missense_Mutation', 'c.1429C>G', 'p.L477V', 'c.(1429-1431)Cta>Gta']], 14865497: [['COLO320', 'g.chr11:14865497C>T', 'PDE3B', 'Silent', 'c.2445C>T', 'p.Y815Y', 'c.(2443-2445)taC>taT']], 17448645: [['COLO320', 'g.chr11:17448645C>G', 'ABCC8', 'Missense_Mutation', 'c.2173G>C', 'p.A725P', 'c.(2173-2175)Gcc>Ccc']], 18010266: [['COLO320', 'g.chr11:18010266T>A', 'SERGEF', 'Missense_Mutation', 'c.722A>T', 'p.H241L', 'c.(721-723)cAt>cTt']], 26695058: [['COLO320', 'g.chr11:26695058T>A', 'SLC5A12', 'Missense_Mutation', 'c.1598A>T', 'p.D533V', 'c.(1597-1599)gAt>gTt']], 28057037: [['COLO320', 'g.chr11:28057037C>A', 'KIF18A', 'Missense_Mutation', 'c.2401G>T', 'p.D801Y', 'c.(2401-2403)Gat>Tat']], 30033995: [['NCIH522', 'g.chr11:30033995C>A', 'KCNA4', 'Missense_Mutation', 'c.231G>T', 'p.Q77H', 'c.(229-231)caG>caT']], 34118200: [['COLO320', 'g.chr11:34118200G>C', 'CAPRIN1', 'Missense_Mutation', 'c.1880G>C', 'p.G627A', 'c.(1879-1881)gGc>gCc']], 34184242: [['COLO320', 'g.chr11:34184242C>A', 'ABTB2', 'Missense_Mutation', 'c.2099G>T', 'p.G700V', 'c.(2098-2100)gGc>gTc']], 35640982: [['HCC827', 'g.chr11:35640982C>T', 'FJX1', 'Silent', 'c.798C>T', 'p.R266R', 'c.(796-798)cgC>cgT']], 36597538: [['NCIH522', 'g.chr11:36597538T>C', 'RAG1', 'Missense_Mutation', 'c.2684T>C', 'p.V895A', 'c.(2683-2685)gTa>gCa']], 45827820: [['NCIH522', 'g.chr11:45827820G>T', 'SLC35C1', 'Silent', 'c.468G>T', 'p.L156L', 'c.(466-468)ctG>ctT']], 46784201: [['NCIH522', 'g.chr11:46784201delA', 'CKAP5', 'Frame_Shift_Del', 'c.4003delT', 'p.S1335fs', 'c.(4003-4005)tccfs']], 46811700: [['HCC827', 'g.chr11:46811700G>T', 'CKAP5', 'Missense_Mutation', 'c.1801C>A', 'p.P601T', 'c.(1801-1803)Ccc>Acc']], 55681588: [['HCC827', 'g.chr11:55681588C>G', 'OR5W2', 'Missense_Mutation', 'c.471G>C', 'p.L157F', 'c.(469-471)ttG>ttC']], 55682025: [['NCIH522', 'g.chr11:55682025_55682026insC', 'OR5W2', 'Frame_Shift_Ins', 'c.33_34insG', 'p.F12fs', 'c.(31-36)gattttfs']], 56230745: [['NCIH522', 'g.chr11:56230745T>C', 'OR5M9', 'Missense_Mutation', 'c.133A>G', 'p.I45V', 'c.(133-135)Att>Gtt']], 56345051: [['HCC827', 'g.chr11:56345051G>C', 'OR5M10', 'Missense_Mutation', 'c.147C>G', 'p.I49M', 'c.(145-147)atC>atG']], 57886597: [['HCC827', 'g.chr11:57886597C>A', 'OR9I1', 'Missense_Mutation', 'c.320G>T', 'p.C107F', 'c.(319-321)tGt>tTt']], 59190311: [['NCIH522', 'g.chr11:59190311G>T', 'OR5A2', 'Missense_Mutation', 'c.116C>A', 'p.T39K', 'c.(115-117)aCg>aAg']], 60235836: [['COLO320', 'g.chr11:60235836A>C', 'MS4A1', 'Missense_Mutation', 'c.789A>C', 'p.E263D', 'c.(787-789)gaA>gaC']], 60640853: [['NCIH522', 'g.chr11:60640853A>G', 'ZP1', 'Missense_Mutation', 'c.1246A>G', 'p.T416A', 'c.(1246-1248)Acc>Gcc']], 60714147: [['COLO320', 'g.chr11:60714147G>T', 'SLC15A3', 'Silent', 'c.705C>A', 'p.G235G', 'c.(703-705)ggC>ggA']], 62289230: [['HCC827', 'g.chr11:62289230G>A', 'AHNAK', 'Missense_Mutation', 'c.12659C>T', 'p.S4220L', 'c.(12658-12660)tCa>tTa']], 62462182: [['COLO320', 'g.chr11:62462182A>C', 'BSCL2', 'Splice_Site', 'c.296T>G', 'p.V99G', 'c.(295-297)gTg>gGg']], 64109520: [['HCC827', 'g.chr11:64109520C>T', 'CCDC88B', 'Missense_Mutation', 'c.730C>T', 'p.P244S', 'c.(730-732)Ccc>Tcc']], 66062444: [['HCC827', 'g.chr11:66062444T>A', 'TMEM151A', 'Missense_Mutation', 'c.727T>A', 'p.F243I', 'c.(727-729)Ttc>Atc']], 67434390: [['NCIH522', 'g.chr11:67434390C>A', 'ALDH3B2', 'Missense_Mutation', 'c.17G>T', 'p.R6L', 'c.(16-18)cGg>cTg']], 69588123: [['HCC827', 'g.chr11:69588123C>T', 'FGF4', 'Missense_Mutation', 'c.575G>A', 'p.R192Q', 'c.(574-576)cGa>cAa']], 69588191: [['COLO320', 'g.chr11:69588191G>A', 'FGF4', 'Silent', 'c.507C>T', 'p.Y169Y', 'c.(505-507)taC>taT']], 73745276: [['HCC827', 'g.chr11:73745276G>A', 'C2CD3', 'Nonsense_Mutation', 'c.5929C>T', 'p.R1977*', 'c.(5929-5931)Cga>Tga']], 74057852: [['COLO320', 'g.chr11:74057852T>C', 'PGM2L1', 'Missense_Mutation', 'c.962A>G', 'p.E321G', 'c.(961-963)gAg>gGg']], 74081972: [['NCIH522', 'g.chr11:74081972_74081973insA', 'PGM2L1', 'Frame_Shift_Ins', 'c.445_446insT', 'p.S149fs', 'c.(445-447)tcafs']], 75439853: [['NCIH522', 'g.chr11:75439853C>A', 'MOGAT2', 'Silent', 'c.669C>A', 'p.I223I', 'c.(667-669)atC>atA']], 76063306: [['COLO320', 'g.chr11:76063306G>A', 'PRKRIR', 'Silent', 'c.888C>T', 'p.F296F', 'c.(886-888)ttC>ttT']], 77066829: [['COLO320', 'g.chr11:77066829G>T', 'PAK1', 'Missense_Mutation', 'c.656C>A', 'p.T219K', 'c.(655-657)aCa>aAa']], 78282480: [['NCIH522', 'g.chr11:78282480G>A', 'NARS2', 'Missense_Mutation', 'c.151C>T', 'p.R51C', 'c.(151-153)Cgt>Tgt']], 85712078: [['NCIH522', 'g.chr11:85712078C>T', 'PICALM', 'Splice_Site', 'c.1017G>A', 'p.K339K', 'c.(1015-1017)aaG>aaA']], 89135573: [['NCIH522', 'g.chr11:89135573G>C', 'NOX4', 'Nonsense_Mutation', 'c.767C>G', 'p.S256*', 'c.(766-768)tCa>tGa']], 89608809: [['COLO320', 'g.chr11:89608809delC', 'TRIM64B', 'Frame_Shift_Del', 'c.377delG', 'p.S126fs', 'c.(376-378)agcfs']], 92539633: [['HCC827', 'g.chr11:92539633T>A', 'FAT3', 'Missense_Mutation', 'c.8749T>A', 'p.S2917T', 'c.(8749-8751)Tca>Aca']], 94134124: [['HCC827', 'g.chr11:94134124A>C', 'GPR83', 'Missense_Mutation', 'c.290T>G', 'p.F97C', 'c.(289-291)tTc>tGc']], 94320399: [['HCC827', 'g.chr11:94320399C>A', 'PIWIL4', 'Silent', 'c.900C>A', 'p.L300L', 'c.(898-900)ctC>ctA']], 99827647: [['NCIH522', 'g.chr11:99827647A>G', 'CNTN5', 'Silent', 'c.783A>G', 'p.T261T', 'c.(781-783)acA>acG']], 108178655: [['NCIH522', 'g.chr11:108178655_108178656insA', 'ATM', 'Frame_Shift_Ins', 'c.5706_5707insA', 'p.K1903fs', 'c.(5707-5709)aaafs']], 108361870: [['COLO320', 'g.chr11:108361870_108361872delACT', 'KDELC2', 'In_Frame_Del', 'c.225_227delAGT', 'p.75_76VV>V', 'c.(223-228)gtagtc>gtc']], 108544196: [['COLO320', 'g.chr11:108544196A>G', 'DDX10', 'Missense_Mutation', 'c.189A>G', 'p.I63M', 'c.(187-189)atA>atG']], 110333090: [['NCIH522', 'g.chr11:110333090_110333091insG', 'FDX1', 'Frame_Shift_Ins', 'c.453_454insG', 'p.C152fs', 'c.(454-456)tgcfs']], 116827730: [['NCIH522', 'g.chr11:116827730C>G', 'SIK3', 'Missense_Mutation', 'c.150G>C', 'p.K50N', 'c.(148-150)aaG>aaC']], 117280540: [['NCIH522', 'g.chr11:117280540G>A', 'CEP164', 'Missense_Mutation', 'c.3955G>A', 'p.A1319T', 'c.(3955-3957)Gcc>Acc']], 118773004: [['COLO320', 'g.chr11:118773004_118773005delGG', 'BCL9L', 'Frame_Shift_Del', 'c.1447_1448delCC', 'p.P483fs', 'c.(1447-1449)ccgfs']], 118899121: [['NCIH522', 'g.chr11:118899121C>A', 'SLC37A4', 'Missense_Mutation', 'c.164G>T', 'p.S55I', 'c.(163-165)aGc>aTc']], 119027723: [['NCIH522', 'g.chr11:119027723C>T', 'ABCG4', 'Splice_Site', 'c.1067C>T', 'p.P356L', 'c.(1066-1068)cCg>cTg']], 120996012: [['COLO320', 'g.chr11:120996012T>G', 'TECTA', 'Splice_Site', 'c.1205T>G', 'p.V402G', 'c.(1204-1206)gTt>gGt']], 121421347: [['COLO320', 'g.chr11:121421347T>A', 'SORL1', 'Missense_Mutation', 'c.2234T>A', 'p.L745Q', 'c.(2233-2235)cTg>cAg']], 124294185: [['NCIH522', 'g.chr11:124294185C>T', 'OR8B4', 'Missense_Mutation', 'c.583G>A', 'p.E195K', 'c.(583-585)Gag>Aag']], 124829794: [['NCIH522', 'g.chr11:124829794G>A', 'CCDC15', 'Silent', 'c.411G>A', 'p.L137L', 'c.(409-411)ttG>ttA']], 128781658: [['HCC827', 'g.chr11:128781658G>C', 'KCNJ5', 'Missense_Mutation', 'c.490G>C', 'p.G164R', 'c.(490-492)Ggg>Cgg']], 130003624: [['HCC827', 'g.chr11:130003624G>C', 'APLP2', 'Splice_Site', '', '', 'c.e12+1']], 130289174: [['NCIH522', 'g.chr11:130289174G>C', 'ADAMTS8', 'Missense_Mutation', 'c.734C>G', 'p.T245R', 'c.(733-735)aCg>aGg']], 134023206: [['NCIH522', 'g.chr11:134023206T>C', 'NCAPD3', 'Silent', 'c.4305A>G', 'p.P1435P', 'c.(4303-4305)ccA>ccG']]}}
    tmp_list = [['NCIH522', 'g.chr1:117127475_117127476insTCT', 'IGSF3', 'In_Frame_Ins', 'c.2639_2640insAGA', 'p.880_880E>EE', 'c.(2638-2640)gag>gaAGAg']
                ,['COLO320', 'g.chr22:29091840_29091841TG>CA', 'CHEK2', 'Missense_Mutation', 'c.1116_1117CA>TG', 'p.K373E', 'c.(1114-1119)tcCAag>tcTGag']
                ,['HCC827', 'g.chr10:62547988_62547989insTTA', 'CDK1', 'Splice_Site', 'c.489_489insTTA', 'p.164_165ins*', 'c.(490-492)gta>gtTTAa']
                ,['HCC827', 'g.chr7:55242466_55242480delGAATTAAGAGAAGCA', 'EGFR', 'In_Frame_Del', 'c.2101_2115delGAATTAAGAGAAGCA', 'p.ELREA701del', 'c.(2101-2115)gaattaagagaagcadel']
                ,['HCC827', 'g.chr17:7578195_7578197delCAC', 'TP53', 'In_Frame_Del', 'c.652_654delGTG', 'p.V218del', 'c.(652-654)gtgdel']
                ,['NCIH522', 'g.chr5:1878626_1878661delGCCCGGCAGTGGCTCCGGCCCGGCCGCCGCGCTGCG', 'IRX4', 'In_Frame_Del', 'c.982_1017delCGCAGCGCGGCGGCCGGGCCGGAGCCACTGCCGGGC', 'p.RSAAAGPEPLPG328del', 'c.(982-1017)cgcagcgcggcggccgggccggagccactgccgggcdel']
                ,['HCC827', 'g.chr3:150612003_150612005delAGA', 'FAM188B2', 'In_Frame_Del', 'c.85_87delTCT', 'p.S29del', 'c.(85-87)tctdel']
                ,['COLO320', 'g.chr1:14105547_14105548insC', 'PRDM2', 'Frame_Shift_Ins', 'c.1257_1258insC', 'p.P420fs', 'c.(1258-1260)cccfs']
                ,['HCC827', 'g.chr1:18554419_18554420insC', 'IGSF21', 'Frame_Shift_Ins', 'c.98_99insC', 'p.LP33fs', 'c.(97-102)ctccccfs']
                ,['COLO320', 'g.chr1:44595414_44595415insC', 'KLF17', 'Frame_Shift_Ins', 'c.471_472insC', 'p.L158fs', 'c.(472-474)ctgfs']
                ,['HCC827', 'g.chr1:222705420_222705421insT', 'HHIPL2', 'Frame_Shift_Ins', 'c.1610_1611insA', 'p.N537fs', 'c.(1609-1611)aacfs']
                ,['COLO320', 'g.chr1:237777657_237777658insA', 'RYR2', 'Frame_Shift_Ins', 'c.5229_5230insA', 'p.K1744fs', 'c.(5230-5232)aaafs']
                ,['COLO320', 'g.chr12:27066991_27066992insT', 'ASUN', 'Frame_Shift_Ins', 'c.1489_1490insA', 'p.T497fs', 'c.(1489-1491)acafs']
                ,['NCIH522', 'g.chr12:40748133_40748145delCTTGGTGCATCTT', 'LRRK2', 'Frame_Shift_Del', 'c.6609_6621delCTTGGTGCATCTT', 'p.ALVHL2203fs', 'c.(6607-6621)gccttggtgcatcttfs']
                ,['NCIH522', 'g.chr11:55682025_55682026insC', 'OR5W2', 'Frame_Shift_Ins', 'c.33_34insG', 'p.F12fs', 'c.(31-36)gattttfs']
                ,['NCIH522', 'g.chr11:74081972_74081973insA', 'PGM2L1', 'Frame_Shift_Ins', 'c.445_446insT', 'p.S149fs', 'c.(445-447)tcafs']
                ,['NCIH522', 'g.chr11:108178655_108178656insA', 'ATM', 'Frame_Shift_Ins', 'c.5706_5707insA', 'p.K1903fs', 'c.(5707-5709)aaafs']
                ,['NCIH522', 'g.chr11:110333090_110333091insG', 'FDX1', 'Frame_Shift_Ins', 'c.453_454insG', 'p.C152fs', 'c.(454-456)tgcfs']
                ,['COLO320', 'g.chr11:118773004_118773005delGG', 'BCL9L', 'Frame_Shift_Del', 'c.1447_1448delCC', 'p.P483fs', 'c.(1447-1449)ccgfs']
                ,['COLO320', 'g.chr2:43458193_43458194insC', 'THADA', 'Frame_Shift_Ins', 'c.5755_5756insG', 'p.A1919fs', 'c.(5755-5757)gccfs']
                ,['COLO320', 'g.chr2:74719525_74719526insG', 'TTC31', 'Frame_Shift_Ins', 'c.1114_1115insG', 'p.W372fs', 'c.(1114-1116)tggfs']
                ,['COLO320', 'g.chr2:128096510_128096511insT', 'MAP3K2', 'Frame_Shift_Ins', 'c.120_121insA', 'p.Q41fs', 'c.(118-123)aaacagfs']
                ,['NCIH522', 'g.chr2:160261539_160261540insT', 'BAZ2B', 'Frame_Shift_Ins', 'c.2763_2764insA', 'p.Q922fs', 'c.(2761-2766)aaacaafs']
                ,['HCC827', 'g.chr2:160667052_160667053insT', 'LY75', 'Frame_Shift_Ins', 'c.4683_4684insA', 'p.L1562fs', 'c.(4681-4686)aaattgfs']
                ,['NCIH522', 'g.chr2:190554522_190554528delAGATTAT', 'ANKAR', 'Frame_Shift_Del', 'c.871_877delAGATTAT', 'p.RLY291fs', 'c.(871-879)agattatatfs']
                ,['HCC827', 'g.chr2:231266128_231266129insT', 'SP140L', 'Frame_Shift_Ins', 'c.1470_1471insT', 'p.F491fs', 'c.(1471-1473)tttfs']
                ,['NCIH522', 'g.chrX:54949009_54949010insG', 'TRO', 'Splice_Site', '', '', 'c.e3-1']
                ,['COLO320', 'g.chr6:71234304_71234305insA', 'FAM135A', 'Frame_Shift_Ins', 'c.1517_1518insA', 'p.TK506fs', 'c.(1516-1521)acaaaafs']
                ,['COLO320', 'g.chr6:106554865_106554866insC', 'PRDM1', 'Frame_Shift_Ins', 'c.1982_1983insC', 'p.K662fs', 'c.(1981-1986)tgcaagfs']
                ,['COLO320', 'g.chr18:30791876_30791877insT', 'CCDC178', 'Frame_Shift_Ins', 'c.2221_2222insA', 'p.T741fs', 'c.(2221-2223)actfs']
                ,['HCC827', 'g.chr7:55242466_55242480delGAATTAAGAGAAGCA', 'EGFR', 'In_Frame_Del', 'c.2101_2115delGAATTAAGAGAAGCA', 'p.ELREA701del', 'c.(2101-2115)gaattaagagaagcadel']
                ,['COLO320', 'g.chr7:131195806_131195807insG', 'PODXL', 'Frame_Shift_Ins', 'c.486_487insC', 'p.S163fs', 'c.(484-489)agcagcfs']
                ,['COLO320', 'g.chr15:43814912_43814913delAA', 'MAP1A', 'Frame_Shift_Del', 'c.1241_1242delAA', 'p.E414fs', 'c.(1240-1242)gaafs']
                ,['NCIH522', 'g.chr17:1373593_1373594insG', 'MYO1C', 'Frame_Shift_Ins', 'c.2296_2297insC', 'p.R766fs', 'c.(2296-2298)cgcfs']
                ,['HCC827', 'g.chr17:4937386_4937389delAGGA', 'SLC52A1', 'Frame_Shift_Del', 'c.395_398delTCCT', 'p.FL132fs', 'c.(394-399)ttcctgfs']
                ,['HCC827', 'g.chr17:7578195_7578197delCAC', 'TP53', 'In_Frame_Del', 'c.652_654delGTG', 'p.V218del', 'c.(652-654)gtgdel']
                ,['COLO320', 'g.chr10:28233229_28233230insG', 'ARMC4', 'Frame_Shift_Ins', 'c.1664_1665insC', 'p.A555fs', 'c.(1663-1665)gcafs']
                ,['NCIH522', 'g.chr8:82583189_82583190insA', 'IMPA1', 'Frame_Shift_Ins', 'c.550_551insT', 'p.C184fs', 'c.(550-552)tgcfs']
                ,['HCC827', 'g.chr8:96275920_96275921insT', 'C8orf37', 'Frame_Shift_Ins', 'c.237_238insA', 'p.P80fs', 'c.(235-240)aaacccfs']
                ,['HCC827', 'g.chr8:105257154_105257155insA', 'RIMS2', 'Frame_Shift_Ins', 'c.3399_3400insA', 'p.K1134fs', 'c.(3400-3402)aaafs']
                ,['NCIH522', 'g.chr8:133150232_133150233insT', 'KCNQ3', 'Frame_Shift_Ins', 'c.1599_1600insA', 'p.F534fs', 'c.(1597-1602)aaattcfs']
                ,['HCC827', 'g.chr8:144921484_144921485insCCGGAG', 'NRBP2', 'Splice_Site', '', '', 'c.e6+1']
                ,['COLO320', 'g.chr9:23701461_23701462insT', 'ELAVL2', 'Frame_Shift_Ins', 'c.628_629insA', 'p.T210fs', 'c.(628-630)accfs']
                ,['COLO320', 'g.chr9:95609089_95609095delTCCAGTG', 'ZNF484', 'Frame_Shift_Del', 'c.1974_1980delCACTGGA', 'p.HTG658fs', 'c.(1972-1980)cacactggafs']
                ,['HCC827', 'g.chr16:773898_773899insT', 'CCDC78', 'Frame_Shift_Ins', 'c.1091_1092insA', 'p.K364fs', 'c.(1090-1092)aagfs']
                ,['COLO320', 'g.chr16:2506811_2506812insT', 'CCNF', 'Frame_Shift_Ins', 'c.2151_2152insT', 'p.L718fs', 'c.(2152-2154)ttgfs']
                ,['NCIH522', 'g.chr16:11362767_11362768insTT', 'TNP2', 'Frame_Shift_Ins', 'c.352_353insAA', 'p.M118fs', 'c.(352-354)atgfs']
                ,['HCC827', 'g.chr16:20996843_20996858delGACCAGGGACATGGGG', 'DNAH3', 'Frame_Shift_Del', 'c.7206_7221delCCCCATGTCCCTGGTC', 'p.APMSLV2402fs', 'c.(7204-7221)gcccccatgtccctggtcfs']
                ,['COLO320', 'g.chr13:32811598_32811599insC', 'FRY', 'Frame_Shift_Ins', 'c.5893_5894insC', 'p.T1965fs', 'c.(5893-5895)accfs']
                ,['NCIH522', 'g.chr5:1878626_1878661delGCCCGGCAGTGGCTCCGGCCCGGCCGCCGCGCTGCG', 'IRX4', 'In_Frame_Del', 'c.982_1017delCGCAGCGCGGCGGCCGGGCCGGAGCCACTGCCGGGC', 'p.RSAAAGPEPLPG328del', 'c.(982-1017)cgcagcgcggcggccgggccggagccactgccgggcdel']
                ,['HCC827', 'g.chr5:37701207_37701208insT', 'WDR70', 'Frame_Shift_Ins', 'c.1240_1241insT', 'p.L414fs', 'c.(1240-1242)cttfs']
                ,['HCC827', 'g.chr3:124538604_124538605insT', 'ITGB5', 'Frame_Shift_Ins', 'c.1019_1020insA', 'p.N340fs', 'c.(1018-1020)aacfs']
                ,['HCC827', 'g.chr3:150612003_150612005delAGA', 'FAM188B2', 'In_Frame_Del', 'c.85_87delTCT', 'p.S29del', 'c.(85-87)tctdel']
                ,['NCIH522', 'g.chr19:11097625_11097626delCC', 'SMARCA4', 'Frame_Shift_Del', 'c.805_806delCC', 'p.P270fs', 'c.(805-807)cccfs']
                ,['COLO320', 'g.chr19:20230027_20230031delATAAG', 'ZNF90', 'Frame_Shift_Del', 'c.1664_1668delATAAG', 'p.HK555fs', 'c.(1663-1668)cataagfs']
                ,['COLO320', 'g.chr19:38380519_38380520insC', 'WDR87', 'Frame_Shift_Ins', 'c.3674_3675insG', 'p.G1225fs', 'c.(3673-3675)ggafs']
                ,['COLO320', 'g.chr19:39760415_39760416insC', 'IFNL2', 'Frame_Shift_Ins', 'c.458_459insC', 'p.R154fs', 'c.(457-462)ggccgcfs']
                ,['NCIH522', 'g.chr19:46201958_46201959insT', 'QPCTL', 'Splice_Site', '', '', 'c.e4+1']
                ,['COLO320', 'g.chr4:27010053_27010054insA', 'STIM2', 'Frame_Shift_Ins', 'c.1153_1154insA', 'p.E385fs', 'c.(1153-1155)gaafs']
                ,['NCIH522', 'g.chr4:57179502_57179503insA', 'KIAA1211', 'Frame_Shift_Ins', 'c.494_495insA', 'p.PK165fs', 'c.(493-498)ccaaaafs']
                ,['HCC827', 'g.chr4:88261656_88261657insA', 'HSD17B11', 'Frame_Shift_Ins', 'c.797_798insT', 'p.L266fs', 'c.(796-798)ttafs']
                ,['COLO320', 'g.chr4:126370100_126370101insG', 'FAT4', 'Frame_Shift_Ins', 'c.7929_7930insG', 'p.G2644fs', 'c.(7930-7932)ggtfs']
                ,['HCC827', 'g.chrM:3902_3906delACCTT', 'MT-ND1', 'Frame_Shift_Del', 'c.596_600delACCTT', 'p.DL199fs', 'c.(595-600)gaccttfs']
                ,['HCC827', 'g.chrM:3908_3909insAAGGT', 'MT-ND1', 'Frame_Shift_Ins', 'c.602_603insAAGGT', 'p.E202fs', 'c.(601-606)gccgaafs']
                ,['HCC827', 'g.chrM:12417_12418insA', 'MT-ND5', 'Frame_Shift_Ins', 'c.81_82insA', 'p.K28fs', 'c.(82-84)aaafs']
                ,['HCC827', 'g.chr1:44445557C>A', 'B4GALT2', 'De_novo_Start_OutOfFrame', '', '', '']
                ,['HCC827', 'g.chr1:74819791delT', 'FPGT-TNNI3K', 'Frame_Shift_Del', 'c.1458delT', 'p.C486fs', 'c.(1456-1458)tgtfs']
                ,['NCIH522', 'g.chr1:152191016delC', 'HRNR', 'Frame_Shift_Del', 'c.3089delG', 'p.G1030fs', 'c.(3088-3090)ggcfs']
                ,['NCIH522', 'g.chr1:165766945G>T', 'TMCO1', 'De_novo_Start_OutOfFrame', '', '', '']
                ,['NCIH522', 'g.chr12:53937233C>A', 'ATF7', 'Splice_Site', '', '', 'c.e4-1']
                ,['NCIH522', 'g.chr12:101869220G>A', 'SPIC', 'De_novo_Start_OutOfFrame', '', '', '']
                ,['NCIH522', 'g.chr11:46784201delA', 'CKAP5', 'Frame_Shift_Del', 'c.4003delT', 'p.S1335fs', 'c.(4003-4005)tccfs']
                ,['COLO320', 'g.chr11:89608809delC', 'TRIM64B', 'Frame_Shift_Del', 'c.377delG', 'p.S126fs', 'c.(376-378)agcfs']
                ,['HCC827', 'g.chr11:130003624G>C', 'APLP2', 'Splice_Site', '', '', 'c.e12+1']
                ,['NCIH522', 'g.chr2:171673225T>G', 'GAD1', 'De_novo_Start_OutOfFrame', '', '', '']
                ,['NCIH522', 'g.chr2:197673937C>A', 'C2orf66', 'Splice_Site', '', '', 'c.e1+1']
                ,['NCIH522', 'g.chr14:93685667delT', 'UBR7', 'Frame_Shift_Del', 'c.920delT', 'p.L307fs', 'c.(919-921)ctgfs']
                ,['COLO320', 'g.chr6:24418688G>C', 'MRS2', 'Splice_Site', '', '', 'c.e9-1']
                ,['COLO320', 'g.chr7:48336827A>T', 'ABCA13', 'Splice_Site', '', '', 'c.e22-1']
                ,['NCIH522', 'g.chr7:117880053G>A', 'ANKRD7', 'Splice_Site', '', '', 'c.e6+1']
                ,['NCIH522', 'g.chr7:150656676delG', 'KCNH2', 'Frame_Shift_Del', 'c.456delC', 'p.T152fs', 'c.(454-456)accfs']
                ,['NCIH522', 'g.chr15:72511412delC', 'PKM', 'Frame_Shift_Del', 'c.27delG', 'p.G9fs', 'c.(25-27)gggfs']
                ,['NCIH522', 'g.chr17:7578277delG', 'TP53', 'Frame_Shift_Del', 'c.572delC', 'p.P191fs', 'c.(571-573)cctfs']
                ,['COLO320', 'g.chr17:7726903delC', 'DNAH2', 'Frame_Shift_Del', 'c.11286delC', 'p.H3762fs', 'c.(11284-11286)cacfs']
                ,['HCC827', 'g.chr17:18218206G>A', 'TOP3A', 'De_novo_Start_OutOfFrame', '', '', '']
                ,['COLO320', 'g.chr17:40346443C>T', 'GHDC', 'De_novo_Start_OutOfFrame', '', '', '']
                ,['HCC827', 'g.chr10:50507279G>A', 'C10orf71', 'Splice_Site', '', '', 'c.e1+1']
                ,['NCIH522', 'g.chr10:118639498G>A', 'ENO4', 'Splice_Site', '', '', 'c.e6+1']
                ,['NCIH522', 'g.chr10:124711543C>G', 'C10orf88', 'Splice_Site', '', '', 'c.e3-1']
                ,['NCIH522', 'g.chr20:57226968A>T', 'STX16', 'De_novo_Start_OutOfFrame', '', '', '']
                ,['NCIH522', 'g.chr8:23294446C>T', 'ENTPD4', 'Splice_Site', '', '', 'c.e10+1']
                ,['HCC827', 'g.chr16:31084206C>T', 'ZNF668', 'Splice_Site', '', '', 'c.e1+1']
                ,['NCIH522', 'g.chr13:51943176delG', 'INTS6', 'Frame_Shift_Del', 'c.2375delC', 'p.S792fs', 'c.(2374-2376)tcafs']
                ,['COLO320', 'g.chr5:102342716G>T', 'PAM', 'Splice_Site', '', '', 'c.e18+1']
                ,['COLO320', 'g.chr5:147261212G>T', 'SCGB3A2', 'Splice_Site', '', '', 'c.e2+1']
                ,['COLO320', 'g.chr5:149901062delC', 'NDST1', 'Frame_Shift_Del', 'c.246delC', 'p.D82fs', 'c.(244-246)gacfs']
                ,['HCC827', 'g.chr19:4702614C>A', 'DPP9', 'Splice_Site', '', '', 'c.e8+1']
                ,['NCIH522', 'g.chr4:25673222G>A', 'SLC34A2', 'Splice_Site', '', '', 'c.e9-1']
                ,['HCC827', 'g.chr21:28315697A>T', 'ADAMTS5', 'Splice_Site', '', '', 'c.e3+1']
                ]
    # tmp_dict = {}
    # for val_list in tmp_list:
    #     tmp_dict.update({val_list[0] + "^" + val_list[1]: val_list})
    # print(tmp_dict)
    result_dict = logic.get_seqs_bfr_aft_trgt_idx(sorted_mut_dict, INITIAL_MAIN)

    # for chr_key, val_dict in result_dict.items():
    #     print(chr_key)
    #     for seq_key, val in val_dict.items():
    #         print(seq_key)
    #         print(val)

    util.make_excel(WORK_DIR + "analyze_hg19", sorted_mut_dict, result_dict)



def test():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    mut_sum_dict = util.read_txt_dvd_by_tab(WORK_DIR + MUT_FILE)
    sorted_mut_dict = logic_prep.sort_dict(mut_sum_dict)

    tmp_list_a = []
    tmp_list_b = []
    for chr_key, val_list in sorted_mut_dict.items():
        # print("chr_key : " + chr_key)
        for seq_key, tmp_val_list in val_list.items():
            # print("seq_key : " + str(seq_key))
            for tmp_val in tmp_val_list:
                Cellline = tmp_val[0]
                Genome_Change = tmp_val[1]
                tmp_key = Cellline + "^" + Genome_Change  # HCC827^g.chr1:915845C>A
                Genomesequence = tmp_val[3]
                cDNA_change = tmp_val[4]
                Codon_Change = tmp_val[6]  # Codon_Change : c.(223-225)Ggt>Tgt
                Codon_Change_arr = Codon_Change.split(">")  # ['c.(223-225)Ggt', 'Tgt']

                # if Genome_Change[-1:] == cDNA_change[-1:]:
                #     print("+")
                # elif cDNA_change == "":
                #     print(Genome_Change[-1:] + " : " + cDNA_change)

                if len(Codon_Change_arr) > 1:
                    pass
                elif "_" in Genome_Change:
                    end_seq = re.findall(r'\d+',Genome_Change[Genome_Change.index(":"):])
                    # if len(Codon_Change_arr) > 1:
                    #     if len(Codon_Change_arr[1]) != 3:
                    #         # print("################## mut_seq's len is not 3 ############################")
                    #         tmp_list_a.append(tmp_val)
                    #     #     # print("[" + tmp_key + "] ")
                    #     #     print(tmp_val)
                    #     #     # print("")
                    # else:
                    #     # print("################## deli with '_' #######################################")
                    #     if 'In_Frame' in Genomesequence:
                    #         in_frm_seq = cDNA_change[cDNA_change.index("del") + len('del'):]
                    #         print("in_frm_seq" + in_frm_seq)
                    #         if in_frm_seq.upper() in tmp_key:
                    #             print("+ : " + end_seq[1])
                    #             print(type(end_seq[1]))
                    #
                    #         else:
                    #             print("- : " + str(end_seq[1]))
                    #         print(tmp_val)
                    #         tmp_list_b.append(tmp_val)
                else:
                    print(tmp_val)
    # print("################## mut_seq's len is not 3 ############################")
    # print(len(tmp_list_a))
    # for arr_a in tmp_list_a:
    #     print(arr_a)
    # print("")
    # print("################## deli with '_' #######################################")
    # print(len(tmp_list_b))
    # for arr_b in tmp_list_b:
    #     print(arr_b)





start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# test()
main_YG()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))


