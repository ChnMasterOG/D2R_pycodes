# a function file
import random

'''
d2r在单个短序列测序中的计算
K_mer_SIZE：K-mer大小，整型必填
S_length：序列长度，整型随机序列生成必填
SEQUENCE：基因序列，字符串，不填根据F_Vector随机生成
F_Vector：碱基出现概率，长度为4的字典，不填默认1/4
output：是否输出待测序列
'''
def D2R_in_single_sequence(K_mer_SIZE, S_length=0, SEQUENCE=None, F_Vector={'A':0.25, 'T':0.25, 'G':0.25, 'C':0.25}, output=True):
    # 生成序列
    if SEQUENCE == None:
        if S_length == 0: return -1
        SEQUENCE = ''
        for i in range(S_length):
            random_num = random.random()
            for j in list(F_Vector.keys()):
                if F_Vector[j] < random_num: random_num -= F_Vector[j]
                else: break
            SEQUENCE += j
    else: S_length = len(SEQUENCE)

    if K_mer_SIZE > S_length: return -1
    if output == True: print('待测序列：' + SEQUENCE)

    # D2R计算
    HashCount = {}
    D2R = 0
    n = S_length - K_mer_SIZE + 1
    for i in range(n):
        k_mer = SEQUENCE[i:i+K_mer_SIZE]
        if k_mer not in HashCount: HashCount[k_mer] = 1
        else: HashCount[k_mer] += 1
    for i in list(HashCount.keys()):
        RawCount = HashCount[i] * (HashCount[i] - 1)
        Uw = 1
        for j in i: Uw *= F_Vector[j]
        Mean = n**2 * Uw**2
        MeanAdjustCount = RawCount - Mean
        D2R += MeanAdjustCount
    D2R = D2R / (n * (n - 1))
    return D2R, HashCount

'''
d2r在测定长序列中的每个局部区域的重复性
K_mer_SIZE：K-mer大小，整型必填
L_length：滑动窗口序列长度，整型必填
G_length：序列长度，整型随机序列生成必填
SEQUENCE：基因序列，字符串，不填根据F_Vector随机生成
F_Vector：碱基出现概率，长度为4的字典，不填默认1/4
'''
def D2R_in_consecutive_sliding_windows_on_a_genome(K_mer_SIZE, L_length, G_length=0, SEQUENCE=None, F_Vector={'A':0.25, 'T':0.25, 'G':0.25, 'C':0.25}, output=True):
    # 生成序列
    if SEQUENCE == None:
        if G_length == 0: return -1
        SEQUENCE = ''
        for i in range(G_length):
            random_num = random.random()
            for j in list(F_Vector.keys()):
                if F_Vector[j] < random_num: random_num -= F_Vector[j]
                else: break
            SEQUENCE += j
    else: G_length = len(SEQUENCE)

    if G_length < L_length: return -1  # 滑动窗口长度应当小于序列长度
    if output == True: print('待测序列：' + SEQUENCE)

    # D2R计算
    D2R_list = []
    m = G_length - L_length + 1
    n = L_length - K_mer_SIZE + 1
    D2R_1, HashCount = D2R_in_single_sequence(K_mer_SIZE, L_length, SEQUENCE[:L_length], F_Vector, False)
    D2R_list.append(D2R_1)
    for i in range(1, m):
        prefix_k_mer = SEQUENCE[i-1:i-1+K_mer_SIZE]
        suffix_k_mer = SEQUENCE[i-K_mer_SIZE+L_length:i+L_length]
        if prefix_k_mer == suffix_k_mer:
            D2R_list.append(D2R_1)
        else:
            if HashCount[prefix_k_mer] > 1:
                deta_w1 = -2*(HashCount[prefix_k_mer]-1)/(n*(n-1))
                HashCount[prefix_k_mer] -= 1
            else:
                Uw = 1
                for j in prefix_k_mer: Uw *= F_Vector[j]
                deta_w1 = n / (n-1) * Uw ** 2
                del HashCount[prefix_k_mer]
            if suffix_k_mer in HashCount:
                deta_w2 = 2*(HashCount[suffix_k_mer])/(n*(n-1))
                HashCount[suffix_k_mer] += 1
            else:
                Uw = 1
                for j in suffix_k_mer: Uw *= F_Vector[j]
                deta_w2 = -n / (n-1) * Uw ** 2
                HashCount[suffix_k_mer] = 1
            D2R_list.append(D2R_1 + deta_w1 + deta_w2)
    return D2R_list

