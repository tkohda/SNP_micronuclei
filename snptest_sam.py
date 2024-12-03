import sys
import re

fname = sys.argv[1].split('.')
i = 0
start = 0
bd = 0
c3h = 0
er = 0
position = []
ref_b = []
alt_b = []

with open(sys.argv[1], 'r') as sam, open(sys.argv[2], 'r') as snpdb, \
     open(f'./{fname[1]}_result.txt', 'w') as resultfile, \
     open(f'./{fname[1]}_result2.txt', 'w') as resultfile2:

    for line2 in snpdb:
        list2 = line2.split('\t')
        position.append(int(list2[1]))
        ref_b.append(list2[3])
        alt_b.append(list2[4])
        i += 1

    snp = 0
    base = "n"
    strain = ""

    for line in sam:
        line_list = line.split()
        seq = line_list[9]
        cigar = line_list[5]
        cigar_num = re.findall(r'\d+', cigar)
        cigar_num_i = [int(s) for s in cigar_num]
        cigar_str = re.findall(r'[^\d]+', cigar)
        b = list(seq)
        seq_len = len(b)
        start_b = int(line_list[3])
        for j in range(start, i):
            if (start_b < position[j] + 1) and (int(line_list[3]) > position[j] - seq_len):
                snp = position[j] - start_b
                snpposition = position[j]
                cigar_len = 0
                m = 0
                for k in cigar_num:
                    cigar_len = cigar_len + cigar_num_i[m]
                    if snp > cigar_len:
                        if cigar_str[m] == "D":
                            snp = snp - cigar_num_i[m]
                        if cigar_str[m] == "I":
                            snp = snp + cigar_num_i[m]
                    m += 1
                if snp > seq_len:
                    break
                base = b[snp]
                if cigar == "*":
                    break
                if base == alt_b[j]:
                    strain = "C3H"
                    c3h += 1
                else:
                    strain = "error"
                    er += 1
                if base == ref_b[j]:
                    strain = "BD"
                    bd += 1
                    er -= 1
                resultfile.write(f"{snpposition}\t{base}\t{ref_b[j]}\t{alt_b[j]}\t{strain}\t{snp}\t{cigar}\n")
                start = j - 1
#                j = i
                if (start_b > position[j] + 1):
                    break
    resultfile2.write(f"BD={bd}\nC3H={c3h}\nerror={er}\n")