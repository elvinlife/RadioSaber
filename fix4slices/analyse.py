from collections import defaultdict
import sys
n_users = 0

def get_exp(fname):
    begin_ts = 9000
    end_ts = 10000
    cumu_bytes = [0 for i in range(n_users)]
    cumu_rbs = [0 for i in range(n_users)]

    with open(fname, "r") as fin:
        for line in fin:
            words = line.split(" ")
            if not words[0].isdigit():
                continue
            if int(words[0]) > begin_ts:
                flow = int(words[2])
                cumu_rbs[flow] = int( words[6] )
                cumu_bytes[flow] = int( words[4] )
            if int(words[0]) > end_ts:
                break
    return cumu_bytes, cumu_rbs

for ue in range(2, 10, 2):
    n_users = ue * 10
    dname = str(ue) + "ues-40ms/"
    print(dname)
    for i in range(3):
        byte, rbs = get_exp(dname + "nvs_pf" + str(i) + ".log")
        print("bytes: sum: {}\t rbs: sum: {}".format( sum(byte)/1000, sum(rbs)/1000 ))
        byte, rbs = get_exp(dname + "oracle_pf" + str(i) + ".log")
        print("bytes: sum: {}\t rbs: sum: {}".format( sum(byte)/1000, sum(rbs)/1000 ))
        byte, rbs = get_exp(dname + "single_pf" + str(i) + ".log")
        print("bytes: sum: {}\t rbs: sum: {}".format( sum(byte)/1000, sum(rbs)/1000 ))
