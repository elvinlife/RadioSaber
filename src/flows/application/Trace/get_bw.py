def get_bw(fname):
    with open(fname + ".dat", "r") as fin:
        cumu_bytes = 0
        for line in fin:
            words = line.split("\t")
            ts = int(words[2])
            if ts % 1000 == 0:
                print("%d bitrate: %d kbps" % (ts, 8 * cumu_bytes / 1000) )
                cumu_bytes = int(words[3])
            else:
                cumu_bytes += int(words[3])

def get_cumu(fname):
    with open(fname + ".dat", "r") as fin:
        cumu_bytes = 0
        for line in fin:
            words = line.split("\t")
            ts = int(words[2])
            if ts > 10000:
                break
            cumu_bytes += int(words[3])
    print(cumu_bytes)


get_bw("foreman_H264_880k")
#get_cumu("foreman_H264_1280k")
