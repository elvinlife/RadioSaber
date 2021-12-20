def fps_time_20(fname, ofname):
    with open(fname+".dat", "r") as fin:
        lines = fin.readlines()
    with open(ofname+".dat", "w") as fout:
        for line in lines:
            words = line.split("\t")
            words[3] = str( int(words[3]) * 10)
            fout.write("%s\n" % ("\t".join(words) ) )

fps_time_20("foreman_H264_128k", "foreman_H264_1280k")
