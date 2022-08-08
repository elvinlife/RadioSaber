import random

all_trace = 158
max_ues = 500
map_array = [ i for i in range(all_trace) ] + [i for i in range(all_trace)] + [i for i in range(all_trace) ]
random.shuffle( map_array )
for i, v in enumerate( map_array ):
    print("%d %d" % (i, v))
