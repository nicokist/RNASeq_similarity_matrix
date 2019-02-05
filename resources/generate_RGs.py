#!/usr/bin/python
import sys


for line in sys.stdin:
    filename=line.split('/')[-1]
    RG_name=filename.split('.bam')[0]
    SM_name=RG_name ## used to split on .sorted.bam, maybe should be changed back
    print('@RG\tID:%s\tSM:%s' % (RG_name, SM_name))

