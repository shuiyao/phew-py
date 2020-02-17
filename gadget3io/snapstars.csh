#!/bin/tcsh

# Generate .stars file from TIPSY binaries and .grp, .sogrp files
#./allstars p50n288zw 58
#./allstars p50n288sw 58
#./allstars p50n288beta2 58
# ./allstars p50n288ezwc 58
#./allstars p50n288dsw 58 all

# Generate .starinfo file from the .stars file and SFRINFO files
# python allstars.py p50n288zw 058
# python allstars.py p50n288sw 058
# python allstars.py p50n288ezwc 058
# python allstars.py p50n288beta2 058
python allstars.py p50n288dsw 058

# ./allstars p50n288sw 108
# python stars_by_mvir.py p50n288sw 108 mh11
# python stars_by_mvir.py p50n288sw 108 mh12
# python stars_by_mvir.py p50n288sw 108 mh13

# ./allstars p50n288zw 108
# python stars_by_mvir.py p50n288zw 108 mh11
# python stars_by_mvir.py p50n288zw 108 mh12
# python stars_by_mvir.py p50n288zw 108 mh13

#./allstars p50n288beta2 108
# python stars_by_mvir.py p50n288beta2 108 mh11
# python stars_by_mvir.py p50n288beta2 108 mh12
# python stars_by_mvir.py p50n288beta2 108 mh13
