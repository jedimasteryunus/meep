import os

os.system("h5topng -t 0:332 -R -Zc dkbluered -a yarg -A cladded_basic-eps-000000.00.h5 cladded_basic-ez.h5")
os.system("convert cladded_basic-ez.t*.png cladded_basic-ez.gif")
os.system("eog -f cladded_basic-ez.gif")
