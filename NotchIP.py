import os

os.system("h5topng -t 0:332 -R -Zc dkbluered -a yarg -A notch-eps-000000.00.h5 notch-ez.h5")
os.system("convert -reverse notch-ez.t*.png notch-ez.gif")
os.system("eog -f notch-ez.gif")
