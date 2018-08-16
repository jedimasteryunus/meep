import os

os.system("h5topng -t 0:332 -R -Zc dkbluered -a yarg -A grating_validation-eps-000000.00.h5 grating_validation-ez.h5")
os.system("convert -reverse grating_validation-ez.t*.png grating_validation-ez.gif")
os.system("eog -f grating_validation-ez.gif")
