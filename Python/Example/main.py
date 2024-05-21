import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import Config.constants as const
from Config import parser
from Config.gravityCoeffs import Cnm, Snm

def main():
    #print(const.MU_EARTH)

    parser.parseEOPdata('../Config/finals.all.csv', '../Config/EOP_All_data.txt')


if __name__ == '__main__':
    main()