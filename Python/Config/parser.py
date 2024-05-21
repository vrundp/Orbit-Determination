import numpy as np
import pandas as pd

def parseGravityCoeffs(filepath, Nmax):

    Cnm = np.zeros((Nmax + 1, Nmax + 1))
    Snm = np.zeros((Nmax + 1, Nmax + 1))

    data = pd.read_csv(filepath, header = None, sep = '\s+')

    max_row_idx = data[(data[0] == Nmax) & (data[1] == Nmax)].index

    for idx, row in data.iterrows():

        n = int(row[0])
        m = int(row[1])

        if idx <= max_row_idx:

            Cnm[n, m] = row[2]
            Snm[n, m] = row[3]

        else:

            break

    return Cnm, Snm


# only need to run once
def exportGravityCoeffs(Cnm, Snm, filename):

    with open(filename, 'w') as file:

        file.write('import numpy as np\n\n')

        file.write('Cnm = np.array([\n')
        for row in Cnm.tolist():
            file.write(f'\t{repr(row)},\n')
        file.write('])\n\n')

        file.write('Snm = np.array([\n')
        for row in Snm.tolist():
            file.write(f'\t{repr(row)},\n')
        file.write('])\n\n')


def parseEOPdata(filepath1, filepath2):
    
    data1 = pd.read_csv(filepath1, sep = ';')
    print(data1.iloc[0])
    
