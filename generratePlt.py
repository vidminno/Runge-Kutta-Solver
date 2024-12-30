""" 

Copyright 2024 Simon Peter

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import numpy as np
import pandas as pd
from matplotlib import rc
import matplotlib.pyplot as plt

if __name__ == "__main__":
    csv_file = r".\wikiPlt.csv"

    # First entry is the time
    df = pd.read_csv(csv_file, header=0)
    a = df.columns

    time = df.iloc[:, 0].to_numpy()
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, df.iloc[:, 1].to_numpy(), label=df.columns[1], color='r',  marker='.')
    plt.plot(time, df.iloc[:, 2].to_numpy(), label=df.columns[2], color='g',  marker='.')
    plt.plot(time, df.iloc[:, 3].to_numpy(), label=df.columns[3], color='b',  marker='.')
    plt.title(r"y' = sin(t)^2 * y")
    plt.ylabel(r"y")
    plt.xlabel("t")
    plt.legend(loc = "upper left")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    i= 5