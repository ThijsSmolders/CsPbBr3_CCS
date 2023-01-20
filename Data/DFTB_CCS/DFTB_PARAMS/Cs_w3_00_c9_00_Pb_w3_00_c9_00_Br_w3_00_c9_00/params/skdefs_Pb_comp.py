atomconfigs = {
    "Pb": AtomConfig(
        znuc=82,
        mass=207.2,
        occupations=[[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],  # s
                     [[3.0, 3.0], [3.0, 3.0],[3.0, 3.0],[3.0, 3.0], [1.0, 1.0]],  # p
                     [[5.0, 5.0], [5.0, 5.0], [5.0, 5.0]],  # d
                     [[7.0, 7.0]],   # f
                    ],
        valenceqns=[[6, ],  # s
                    [6, ],  # p
                    [], 
                    []
                   ],
        relativistic=True 
        ),
	}

skbases = {
   "Pb": SlaterBasis(
        exponents=[[0.83968, 2.0992, 5.248, 13.12, 32.8, 82.], # s
                   [0.83968, 2.0992, 5.248, 13.12, 32.8, 82.], # p
                   [0.83968, 2.0992, 5.248, 13.12, 32.8, 82.], # d
                   [0.83968, 2.0992, 5.248, 13.12, 32.8, 82.], # f
                   ],
        maxpowers=[ 3,  # s
                    3,  # p
                    3,  # d 
                    3,  # f 
                    ], 
        ), 
	}
compressions = {
	 "Pb": Compression(
		potcomp='potential',
		potcomp_parameters=[(2, 9.0), (2, 9.0), (2, 9.0), (2, 9.0)],
		wavecomp='potential',
		wavecomp_parameters=[(2, 3.0), (2, 3.0), (2, 3.0), (2, 3.0)],
		)
	}