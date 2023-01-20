atomconfigs = {
    "Cs": AtomConfig(
        znuc=55,
        mass=132.90545196,
        occupations=[[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 0.0]],  # s
                     [[3.0, 3.0], [3.0, 3.0],[3.0, 3.0],[3.0, 3.0], [0.0, 0.0]],  # p
                     [[5.0, 5.0], [5.0, 5.0]],  # d
                    ],
        valenceqns=[[6, ],  # s
                    [],  # p
                    [], 
                   ],
        relativistic=True 
        ),
	}

skbases = {
   "Cs": SlaterBasis(
        exponents=[[0.5632 , 1.408, 3.52, 8.8, 22., 55.], # s
                   [0.5632 , 1.408, 3.52, 8.8, 22., 55.], # p
                   [0.5632 , 1.408, 3.52, 8.8, 22., 55.], # d 
                  ],
        maxpowers=[ 3,  # s
                    3,  # p
                    3,  # d 
                    ], 
        ), 
	}
compressions = {
	 "Cs": Compression(
		potcomp='potential',
		potcomp_parameters=[(2, 9.0), (2, 9.0), (2, 9.0)],
		wavecomp='potential',
		wavecomp_parameters=[(2, 3.0), (2, 3.0), (2, 3.0)],
		)
	}