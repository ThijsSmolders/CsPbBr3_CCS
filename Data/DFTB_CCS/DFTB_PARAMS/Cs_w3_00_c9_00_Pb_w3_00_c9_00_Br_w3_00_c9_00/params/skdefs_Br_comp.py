atomconfigs = {
    "Br": AtomConfig(
       znuc=35,
       mass=79.904,
       occupations=[
           [ [ 1.0, 1.0 ], [ 1.0, 1.0 ], [ 1.0, 1.0 ], [ 1.0, 1.0 ] ],  # s
           [               [ 3.0, 3.0 ], [ 3.0, 3.0 ], [ 3.0, 2.0 ] ],  # p
           [                             [ 5.0, 5.0 ], [ 0.0, 0.0 ] ],  # d suppress_in_twocnt
           ],
       valenceqns=[
           [4, ],  # s
           [4, ],  # p
           [ ],  # d suppress_in_twocnt
           ],
       relativistic=True,
       ),

	}

skbases = {
   "Br": SlaterBasis(
        exponents=[[0.896, 2.24, 5.6, 14., 35.], # s
                   [0.896, 2.24, 5.6, 14., 35.], # p
                   [0.896, 2.24, 5.6, 14., 35.], # d                   
                   ],
        maxpowers=[ 3,  # s
                    3,  # p
                    3,  # d 
                    ]
        ),
	}
compressions = {
	 "Br": Compression(
		potcomp='potential',
		potcomp_parameters=[(2, 9.0), (2, 9.0), (2, 9.0)],
		wavecomp='potential',
		wavecomp_parameters=[(2, 3.0), (2, 3.0), (2, 3.0)],
		)
	}