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
