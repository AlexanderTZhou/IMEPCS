# IMEPCS

AUTHOR: Alexander Zhou

INSTITUTION: Hong Kong University of Science and Technology

E-MAIL: atzhou@cse.ust.hk

# Arguments

argv[0] = ./epcli

argv[1] = Process Type

	*	1 : MEPCS
	*	2 : IMEPCS

argv[2] = Edge File

	*	tsv formal (ID1	ID2 Sign[-1 or 1])

argv[3] = Query Node
 
argv[4] = Epsilon Threshold Value

	*	Between 0.0 and 1.0 (Recommended >0.5)
	
argv[5] = Phi Threshold Value

	*	Between 0.0 and 1.0 (Recommended <0.5) 
	
argv[6] = Output File (optional)
	
	*	Prints the runtime and size into a file
