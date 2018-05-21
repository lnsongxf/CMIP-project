**************************************************************
README.txt: se presenta la notación usada para las figuras
**************************************************************

**************************
[I] Notación general:
**************************
Todas las figuras de las simulaciones han seguido la siguiente estructura:
           
       'fig' + 'simul_name' + 'number' .pdf

**************************
[II] 'simul_name':
**************************
Los ejercicios presentados en el slides son los siguientes (en el mismo orden): 
		
	fp, rpa2, bpa2, bp_zlb, bp_zlb4, y bp_zlb3. 

Las simulaciones que se han realizado hasta ahora son:

	[1] Fisherian Doc:

	'fp'     : Fisher's doc and credit crunch.

	[2] Fiscal Policy:

	'rpa1'   : Unanticipated fiscal policy.
	'rpa2'   : Anticipated fiscal policy. 

	[3] Credit channel:

	'bpa1'   : Unanticipated decrease in real spread. 
	'bpa2'   : Anticipated decrease in real spread. 
	'bpb1'   : Unanticipated increase in real spread. 
	'bpb2'   : Anticipated increase in real spread. 


	[4] Credit Crunch, MP reaction, and ZLB:

	'bp_zlb1' : Unanticipated credit crunch and ZLB with MP reaction.
	'bp_zlb2' : Anticipated credit crunch and ZLB with MP reaction. 
	'bp_zlb3' : Unanticipated negative interest rate on reserves.
	'bp_zlb4' : Anticipated negative interest rate on reserves. 


**************************
[III] 'number':
**************************
En todos los ejercicios se hace usa la siguiente notación

	For simulations with zero real spread (e.g., 'fp' and 'rpa1'):
	  	1: Nominal lending and deposit rates at the stationary equilibrium.
		2: Price index path. (constant aggregate money MP).
		3: Fisher decomposition (constant aggregate money MP).
		4: Price index path. (Taylor rule).
		5: Fisher decomposition (Taylor rule).
	
	For simulations with positive real spread target (e.g., 'bpa1' and 'bp_zlb2'):
	  	1: Nominal lending and deposit rates at the stationary equilibrium.
		2: Real spreads at the stationary equilibrium.
		3: Price index path.
		4: Currency and reserves.
		5: Fisher decomposition.
	
	6: Certainty equivalent.
	7: Lending and deposit rates.
	8: Government's equity.
	9: Real transfers.
       10: Real activity.
       11: Credit and Deposit response. 
       12: Steady-state distribution vs ex-ante distribution
       13: Mass of restricted agents.
       14: 3-D wealth's distribution.	 	
       15: 2-D wealth's distribution (only surface, 'heatmap').
       16: 3-D wealth's distribution (mesh).

En los slides se presenta los resultados siguiendo el siguiente orden:	
	[1] Real variables' evolution:	15,11,7,10	
		
	[2] Nominal variables' evolution (*): 5,4,3

	(*) En el ejercicio de Fisher doc ('fp'): 3,2,5,4

*******************************
[IV] Observación:   
*******************************
Se debe notar los indices de cada ejercicio, por ejemplo:

	'fig' + 'rpa1' + '15' == 'figrpa115.pdf'

