awk '   BEGIN { 
		CA["N"]=-0.55;CA["D"]=-0.55;CA["S"]=-0.55;CA["H"]=-0.55;CA["F"]=-0.55;CA["W"]=-0.55;CA["Y"]=-0.55;CA["C"]=-0.55;
		CA["K"]=-0.69;CA["R"]=-0.69;CA["P"]=-0.69;
		CA["Q"]=-0.69;CA["E"]=-0.69;CA["M"]=-0.69;
		CA["A"]=-0.68;
		CA["I"]=-0.77;
		CA["L"]=-0.62;
		CA["T"]=-0.63;
		CA["V"]=-0.84;
		CB["N"]=-0.71;CB["D"]=-0.71;CB["S"]=-0.71;CB["H"]=-0.71;CB["F"]=-0.71;CB["W"]=-0.71;CB["Y"]=-0.71;CB["C"]=-0.71;
                CB["K"]=-1.11;CB["R"]=-1.11;CB["P"]=-1.11;
                CB["Q"]=-0.97;CB["E"]=-0.97;CB["M"]=-0.97;
                CB["A"]=-1.00;
                CB["I"]=-1.28;
                CB["L"]=-1.26;
                CB["T"]=-0.81;
                CB["V"]=-1.20;	
	}
	NF~4 && $1!="DATA"{
		if( $3=="CA") printf("%4d %1s %4s %8.3f\n", $1,$2,$3,$4-CA[$2]);
		else if( $3=="CB") printf("%4d %1s %4s %8.3f\n", $1,$2,$3,$4-CB[$2]);
		else printf("%4d %1s %4s %8.3f\n", $1,$2,$3,$4);
	}
	NF!=4 || $1=="DATA"{ print };
' $1

# add corrections from deuteration effects for CA and CB shifts
