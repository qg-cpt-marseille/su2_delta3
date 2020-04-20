/* Copyright 2019 Giorgio Sarno, Pietro Donà and Francesco Gozzini */

/* sl2cfoam is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   sl2cfoam is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.*/

#include "delta3_4d.h"
#include "sl2cfoam.h"
#include "error.h"
#include "dbg.h"
#include "utilities.h"
#include "wigxjpf.h"
#include "jsymbols.h"
#include "common.h"



void j15_symbol_hash_delta3(struct data_folders* fd, 
							dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5, 
						    dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10, 
						    dspin two_i1, dspin two_i2, dspin two_i3) {

    char filename[1024];
    char path_6j[1024];

    // build path for 6j symbols
    sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.6j", two_j2, two_j3, two_j5, two_j4, two_j6, two_i1);
    strcpy(path_6j, fd->foursimp_hashtable6j_6j);
    strcat(path_6j, filename);

    //////////////////// Hash Table initialization ///////////////////

    // check for previously computed tables.
    khash_t(HashTableJ6) *h = NULL;

	if (file_exist(path_6j) != 0) {
		h = kh_load(HashTableJ6, path_6j);
	}

    if (h == NULL) {
        h = kh_init(HashTableJ6);
    }

    //////////////////// Initialize pointer for J6 symbols ////////////////////

    dspin* A = malloc(6 * sizeof(dspin));

    //////////////////// Start cycling all possible L combinations ////////////////////

    dspin two_limk4_1, two_limk4_2, 
		  two_limk5_1, two_limk5_2;

    two_limk4_1 = max(abs(two_j5-two_j6),abs(two_j10-two_j9));
    two_limk4_2 = min(two_j5+two_j6,two_j10+two_j9);

    two_limk5_1 = max(abs(two_j4-two_j6),abs(two_j7-two_j8));
    two_limk5_2 = min(two_j4+two_j6,two_j7+two_j8);

	///////////////////////////////////////
	// using IRREDUCIBLE coupling
	///////////////////////////////////////

	dspin two_limx1, two_limx2;

	for(dspin two_k4 = two_limk4_1; two_k4 <= two_limk4_2; two_k4 +=  2){
		for(dspin two_k5 = two_limk5_1; two_k5 <= two_limk5_2; two_k5 += 2) {
			
			two_limx1 = max(abs(two_i1-two_j7), abs(two_j5-two_k5));
			two_limx1 = max(two_limx1, abs(two_k4-two_j8));
			two_limx1 = max(two_limx1, abs(two_j10-two_i3));
			two_limx1 = max(two_limx1, abs(two_i2-two_j3));

			two_limx2 = min(two_i1+two_j7, two_j5+two_k5);
			two_limx2 = min(two_limx2, two_k4+two_j8);
			two_limx2 = min(two_limx2, two_j10+two_i3);
			two_limx2 = min(two_limx2, two_i2+two_j3);

			if (two_limx2 < two_limx1) {
				continue;
			}

			for (dspin two_x = two_limx1; two_x <= two_limx2; two_x += 2) {

				J6Symbol_Hash(h, &A,
							  two_i1, two_j7, two_x, 
							  two_k5, two_j5, two_j4 );
				J6Symbol_Hash(h, &A,
							  two_j5, two_k5, two_x, 
							  two_j8, two_k4, two_j6);                       
				J6Symbol_Hash(h, &A,
							  two_k4, two_j8, two_x, 
							  two_i3, two_j10, two_j9);             
				J6Symbol_Hash(h, &A,
							  two_j10, two_i3, two_x, 
							  two_j3, two_i2, two_j1);
				J6Symbol_Hash(h, &A,
							  two_i2, two_j3, two_x, 
							  two_i1, two_j7, two_j2);
			}
		}
	}

    free(A);

    // write hashtables to disk
    kh_write(HashTableJ6, h, path_6j);

    // free memory
    kh_destroy(HashTableJ6, h);

}

void j15_symbol_hash_delta3_2(kh_HashTableJ6_t *h,
							  dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5, 
						      dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10, 
						      dspin two_i1, dspin two_i2, dspin two_i3) {


    //////////////////// Initialize pointer for J6 symbols ////////////////////

    dspin* A = malloc(6 * sizeof(dspin));

    //////////////////// Start cycling all possible L combinations ////////////////////

    dspin two_limk4_1, two_limk4_2, 
		  two_limk5_1, two_limk5_2;

    two_limk4_1 = max(abs(two_j5-two_j6),abs(two_j10-two_j9));
    two_limk4_2 = min(two_j5+two_j6,two_j10+two_j9);

    two_limk5_1 = max(abs(two_j4-two_j6),abs(two_j7-two_j8));
    two_limk5_2 = min(two_j4+two_j6,two_j7+two_j8);

	///////////////////////////////////////
	// using IRREDUCIBLE coupling
	///////////////////////////////////////

	dspin two_limx1, two_limx2;

	for(dspin two_k4 = two_limk4_1; two_k4 <= two_limk4_2; two_k4 +=  2){
		for(dspin two_k5 = two_limk5_1; two_k5 <= two_limk5_2; two_k5 += 2) {
			
			two_limx1 = max(abs(two_i1-two_j7), abs(two_j5-two_k5));
			two_limx1 = max(two_limx1, abs(two_k4-two_j8));
			two_limx1 = max(two_limx1, abs(two_j10-two_i3));
			two_limx1 = max(two_limx1, abs(two_i2-two_j3));

			two_limx2 = min(two_i1+two_j7, two_j5+two_k5);
			two_limx2 = min(two_limx2, two_k4+two_j8);
			two_limx2 = min(two_limx2, two_j10+two_i3);
			two_limx2 = min(two_limx2, two_i2+two_j3);

			if (two_limx2 < two_limx1) {
				continue;
			}

			for (dspin two_x = two_limx1; two_x <= two_limx2; two_x += 2) {

				J6Symbol_Hash(h, &A,
							  two_i1, two_j7, two_x, 
							  two_k5, two_j5, two_j4 );
				J6Symbol_Hash(h, &A,
							  two_j5, two_k5, two_x, 
							  two_j8, two_k4, two_j6);                       
				J6Symbol_Hash(h, &A,
							  two_k4, two_j8, two_x, 
							  two_i3, two_j10, two_j9);             
				J6Symbol_Hash(h, &A,
							  two_j10, two_i3, two_x, 
							  two_j3, two_i2, two_j1);
				J6Symbol_Hash(h, &A,
							  two_i2, two_j3, two_x, 
							  two_i1, two_j7, two_j2);
			}
		}
	}

    free(A);

}

void sl2cfoam_hash_four_ampl_BF_delta3(struct data_folders* fd,
									   dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5, 
									   dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10, 
									   dspin two_i1, dspin two_i2, dspin two_i3) {                      

	// build path for 6js
	char filename[1024];
    char path_6j[1024];
    sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.6j", two_j2, two_j3, two_j5, two_j4, two_j6, two_i1);
    strcpy(path_6j, fd->foursimp_hashtable6j_6j);
    strcat(path_6j, filename);

	// load hash tables for 6Js
	khash_t(HashTableJ6) *h = kh_load(HashTableJ6, path_6j);

	char path_ampl[1024];

	// build path for amplitude values
	sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.bf", two_j2, two_j3, two_j5, two_j4, two_j6, two_i1);
	strcpy(path_ampl, fd->foursimp_bf_hashtableampl_ampl);
	strcat(path_ampl, filename);

	// initialize or load hash table of amplitudes on disk
	kh_HashTableEPRL_t *h3 = NULL;
	
	if (file_exist(path_ampl)!= 0) return;

	// table not found, initialize it
	h3 = kh_init(HashTableEPRL);
		
    //////////////////////////////////////////////////////////////////////////////
    // loop over i4 and i5
    //////////////////////////////////////////////////////////////////////////////

    dspin two_i4, two_i5;

	dspin two_limi4_1, two_limi4_2;
	dspin two_limi5_1, two_limi5_2;

    two_limi4_1 = max(abs(two_j5-two_j6),abs(two_j10-two_j9));
    two_limi4_2 = min(two_j5+two_j6,two_j10+two_j9);

    two_limi5_1 = max(abs(two_j4-two_j6),abs(two_j7-two_j8));
    two_limi5_2 = min(two_j4+two_j6,two_j7+two_j8);
        
    for (two_i4 = two_limi4_1; two_i4 <= two_limi4_2; two_i4 += 2) {
    for (two_i5 = two_limi5_1; two_i5 <= two_limi5_2; two_i5 += 2) {
		
        int ret; // return codes from khash

		// key for disk cache of current amplitude
		HashTableEPRL_key_t keyEPRL = {two_j1, two_j6, two_j7, two_j8, two_j9, two_j10,
									   two_i2, two_i3, two_i4, two_i5};

		// check if value has been already computed and stored
		if (kh_get(HashTableEPRL, h3, keyEPRL) != kh_end(h3)) {

			// value found!
			// go to next intertwiner set
			continue;

		}

		//////////////////////////////////////////////////////////////////////
		// value was not found or table was initialized
		// compute it then
		//////////////////////////////////////////////////////////////////////

		double ampl = 0.0;

		dspin *A = malloc(6 * sizeof(dspin));

		// Compute 15j Symbol
		ampl = J15Symbol(h, &A,
						two_j1, two_j2, two_j3, two_j4, two_j5,
						two_j6, two_j7, two_j8, two_j9, two_j10,
						two_i1, two_i2, two_i3, two_i4, two_i5);

		free(A);

		// normalize amplitude with dimensions of all 5 intertwiners
		ampl = ampl *  d(two_i1) * d(two_i2) * d(two_i3) * d(two_i4) * d(two_i5);
		
		////////////////////////////////////////////////////////
		// put key into amplitudes cache to be written to disk
		////////////////////////////////////////////////////////

		khint_t s = kh_put(HashTableEPRL, h3, keyEPRL, &ret);

		if (ret == -1) {
			error("error inserting key into EPRL hash table");
		}

		if (ret == 1) {

			kh_val(h3, s) = ampl;

			// check if insertion worked
			s = kh_get(HashTableEPRL, h3, keyEPRL);
			if (kh_val(h3, s) == 0 && ampl != 0) {

				// retry
				kh_val(h3, s) = ampl;

				// check again and now if not fail
				if (kh_val(h3, s) == 0 && ampl != 0) {
					error("error inserting value into EPRL hashtable");
				}

			}

		}

        } // i4
        } // i5 

        // write hash table for this x,i1 to disk
        kh_write(HashTableEPRL, h3, path_ampl);
        // free memory
        kh_destroy(HashTableJ6, h);
        kh_destroy(HashTableEPRL, h3);

}

// similar to the previous function but 
// one has to provide the hashtables to it
// usefull when itereting 
// over boundary intertwiners
void sl2cfoam_hash_four_ampl_BF_delta3_2(kh_HashTableJ6_t *h,
									 	 kh_HashTableEPRL_t *h3,
										 dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5, 
										 dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10, 
										 dspin two_i1, dspin two_i2, dspin two_i3) {                      

	
    //////////////////////////////////////////////////////////////////////////////
    // loop over i4 and i5
    //////////////////////////////////////////////////////////////////////////////

    dspin two_i4, two_i5;

	dspin two_limi4_1, two_limi4_2;
	dspin two_limi5_1, two_limi5_2;

    two_limi4_1 = max(abs(two_j5-two_j6),abs(two_j10-two_j9));
    two_limi4_2 = min(two_j5+two_j6,two_j10+two_j9);

    two_limi5_1 = max(abs(two_j4-two_j6),abs(two_j7-two_j8));
    two_limi5_2 = min(two_j4+two_j6,two_j7+two_j8);
        
    for (two_i4 = two_limi4_1; two_i4 <= two_limi4_2; two_i4 += 2) {
    for (two_i5 = two_limi5_1; two_i5 <= two_limi5_2; two_i5 += 2) {
		
        int ret; // return codes from khash

		// key for disk cache of current amplitude
		HashTableEPRL_key_t keyEPRL = {two_j1, two_j6, two_j7, two_j8, two_j9, two_j10,
									   two_i2, two_i3, two_i4, two_i5};

		// check if value has been already computed and stored
		if (kh_get(HashTableEPRL, h3, keyEPRL) != kh_end(h3)) {

			// value found!
			// go to next intertwiner set
			continue;

		}

		dspin *A = malloc(6 * sizeof(dspin));


		//////////////////////////////////////////////////////////////////////
		// value was not found or table was not initialized
		// compute it then
		//////////////////////////////////////////////////////////////////////

		double ampl = 0.0;

		// Compute 15j Symbol
		ampl = J15Symbol(h, &A,
						two_j1, two_j2, two_j3, two_j4, two_j5,
						two_j6, two_j7, two_j8, two_j9, two_j10,
						two_i1, two_i2, two_i3, two_i4, two_i5);


		// normalize amplitude with dimensions of all 5 intertwiners
		ampl = ampl *  d(two_i1) * d(two_i2) * d(two_i3) * d(two_i4) * d(two_i5);
		
		////////////////////////////////////////////////////////
		// put key into amplitudes cache to be written to disk
		////////////////////////////////////////////////////////

		khint_t s = kh_put(HashTableEPRL, h3, keyEPRL, &ret);

		if (ret == -1) {
			error("error inserting key into EPRL hash table");
		}

		if (ret == 1) {

			kh_val(h3, s) = ampl;

			// check if insertion worked
			s = kh_get(HashTableEPRL, h3, keyEPRL);
			if (kh_val(h3, s) == 0 && ampl != 0) {

				// retry
				kh_val(h3, s) = ampl;

				// check again and now if not fail
				if (kh_val(h3, s) == 0 && ampl != 0) {
					error("error inserting value into EPRL hashtable");
				}

			}

		}

		free(A);

        } // i4
        } // i5 

}

// TODO Check for the case tj2==tj3 and tj1 > tj3
dspin bound_x(dspin tj1, dspin tj2, dspin tj3) {

	if ( tj1 == tj2 && tj1 == tj3 && tj2 == tj3) return 0;
	if ( tj2 == tj3 && tj1 < tj3) return tj1 - tj2 - tj3;
	if ( tj2 > tj3 && tj1 >= tj2) return tj2 - tj1 - tj3;
	if ( tj2 > tj3 && tj1 <  tj2) return tj1 - tj2 - tj3;
	if ( tj2 < tj3 && tj1 >= tj3) return tj3 - tj1 - tj2;
	if ( tj2 < tj3 && tj1 < tj3 ) return tj1 - tj2 - tj3;

	error("error in internal spin bound");
	
}

/**********************************************************************/

// to assembly the delta3 amplitude in the 4d BF case
// NB one has to precompute all symbols
// before using it
double delta3_bf4_creator (struct data_folders* fd, 
						   dspin tj_124, dspin tj_245, dspin tj_234, dspin tj_235, dspin tj_345, dspin tx,
						   dspin tj_125, dspin tj_123, dspin tj_134, dspin tj_145, dspin tj_256, dspin tj_236,
						   dspin tj_126, dspin tj_136, dspin tj_356, dspin tj_156, dspin tj_346, dspin tj_146, dspin tj_456,
						   dspin ti_1, dspin ti_2, dspin ti_3, dspin ti_4, dspin ti_5,
						   dspin ti_6, dspin ti_7, dspin ti_8, dspin ti_9) {

	// here we load hashed 15j symbols
	char local_filename[128];
	char path_ampl_1[1024], path_ampl_2[1024], path_ampl_3[1024];

	// build path for amplitude values
	sprintf(local_filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_245, tj_234, tj_345, tj_235, tx, ti_1);
	strcpy(path_ampl_1, fd->foursimp_bf_hashtableampl_ampl);
	strcat(path_ampl_1, local_filename);

	sprintf(local_filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_236, tj_126, tj_123, tj_136, tx, ti_4);
	strcpy(path_ampl_2, fd->foursimp_bf_hashtableampl_ampl);
	strcat(path_ampl_2, local_filename);

	sprintf(local_filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_146, tj_456, tj_156, tj_145, tx, ti_7);
	strcpy(path_ampl_3, fd->foursimp_bf_hashtableampl_ampl);
	strcat(path_ampl_3, local_filename);

	khash_t(HashTableBF4_1)* tablebf_1 = NULL;
	khash_t(HashTableBF4_2)* tablebf_2 = NULL;
	khash_t(HashTableBF4_3)* tablebf_3 = NULL;

	// finally load tables
	if (file_exist(path_ampl_1) && file_exist(path_ampl_2) && file_exist(path_ampl_3)) {

		tablebf_1 = kh_load(HashTableBF4_1, path_ampl_1);
		tablebf_2 = kh_load(HashTableBF4_2, path_ampl_2);
		tablebf_3 = kh_load(HashTableBF4_3, path_ampl_3);

	} else {
		error("files %s - %s - %s not found, table for given spins not precomputed", path_ampl_1, path_ampl_2, path_ampl_3);
	}		  		

	// bounds for internal intertwiners
	dspin tk3_min, tk3_max;

	tk3_min = (dspin) max(abs(tj_356-tj_136), abs(tx-tj_156));
	tk3_max = (dspin) min(tj_356 + tj_136, tx + tj_156); 
	
	dspin tk1_min, tk1_max;

	tk1_min = (dspin) max(abs(tj_134-tj_145), abs(tx-tj_345));
	tk1_max = (dspin) min(tj_134 + tj_145, tx + tj_345);

	dspin tk2_min, tk2_max;

	tk2_min = (dspin) max(abs(tj_125-tj_235), abs(tx-tj_123));
	tk2_max = (dspin) min(tj_125 + tj_235, tx + tj_123);

	double delta3_x = 0.0;
	double err = 0.0;

	for (dspin tk_1 = tk1_min; tk_1 <= tk1_max; tk_1 += 2) {
	for (dspin tk_2 = tk2_min; tk_2 <= tk2_max; tk_2 += 2) {
	for (dspin tk_3 = tk3_min; tk_3 <= tk3_max; tk_3 += 2) { 

		// now get values from tables

		double simplex_1, simplex_2, simplex_3;

		// get first 4 simplex from table
		HashTableBF4_1_key_t keybf_1 = {tj_124, tx, tj_125, tj_123, tj_134, tj_145,
										ti_2, ti_3, tk_1, tk_2};
		khint_t s_1 = kh_get(HashTableBF4_1, tablebf_1, keybf_1);

		if (s_1 != kh_end(tablebf_1)) {		
			simplex_1 = kh_val(tablebf_1, s_1);
		} else {
			error("value for delta3 not found");
		}

		// get second 4 simplex from table
		HashTableBF4_2_key_t keybf_2 = {tj_256, tx, tj_356, tj_156, tj_125, tj_235,
										ti_5, ti_6, tk_2, tk_3};
		khint_t s_2 = kh_get(HashTableBF4_2, tablebf_2, keybf_2);

		if (s_2 != kh_end(tablebf_2)) {		
			simplex_2 = kh_val(tablebf_2, s_2);
		} else {
			error("value for delta3 not found");
		}

		// get third 4 simplex from table
		HashTableBF4_3_key_t keybf_3 = {tj_346, tx, tj_134, tj_345, tj_356, tj_136,
										ti_8, ti_9, tk_3, tk_1};
		khint_t s_3 = kh_get(HashTableBF4_3, tablebf_3, keybf_3);

		if (s_3 != kh_end(tablebf_3)) {		
			simplex_3 = kh_val(tablebf_3, s_3);
		} else {
			error("value for delta3 not found");
		}

		mpfr_t simplexMPFR_1, simplexMPFR_2, simplexMPFR_3;

		mpfr_init_set_d(simplexMPFR_1, simplex_1, MPFR_RNDN);
		mpfr_init_set_d(simplexMPFR_2, simplex_2, MPFR_RNDN);
		mpfr_init_set_d(simplexMPFR_3, simplex_3, MPFR_RNDN);

		// multiply simplices as MPFR variables
		mpfr_mul(simplexMPFR_1, simplexMPFR_1, simplexMPFR_2, MPFR_RNDN);
		mpfr_mul(simplexMPFR_1, simplexMPFR_1, simplexMPFR_3, MPFR_RNDN);

		// multiply for the virtual spin dimension and phase
		mpfr_mul_si(simplexMPFR_1, simplexMPFR_1, d(tx)*real_negpow(3*tx + 3*(tk_1 + tk_2 + tk_3)), MPFR_RNDN);
		mpfr_div_si(simplexMPFR_1, simplexMPFR_1, d(tk_1), MPFR_RNDN);
		mpfr_div_si(simplexMPFR_1, simplexMPFR_1, d(tk_2), MPFR_RNDN);
		mpfr_div_si(simplexMPFR_1, simplexMPFR_1, d(tk_3), MPFR_RNDN);

		compsum_mpfr(&err, &delta3_x, simplexMPFR_1);  

		mpfr_clears(simplexMPFR_1, simplexMPFR_2, simplexMPFR_3, NULL);   

	} //k3
	} //k2
	} //k1

	kh_destroy(HashTableBF4_1, tablebf_1);
	kh_destroy(HashTableBF4_2, tablebf_2);
	kh_destroy(HashTableBF4_3, tablebf_3);

	return delta3_x;
}

/**********************************************************************/

// function to hash a single 4 simplex
// with three coherent states
// we resum over three coherent states
double complex delta3_simplex_coherent_bf4_hashing (struct data_folders* fd,
													kh_HashTableCoherent_t *hcoherent,
													kh_HashTableCoherentBF_t *hbfcoherent,
													dspin tj_124, dspin tj_245, dspin tj_234, dspin tj_235, dspin tj_345,
													dspin tx, dspin tj_125, dspin tj_123, dspin tj_134, dspin tj_145, 
													dspin tr_i1[2], dspin tr_i2[2], dspin tr_i3[2],
													dspin tk_1, dspin tk_2) {
	// threshold to avoid zeros
	double threshold = 1e-40;
	// bounds for internal intertwiners

	double complex coherent_simplex = 0.0 + 0.0*I;
	double complex err  = 0.0 + 0.0*I;
	
	for (dspin ti_1 = tr_i1[0]; ti_1 <= tr_i1[1]; ti_1 +=2) {

		char local_filename[128];
		char path_ampl_1[1024];
		
		// build path for amplitude values
		sprintf(local_filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_245, tj_234, tj_345, tj_235, tx, ti_1);
		strcpy(path_ampl_1, fd->foursimp_bf_hashtableampl_ampl);
		strcat(path_ampl_1, local_filename);

		// load the table 
		khash_t(HashTableBF4_1)* tablebf_1 = NULL;
		
		// check if the file exists and load
		if (file_exist(path_ampl_1) != 0) {
			tablebf_1 = kh_load(HashTableBF4_1, path_ampl_1);
		} else {
			error("files %s not found, table for given spins not precomputed", path_ampl_1);
		}

	for (dspin ti_2 = tr_i2[0]; ti_2 <= tr_i2[1]; ti_2 +=2) { 
	for (dspin ti_3 = tr_i3[0]; ti_3 <= tr_i3[1]; ti_3 +=2) { 
		
		// get coherent state 
		khint_t c[3];
		HashTableCoherent_key_t coh1 = {0, ti_1};
		c[0] = kh_get(HashTableCoherent, hcoherent, coh1);

		HashTableCoherent_key_t coh2 = {1, ti_2};
		c[1] = kh_get(HashTableCoherent, hcoherent, coh2);

		HashTableCoherent_key_t coh3 = {2, ti_3};
		c[2] = kh_get(HashTableCoherent, hcoherent, coh3);

		double complex coherent = 0.0 + 0.0*I;

		coherent = mpfr_complex_mul(kh_value(hcoherent, c[0]), kh_value(hcoherent, c[1]));
		coherent = mpfr_complex_mul(coherent, kh_value(hcoherent, c[2]));

		double simplex_1 = 0.0 + 0.0*I;

		// get first 4 simplex from table
		HashTableBF4_1_key_t keybf_1 = {tj_124, tx, tj_125, tj_123, tj_134, tj_145,
										ti_2, ti_3, tk_1, tk_2};
		khint_t s_1 = kh_get(HashTableBF4_1, tablebf_1, keybf_1);

		if (s_1 != kh_end(tablebf_1)) {		
			simplex_1 = kh_val(tablebf_1, s_1);
		} else {
			error("value for delta3 not found");
		}

		mpfr_t  valMPFR;
		mpfr_init_set_d (valMPFR, simplex_1, MPFR_RNDN);

		mpfr_t  coherentMPFR_Re, coherentMPFR_Im;

		mpfr_init_set_d (coherentMPFR_Re, creal(coherent), MPFR_RNDN);
		mpfr_init_set_d (coherentMPFR_Im, cimag(coherent), MPFR_RNDN);

		mpfr_mul (coherentMPFR_Re, coherentMPFR_Re, valMPFR, MPFR_RNDN);
		mpfr_mul (coherentMPFR_Im, coherentMPFR_Im, valMPFR, MPFR_RNDN);

		compsum_mpfr_complex (&err, &coherent_simplex, coherentMPFR_Re, coherentMPFR_Im );

		mpfr_clears ( valMPFR, coherentMPFR_Re, coherentMPFR_Im, NULL);

	} // i3
	} // i2

		kh_destroy(HashTableBF4_1, tablebf_1);

	} // i1

	return coherent_simplex;
 }



// to assembly the delta3 amplitute in the 4d BF case
// and compute the angle
// NB one has to precompute all symbols 
double delta3_angles_bf4_creator (struct data_folders* fd, 
								  dspin tj_124, dspin tj_245, dspin tj_234, dspin tj_235, dspin tj_345, dspin tx,
								  dspin tj_125, dspin tj_123, dspin tj_134, dspin tj_145, dspin tj_256, dspin tj_236,
								  dspin tj_126, dspin tj_136, dspin tj_356, dspin tj_156, dspin tj_346, dspin tj_146, dspin tj_456,
								  dspin ti1_min, dspin ti1_max, dspin ti2_min, dspin ti2_max, dspin ti3_min, dspin ti3_max,
								  dspin ti4_min, dspin ti4_max, dspin ti5_min, dspin ti5_max, dspin ti6_min, dspin ti6_max,
								  dspin ti7_min, dspin ti7_max, dspin ti8_min, dspin ti8_max, dspin ti9_min, dspin ti9_max) {

	double delta3 = 0.0;
	double errd3  = 0.0;
	// threshold to avoid zeros
	double threshold = 1e-40;

	// bounds for internal intertwiners

	dspin tk1_min, tk1_max;

	tk1_min = (dspin) max(abs(tj_134-tj_145), abs(tx-tj_345));
	tk1_max = (dspin) min(tj_134 + tj_145, tx + tj_345);

	dspin tk2_min, tk2_max;

	tk2_min = (dspin) max(abs(tj_125-tj_235), abs(tx-tj_123));
	tk2_max = (dspin) min(tj_125 + tj_235, tx + tj_123);

	dspin tk3_min, tk3_max;

	tk3_min = (dspin) max(abs(tj_356-tj_136), abs(tx-tj_156));
	tk3_max = (dspin) min(tj_356 + tj_136, tx + tj_156); 

	// start cycling intertwiners	
	for (dspin ti_1 = ti1_min; ti_1 <= ti1_max; ti_1 +=2) { 

		
		char local_filename[128];
		char path_ampl_1[1024];
		
		// build path for amplitude values
		sprintf(local_filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_245, tj_234, tj_345, tj_235, tx, ti_1);
		strcpy(path_ampl_1, fd->foursimp_bf_hashtableampl_ampl);
		strcat(path_ampl_1, local_filename);

		// load the table 
		khash_t(HashTableBF4_1)* tablebf_1 = NULL;
		
		// check if the file exists and load
		if (file_exist(path_ampl_1)) {
			tablebf_1 = kh_load(HashTableBF4_1, path_ampl_1);
		} else {
			error("files %s not found, table for given spins not precomputed", path_ampl_1);
		}		 


	for (dspin ti_4 = ti4_min; ti_4 <= ti4_max; ti_4 +=2) { 

		
		char path_ampl_2[1024];

		sprintf(local_filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_236, tj_126, tj_123, tj_136, tx, ti_4);
		strcpy(path_ampl_2, fd->foursimp_bf_hashtableampl_ampl);
		strcat(path_ampl_2, local_filename);

		khash_t(HashTableBF4_2)* tablebf_2 = NULL;
		
		if (file_exist(path_ampl_2)) {
			tablebf_2 = kh_load(HashTableBF4_2, path_ampl_2);
		} else {
			error("files %s not found, table for given spins not precomputed", path_ampl_2);
		}	

	for (dspin ti_7 = ti7_min; ti_7 <= ti7_max; ti_7 +=2) { 


		char path_ampl_3[1024];

		sprintf(local_filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_146, tj_456, tj_156, tj_145, tx, ti_7);
		strcpy(path_ampl_3, fd->foursimp_bf_hashtableampl_ampl);
		strcat(path_ampl_3, local_filename);

		khash_t(HashTableBF4_3)* tablebf_3 = NULL;
		
		if (file_exist(path_ampl_3)) {
			tablebf_3 = kh_load(HashTableBF4_3, path_ampl_3);
		} else {
			error("files %s not found, table for given spins not precomputed", path_ampl_3);
		}		 

	for (dspin ti_2 = ti2_min; ti_2 <= ti2_max; ti_2 +=2) { 
	for (dspin ti_3 = ti3_min; ti_3 <= ti3_max; ti_3 +=2) { 
	for (dspin ti_5 = ti5_min; ti_5 <= ti2_max; ti_5 +=2) { 
	for (dspin ti_6 = ti6_min; ti_6 <= ti3_max; ti_6 +=2) { 
	for (dspin ti_8 = ti8_min; ti_8 <= ti8_max; ti_8 +=2) { 
	for (dspin ti_9 = ti9_min; ti_9 <= ti9_max; ti_9 +=2) { 

		double delta3_x = 0.0;
		double err = 0.0;
		
		for (dspin tk_1 = tk1_min; tk_1 <= tk1_max; tk_1 += 2) {
		for (dspin tk_2 = tk2_min; tk_2 <= tk2_max; tk_2 += 2) {
		for (dspin tk_3 = tk3_min; tk_3 <= tk3_max; tk_3 += 2) { 

			// now get values from tables

			double simplex_1, simplex_2, simplex_3;

			// get first 4 simplex from table
			HashTableBF4_1_key_t keybf_1 = {tj_124, tx, tj_125, tj_123, tj_134, tj_145,
											ti_2, ti_3, tk_1, tk_2};
			khint_t s_1 = kh_get(HashTableBF4_1, tablebf_1, keybf_1);

			if (s_1 != kh_end(tablebf_1)) {		
				simplex_1 = kh_val(tablebf_1, s_1);
			} else {
				error("value for delta3 not found");
			}

			if(fabs(simplex_1) < threshold) continue;

			// get second 4 simplex from table
			HashTableBF4_2_key_t keybf_2 = {tj_256, tx, tj_356, tj_156, tj_125, tj_235,
											ti_5, ti_6, tk_2, tk_3};
			khint_t s_2 = kh_get(HashTableBF4_2, tablebf_2, keybf_2);

			if (s_2 != kh_end(tablebf_2)) {		
				simplex_2 = kh_val(tablebf_2, s_2);
			} else {
				error("value for delta3 not found");
			}

			if(fabs(simplex_2) < threshold) continue;

			// get third 4 simplex from table
			HashTableBF4_3_key_t keybf_3 = {tj_346, tx, tj_134, tj_345, tj_356, tj_136,
											ti_8, ti_9, tk_3, tk_1};
			khint_t s_3 = kh_get(HashTableBF4_3, tablebf_3, keybf_3);

			if (s_3 != kh_end(tablebf_3)) {		
				simplex_3 = kh_val(tablebf_3, s_3);
			} else {
				error("value for delta3 not found");
			}
				
			if(fabs(simplex_3) < threshold) continue;

			mpfr_t simplexMPFR_1, simplexMPFR_2, simplexMPFR_3;

			mpfr_init_set_d(simplexMPFR_1, simplex_1, MPFR_RNDN);
			mpfr_init_set_d(simplexMPFR_2, simplex_2, MPFR_RNDN);
			mpfr_init_set_d(simplexMPFR_3, simplex_3, MPFR_RNDN);

			// multiply simplices as MPFR variables
			mpfr_mul(simplexMPFR_1, simplexMPFR_1, simplexMPFR_2, MPFR_RNDN);
			mpfr_mul(simplexMPFR_1, simplexMPFR_1, simplexMPFR_3, MPFR_RNDN);

			// multiply for the virtual spin dimension and phase
			mpfr_mul_si(simplexMPFR_1, simplexMPFR_1, d(tx)*real_negpow(3*tx + 3*(tk_1 + tk_2 + tk_3)), MPFR_RNDN);
			mpfr_div_si(simplexMPFR_1, simplexMPFR_1, d(tk_1), MPFR_RNDN);
			mpfr_div_si(simplexMPFR_1, simplexMPFR_1, d(tk_2), MPFR_RNDN);
			mpfr_div_si(simplexMPFR_1, simplexMPFR_1, d(tk_3), MPFR_RNDN);

			compsum_mpfr(&err, &delta3_x, simplexMPFR_1);  

			if (fabs(delta3_x) > 100) {
				//printf("%i %i %i %i %17g %17g %17g %17g \n", tx, tk_1,tk_2,tk_3, 
				//		simplex_1,simplex_2,simplex_3,delta3_x);
				error("non consistent value for delta3_x");
			}

			mpfr_clears(simplexMPFR_1, simplexMPFR_2, simplexMPFR_3, NULL);   

		} //k3
		} //k2
		} //k1

		// TODO: METTERE QUI OPERAZIONI CON DELTA3_X

	} // i9
	} // i8
	} // i6
	} // i5
	} // i3
	} // i2
		kh_destroy(HashTableBF4_3, tablebf_3);
	} // i7
		kh_destroy(HashTableBF4_2, tablebf_2);
	} // i4		
		kh_destroy(HashTableBF4_1, tablebf_1);
	} // i1

	return delta3;

}




/**********************************************************************/

// EPRL

// to assembly the delta3 amplitude in the 4d EPRL case
// NB one has to precompute all symbols boosters and 
// 4 simplices before using it
double delta3_eprl4_creator (struct data_folders* fd, 
                            dspin tj_124, dspin tj_245, dspin tj_234, dspin tj_235, dspin tj_345, dspin tx,
                            dspin tj_125, dspin tj_123, dspin tj_134, dspin tj_145, dspin tj_256, dspin tj_236,
                            dspin tj_126, dspin tj_136, dspin tj_356, dspin tj_156, dspin tj_346, dspin tj_146, dspin tj_456,
                            dspin ti_1, dspin ti_2, dspin ti_3, dspin ti_4, dspin ti_5,
                            dspin ti_6, dspin ti_7, dspin ti_8, dspin ti_9,
							dspin td, float immirzi) {

	// bounds for internal intertwiners
	dspin tk3_min, tk3_max;

	tk3_min = (dspin) max(abs(tj_356-tj_136), abs(tx-tj_156));
	tk3_max = (dspin) min(tj_356 + tj_136, tx + tj_156); 
	
	dspin tk1_min, tk1_max;

	tk1_min = (dspin) max(abs(tj_134-tj_145), abs(tx-tj_345));
	tk1_max = (dspin) min(tj_134 + tj_145, tx + tj_345);

	dspin tk2_min, tk2_max;

	tk2_min = (dspin) max(abs(tj_125-tj_235), abs(tx-tj_123));
	tk2_max = (dspin) min(tj_125 + tj_235, tx + tj_123);

	double delta3_x = 0.0;
	double err = 0.0;
	
	// threshold to avoid zeros
	double threshold = 1e-40;
	
	// initialize tables 
	khash_t(HashTableEPRL_1)* tableeprl_1 = NULL;		
	khash_t(HashTableEPRL_2)* tableeprl_2 = NULL;
	khash_t(HashTableEPRL_3)* tableeprl_3 = NULL;

	for (dspin tk_1 = tk1_min; tk_1 <= tk1_max; tk_1 += 2) {

		// here we load hashed 15j symbols
		char local_filename[128];
		char path_ampl_1[1024], path_ampl_2[1024], path_ampl_3[1024];

		// build path for amplitude values, first 4 simplex
		sprintf(local_filename, "%i.%i.%i.%i_%i_%i.eprl", tj_345, tx, tj_134, tj_145, tk_1, td);
		strcpy(path_ampl_1, fd->foursimp_imm_hashtableampl_ampl);
		strcat(path_ampl_1, local_filename);
		
		// finally load tables
		if (file_exist(path_ampl_1)) {
			tableeprl_1 = kh_load(HashTableEPRL_1, path_ampl_1);
		} else {
			error("files %s not found, table for given spins not precomputed", path_ampl_1);
		}	

	for (dspin tk_2 = tk2_min; tk_2 <= tk2_max; tk_2 += 2) {
		
		// load second 4 simplex
		sprintf(local_filename, "%i.%i.%i.%i_%i_%i.eprl", tj_123, tx, tj_125, tj_235, tk_2, td);
		strcpy(path_ampl_2, fd->foursimp_imm_hashtableampl_ampl);
		strcat(path_ampl_2, local_filename);

		if (file_exist(path_ampl_2)) {
			tableeprl_2 = kh_load(HashTableEPRL_2, path_ampl_2);
		} else {
			error("files %s not found, table for given spins not precomputed", path_ampl_2);
		}	

		double simplex_1, simplex_2, simplex_3;

		// we already have the value for the first simplex
		// get first 4 simplex from table
		HashTableEPRL_1_key_t keyeprl_1 = {tj_124, tj_245, tj_234, tj_235, tj_125, tj_123,
											ti_1, ti_2, ti_3, tk_2};
		
		khint_t s_1 = kh_get(HashTableEPRL_1, tableeprl_1, keyeprl_1);
		if (s_1 != kh_end(tableeprl_1)) {		
			simplex_1 = kh_val(tableeprl_1, s_1);
		} else {
			error("value for delta3 not found");
		}

		if (fabs(simplex_1) < threshold) continue;
		
	for (dspin tk_3 = tk3_min; tk_3 <= tk3_max; tk_3 += 2) { 
		
		// load third 4 simplex
		sprintf(local_filename, "%i.%i.%i.%i_%i_%i.eprl", tj_156, tx, tj_356, tj_136, tk_3, td);
		strcpy(path_ampl_3, fd->foursimp_imm_hashtableampl_ampl);
		strcat(path_ampl_3, local_filename);

		if (file_exist(path_ampl_3)) {
			tableeprl_3 = kh_load(HashTableEPRL_3, path_ampl_3);
		} else {
			error("files %s not found, table for given spins not precomputed", path_ampl_3);
		}	

		// get second 4 simplex from table
		HashTableEPRL_2_key_t keyeprl_2 = {tj_256, tj_236, tj_126, tj_136, tj_356, tj_156,
											ti_4, ti_5, ti_6, tk_3};
		
		khint_t s_2 = kh_get(HashTableEPRL_2, tableeprl_2, keyeprl_2);
		if (s_2 != kh_end(tableeprl_2)) {		
			simplex_2 = kh_val(tableeprl_2, s_2);
		} else {
			error("value for delta3 not found");
		}

		if (fabs(simplex_2) < threshold) continue;

		// get third 4 simplex from table
		HashTableEPRL_3_key_t keyeprl_3 = {tj_346, tj_146, tj_456, tj_145, tj_134, tj_345,
											ti_7, ti_8, ti_9, tk_1};
		
		khint_t s_3 = kh_get(HashTableEPRL_3, tableeprl_3, keyeprl_3);
		if (s_3 != kh_end(tableeprl_3)) {		
			simplex_3 = kh_val(tableeprl_3, s_3);
		} else {
			error("value for delta3 not found");
		}

		if (fabs(simplex_3) < threshold) continue;

		// assembly everything at arbitraty precision

		mpfr_t simplexMPFR_1, simplexMPFR_2, simplexMPFR_3;

		mpfr_init_set_d(simplexMPFR_1, simplex_1, MPFR_RNDN);
		mpfr_init_set_d(simplexMPFR_2, simplex_2, MPFR_RNDN);
		mpfr_init_set_d(simplexMPFR_3, simplex_3, MPFR_RNDN);

		// multiply simplices as MPFR variables
		mpfr_mul(simplexMPFR_1, simplexMPFR_1, simplexMPFR_2, MPFR_RNDN);
		mpfr_mul(simplexMPFR_1, simplexMPFR_1, simplexMPFR_3, MPFR_RNDN);

		// multiply for the virtual spin dimension and phase
		mpfr_mul_si(simplexMPFR_1, simplexMPFR_1, d(tx)*real_negpow(3*tx + 3*(tk_1 + tk_2 + tk_3)), MPFR_RNDN);

		// divide for additional internal int dimensions
		mpfr_div_si(simplexMPFR_1, simplexMPFR_1, d(tk_1), MPFR_RNDN);
		mpfr_div_si(simplexMPFR_1, simplexMPFR_1, d(tk_2), MPFR_RNDN);
		mpfr_div_si(simplexMPFR_1, simplexMPFR_1, d(tk_3), MPFR_RNDN);

		compsum_mpfr(&err, &delta3_x, simplexMPFR_1);  

		mpfr_clears(simplexMPFR_1, simplexMPFR_2, simplexMPFR_3, NULL);   
		
		kh_destroy(HashTableEPRL_3, tableeprl_3);

	} // k3
		kh_destroy(HashTableEPRL_2, tableeprl_2);
	} // k2
		kh_destroy(HashTableEPRL_1, tableeprl_1);
	} // k1

	return delta3_x;

}