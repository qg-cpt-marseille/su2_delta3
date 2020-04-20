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

#ifndef __SL2CFOAM_DELTA3_BF4_H__
#define __SL2CFOAM_DELTA3_BF4_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "khash.h"
#include "config.h"
#include "hashing.h"

#define MAX_J_WIG 1000

// auxiliary tables for the 3 simplices amplitudes
ARRAY_TABLE_INIT(HashTableBF4_1, 10, double);
ARRAY_TABLE_INIT(HashTableBF4_2, 10, double);
ARRAY_TABLE_INIT(HashTableBF4_3, 10, double);

// 2 spins key for coherent states
ARRAY_TABLE_INIT(HashTableCoherent, 2, double complex);

// 3 key for coherent 4 simplex
ARRAY_TABLE_INIT(HashTableCoherentBF, 4 , double complex);


////////////////////////////////////////////////////////////////
// Functions to compute the 4 simplices delta 3 amplitude BF
////////////////////////////////////////////////////////////////

/**********************************************************************/

// hash all 6j symbols to compose the 
// 4d delta3, namely to construct 15j
void j15_symbol_hash_delta3(struct data_folders* fd, 
							dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5, 
						    dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10, 
						    dspin two_i1, dspin two_i2, dspin two_i3); 

// similar to the previous function but 
// one has to provide the hashtables to it
// usefull when itereting 
// over boundary intertwiners	
void j15_symbol_hash_delta3_2(kh_HashTableJ6_t *h,
							  dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5, 
						      dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10, 
						      dspin two_i1, dspin two_i2, dspin two_i3);							

// hash all 15j symbols involved
// in the delta3 computation
void sl2cfoam_hash_four_ampl_BF_delta3(struct data_folders* fd,
									   dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5, 
									   dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10, 
									   dspin two_i1, dspin two_i2, dspin two_i3);

// similar to the previous function but 
// one has to provide the hashtables to it
// usefull when itereting 
// over boundary intertwiners									   
void sl2cfoam_hash_four_ampl_BF_delta3_2(kh_HashTableJ6_t *h,
									 	 kh_HashTableEPRL_t *h3,
										dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4, dspin two_j5, 
										dspin two_j6, dspin two_j7, dspin two_j8, dspin two_j9, dspin two_j10, 
										dspin two_i1, dspin two_i2, dspin two_i3);									                        

// utility to find the douns
// over the internal spin
dspin bound_x(dspin tj1, dspin tj2, dspin tj3);


// to assembly the delta3 amplitu in the 4d BF case
// NB one has to precompute all symbols
// before using it
double delta3_bf4_creator (struct data_folders* fd, 
						   dspin tj_124, dspin tj_245, dspin tj_234, dspin tj_235, dspin tj_345, dspin tx,
						   dspin tj_125, dspin tj_123, dspin tj_134, dspin tj_145, dspin tj_256, dspin tj_236,
						   dspin tj_126, dspin tj_136, dspin tj_356, dspin tj_156, dspin tj_346, dspin tj_146, dspin tj_456,
						   dspin ti_1, dspin ti_2, dspin ti_3, dspin ti_4, dspin ti_5,
						   dspin ti_6, dspin ti_7, dspin ti_8, dspin ti_9);

// function to hash a single 4 simplex
// with three coherent states
// we resum over three coherent states
double complex delta3_simplex_coherent_bf4_hashing (struct data_folders* fd,
													kh_HashTableCoherent_t *hcoherent,
													kh_HashTableCoherentBF_t *hbfcoherent,
													dspin tj_124, dspin tj_245, dspin tj_234, dspin tj_235, dspin tj_345,
													dspin tx, dspin tj_125, dspin tj_123, dspin tj_134, dspin tj_145, 
													dspin tr_i1[2], dspin tr_i2[2], dspin tr_i3[2],
													dspin tk_1, dspin tk_2);

// to assembly the delta3 amplitute in the 4d BF case
// and compute the angle
// NB one has to precompute all symbols 
double delta3_angles_bf4_creator (struct data_folders* fd, 
								  dspin tj_124, dspin tj_245, dspin tj_234, dspin tj_235, dspin tj_345, dspin tx,
								  dspin tj_125, dspin tj_123, dspin tj_134, dspin tj_145, dspin tj_256, dspin tj_236,
								  dspin tj_126, dspin tj_136, dspin tj_356, dspin tj_156, dspin tj_346, dspin tj_146, dspin tj_456,
								  dspin ti1_min, dspin ti1_max, dspin ti2_min, dspin ti2_max, dspin ti3_min, dspin ti3_max,
								  dspin ti4_min, dspin ti4_max, dspin ti5_min, dspin ti5_max, dspin ti6_min, dspin ti6_max,
								  dspin ti7_min, dspin ti7_max, dspin ti8_min, dspin ti8_max, dspin ti9_min, dspin ti9_max); 


////////////////////////////////////////////////////////////////
// Functions to compute the 4 simplices delta 3 amplitude EPRL
////////////////////////////////////////////////////////////////

// auxiliary tables for the 3 simplices amplitudes
ARRAY_TABLE_INIT(HashTableEPRL_1, 10, double);
ARRAY_TABLE_INIT(HashTableEPRL_2, 10, double);
ARRAY_TABLE_INIT(HashTableEPRL_3, 10, double);

// to assembly the delta3 amplitude in the 4d EPRL case
// NB one has to precompute all symbols boosters and 
// 4 simplices before using it
double delta3_eprl4_creator (struct data_folders* fd, 
                            dspin tj_124, dspin tj_245, dspin tj_234, dspin tj_235, dspin tj_345, dspin tx,
                            dspin tj_125, dspin tj_123, dspin tj_134, dspin tj_145, dspin tj_256, dspin tj_236,
                            dspin tj_126, dspin tj_136, dspin tj_356, dspin tj_156, dspin tj_346, dspin tj_146, dspin tj_456,
                            dspin ti_1, dspin ti_2, dspin ti_3, dspin ti_4, dspin ti_5,
                            dspin ti_6, dspin ti_7, dspin ti_8, dspin ti_9,
                            dspin td, float immirzi);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_DELTA3_BF4_H__*/