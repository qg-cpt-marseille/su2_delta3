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


#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <libconfig.h>   
#include <gsl/gsl_statistics_double.h>         

#include "sl2cfoam.h"
#include "common.h"
#include "error.h"
#include "dbg.h"
#include "utilities.h"
#include "wigxjpf.h"
#include "jsymbols.h"
#include "config.h"
#include "khash.h"
#include "delta3_4d.h"
#include "coherentstates.h"

#define COHERENT_ARRAY 200

// utility to get angles from 
// configuration file
void get_angles (config_t *cf, char *name_config, 
				 double f[9][4], double t[9][4]) {

		
		const config_setting_t *angles;	
		char config_angles[1024];


		for (int a = 0; a <= 8; a ++) {

			char str[12];
			sprintf(str, "%i", a);
			sprintf(config_angles, "%s.angles%s", name_config, str);

			angles = config_lookup(cf, config_angles);
			f[a][0] = config_setting_get_float_elem(angles, 0);
			t[a][0] = config_setting_get_float_elem(angles, 1);
			f[a][1] = config_setting_get_float_elem(angles, 2);
			t[a][1] = config_setting_get_float_elem(angles, 3);
			f[a][2] = config_setting_get_float_elem(angles, 4);
			t[a][2] = config_setting_get_float_elem(angles, 5);
			f[a][3] = config_setting_get_float_elem(angles, 6);
			t[a][3] = config_setting_get_float_elem(angles, 7);
			
		}

	
}


void coherent_hash (kh_HashTableCoherent_t *hdcoherent, dspin tj_1, dspin tj_2, dspin tj_3, dspin tj_4, 
			 	    double f[3][4],double t[3][4],	
			        dspin ti_min, dspin ti_max, dspin num_t,
					int sign1, int sign2, int sign3, int sign4) {
	
	double f1,f2,f3,f4,t1,t2,t3,t4;					
	
	// fix the right convention for
	// orientation at each 4 tetrahedra
	if(sign1 == -1) {
		f1 = M_PI + f[num_t][0];
		t1 = M_PI - t[num_t][0];
	} else {
		f1 = f[num_t][0];
		t1 = t[num_t][0];
	}	

	if(sign2 == -1) {
		f2 = M_PI + f[num_t][1];
		t2 = M_PI - t[num_t][1];
	} else {
		f2 = f[num_t][1];
		t2 = t[num_t][1];
	}			

	if(sign3 == -1) {
		f3 = M_PI + f[num_t][2];
		t3 = M_PI - t[num_t][2];
	} else {
		f3 = f[num_t][2];
		t3 = t[num_t][2];
	}			

	if(sign4 == -1) {
		f4 = M_PI + f[num_t][3];
		t4 = M_PI - t[num_t][3];
	} else {
		f4 = f[num_t][3];
		t4 = t[num_t][3];
	}							

	// compute
	for (dspin ti = ti_min; ti <= ti_max; ti += 2) {

		khint_t k;
		int ret;
		double complex coherentStatei = 0. + 0.*I;
		HashTableCoherent_key_t key = {num_t , ti};
		k = kh_put(HashTableCoherent, hdcoherent, key, &ret);

		if ( ret == 1 ){
	
			coherentStatei = CoherentStateI(ti, tj_1, tj_2, tj_3, tj_4,
											f1, t1, f2, t2, f3, t3, f4, t4,
											sign1, sign2, sign3, sign4);
			kh_value(hdcoherent, k) = coherentStatei;
			k = kh_get(HashTableCoherent, hdcoherent, key);
			//printf("{%i, %17g + I*%17g}, \n", ti, creal(coherentStatei), cimag(coherentStatei));

			if( kh_val(hdcoherent,k) == 0 || kh_val(hdcoherent,k) != coherentStatei){
			kh_val(hdcoherent,k) = coherentStatei;
			}
		}
	}

}

// tool to rearrange spins in array
// use anticlockwise convention
void rearrange_spins (dspin tj[5][4],
					  dspin tj_1, dspin tj_2, dspin tj_3, dspin tj_4, dspin tj_5,
				 	  dspin tj_6, dspin tj_7, dspin tj_8, dspin tj_9, dspin tj_10){
		
		//first spin of tetrahedra 1				   
		tj[0][0] = tj_2;
		tj[0][1] = tj_3;
		tj[0][2] = tj_5;
		tj[0][3] = tj_4;

		// tetrahedra 2
		tj[1][0] = tj_1;
		tj[1][1] = tj_10;
		tj[1][2] = tj_7;
		tj[1][3] = tj_2;

		tj[2][0] = tj_9;
		tj[2][1] = tj_8;
		tj[2][2] = tj_3;
		tj[2][3] = tj_1;

		tj[3][0] = tj_6;
		tj[3][1] = tj_5;
		tj[3][2] = tj_10;
		tj[3][3] = tj_9;

		tj[4][0] = tj_4;
		tj[4][1] = tj_7;
		tj[4][2] = tj_8;
		tj[4][3] = tj_6;

}


void rearrange_spins_boundary_delta3 (dspin tj[3][4],
					                  dspin tj_1, dspin tj_2, dspin tj_3, dspin tj_4, dspin tj_5,
				 	         		  dspin tj_7, dspin tj_8, dspin tj_9, dspin tj_10){
		
		// tetrahedra 1				   
		tj[0][0] = tj_2;
		tj[0][1] = tj_3;
		tj[0][2] = tj_5;
		tj[0][3] = tj_4;

		// tetrahedra 2
		tj[1][0] = tj_1;
		tj[1][1] = tj_10;
		tj[1][2] = tj_7;
		tj[1][3] = tj_2;

		// tetrahedra 3				   
		tj[2][0] = tj_9;
		tj[2][1] = tj_8;
		tj[2][2] = tj_3;
		tj[2][3] = tj_1;

}

void coherent_cut (  const kh_HashTableCoherent_t *h1, dspin *range,
                     dspin two_j5, dspin two_j6, dspin two_j9, dspin two_j10,
                     dspin cutoff_gaussian, dspin num_t){

  dspin two_limi4_1, two_limi4_2;

  two_limi4_1 = max(abs(two_j5-two_j6),abs(two_j10-two_j9));
  two_limi4_2 = min(two_j5+two_j6,two_j10+two_j9);

  if (cutoff_gaussian == 0 || two_j5 < 5 || two_j6 < 5 || 
  							  two_j9 < 5 || two_j10 < 5 ) {
    range[0] = two_limi4_1;
    range[1] = two_limi4_2;
  }
  else
  {

    khint_t cGauss;
    dspin two_range_i4 = two_limi4_2 - two_limi4_1;
    double dataGauss[two_range_i4];

    //#pragma omp parallel for
    for(dspin two_i4 = two_limi4_1; two_i4 <= two_limi4_2; two_i4 += 2){

              HashTableCoherent_key_t cohGauss = {num_t, two_i4};
              cGauss = kh_get(HashTableCoherent, h1, cohGauss);
              dataGauss[two_i4] = cabsl(kh_value(h1, cGauss));
    }

    double standardDeviation =  gsl_stats_sd (dataGauss, 2, two_range_i4/2);
    double standardDeviationErr= standardDeviation/cutoff_gaussian;

    for(dspin two_i4 = two_limi4_1; two_i4 <= two_limi4_2; two_i4 += 2){

      if (dataGauss[two_i4] > standardDeviationErr){

        range[0] = two_i4;
        break;
      }
    }

    for(unsigned int two_i4 = range[0]; two_i4 <= two_limi4_2; two_i4+=2){
      if (standardDeviationErr > dataGauss[two_i4] ){
        range[1] = two_i4-2;
        break;
      }
    }
  }

}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/


int main(int argc, char **argv) {


	if (argc < 1) {
        fprintf(stderr, "Usage: %s [-v] [-s/-S] [tj/dspin] [-x] -o [out folder] [-d] [-c] [-g]\n", argv[0]);
        exit(EXIT_FAILURE);
    }
		  	   
    // verbose flag
    bool verbose;
    verbose = false;

	// symmetry flag
    bool do_sym;
    do_sym = false;

	// filepath 
	char filename[128];
	char filepath[1024];


 	// default directory	
	strcpy(filepath, "./delta3/data_delta3/");
	struct data_folders* fd = get_data_folders(0);

	// config directory	
	char filepath_config[1024];
	char name_config[1024];
	strcpy(filepath_config, "./delta3/config/");

    /////////////////////////////////////////////////////////////////////
    // flag parsing
    /////////////////////////////////////////////////////////////////////

	bool spins_initialized = false;
	bool conf_initialized = false;
	bool cutoff_initialized = false;
    dspin tjs[18];
	dspin tj;
	dspin tx_min, tx_max;

	// we fix a parameter to control the error 
	// performed while cutting the int sums
	int cutoff_gaussian;

    int opt;
    while ((opt = getopt(argc, argv, "vds:S:x:o:c:g:")) != -1) {
        switch (opt) {
        case 'v': verbose = true; break;
		case 'd': do_sym = true; break;
		case 's': // set spin individually
            sscanf(optarg, "%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu,%hu", 
                    &tjs[0],  &tjs[1],  &tjs[2],  &tjs[3], &tjs[4],
		 		    &tjs[5],  &tjs[6],  &tjs[7],  &tjs[8], &tjs[9],
				    &tjs[10], &tjs[11], &tjs[12], &tjs[13],
			        &tjs[14], &tjs[15], &tjs[16], &tjs[17]);
			spins_initialized = true;
			break;
        case 'S': // set all equal spins
            if (!(tj = (dspin)atoi(optarg))) {
                error("error parsing spin argument");
            };
            int i;
            for (i = 0; i < 18; i++) {
                tjs[i] = tj;   
            }
			spins_initialized = true;
			break;
		case 'x':
 			sscanf(optarg, "%hu,%hu", &tx_min,  &tx_max);
			break;	
		case 'o':
            strcpy(filepath, optarg);
            strcat(filepath, "/");
			break;	
		case 'c':
            strcpy(name_config, optarg);
			conf_initialized = true;
			break;	
		case 'g':
            cutoff_gaussian = (dspin)atoi(optarg);
			cutoff_initialized = true;
			break;	
        default:
            fprintf(stderr, "Usage: %s [-v] [-s/-S] [tj/dspin] -o [out folder] [-d] [-c] [-g] \n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }
	

	if (!(spins_initialized) ) {
		error("spins not set with command line arguments")
	}
	if (!(conf_initialized) ) {
		error("configuration not set with command line arguments")
	}
		if (!(cutoff_initialized) ) {
		error("gaussian cutoff not set with command line arguments")
	}

  	/////////////////////////////////////////////////////////////////////
    // init library and setup spins
    /////////////////////////////////////////////////////////////////////

    // initialize library
    struct sl2cfoam_config libconf;
    libconf.data_folder = "../data_sl2cfoam/";
    libconf.verbosity = 0;
	libconf.coupling = SL2CFOAM_COUPLING_IRREDUCIBLE;
    
    sl2cfoam_init_conf(&libconf);

	// check data folders
	check_data_4simplex(0.0);

	// spins
	dspin tj_124, tj_245, tj_234, tj_235, tj_345,
		  tj_125, tj_123, tj_134, tj_145, tj_256,
		  tj_236, tj_126, tj_136, tj_356,
		  tj_156, tj_346, tj_146, tj_456;

	tj_124 = tjs[0];
	tj_245 = tjs[1];
	tj_234 = tjs[2]; 
	tj_235 = tjs[3]; 
	tj_345 = tjs[4];
	tj_125 = tjs[5]; 
	tj_123 = tjs[6]; 
	tj_134 = tjs[7]; 
	tj_145 = tjs[8]; 
	tj_256 = tjs[9];
	tj_236 = tjs[10]; 
	tj_126 = tjs[11];
	tj_136 = tjs[12]; 
	tj_356 = tjs[13];
	tj_156 = tjs[14]; 
	tj_346 = tjs[15]; 
	tj_146 = tjs[16]; 
	tj_456 = tjs[17];

	sprintf(filename, "d3-dist4-bf_%i.%i.%i.%i.%i.%i.%i.%i.%i.%i.%i.%i.%i.%i.%i.%i.%i.%i_%i.%i_%i.%s.dat", 
			tj_124, tj_245, tj_234, tj_235, tj_345,
		 	tj_125, tj_123, tj_134, tj_145, tj_256,
			tj_236, tj_126, tj_136, tj_356,
			tj_156, tj_346, tj_146, tj_456,
			tx_min, tx_max, cutoff_gaussian,
			name_config);
	
	// open file
	strcpy(filepath, filepath);
	strcat(filepath, filename);

	FILE *file = fopen(filepath, "w");
	if (file == NULL){
		error("error opening file for writing");
	}		

	/////////////////////////////////////////////////////////////////////
    // get coherent states information from config file
    /////////////////////////////////////////////////////////////////////

	//////////////////// initialize configuration file ////////////////////

	config_t cfg, *cf;

	cf = &cfg;
  	config_init(cf);

 	//////////////////// check configuration file ////////////////////
	strcpy(filepath_config, "./delta3/config/d3-config.txt");
	
	if (!config_read_file(cf, filepath_config))
	{
		fprintf(stderr, "%s:%d - %s\n",
				config_error_file(cf),
				config_error_line(cf),
				config_error_text(cf));
		config_destroy(cf);
		return(EXIT_FAILURE);
	}


    //////////////////// initialize configuration angles ////////////////////
    /////////////// phi and Theta for 9 boundary tetrahedra //////////////

	const config_setting_t *angles;
	double f[9][4] ;
	double t[9][4] ;

	// get angles from config file
	get_angles (cf, name_config, f, t );

	// we rearrange spins using arrays to be able
	// to better check what we are doing

	// spins for three tetrahedra in three 4-simplices
	dspin ts_1[3][4], ts_2[3][4], ts_3[3][4];

	rearrange_spins_boundary_delta3 (ts_1,
					 				 tj_124, tj_245, tj_234, tj_235, tj_345,
							         tj_125, tj_123, tj_134, tj_145); 

	rearrange_spins_boundary_delta3 (ts_2,
							  		 tj_256, tj_236, tj_126, tj_136, tj_123,
							   		 tj_356, tj_156, tj_125, tj_235);		   

	rearrange_spins_boundary_delta3 (ts_3,
							  		 tj_346, tj_146, tj_456, tj_145, tj_156,
							 		 tj_134, tj_345, tj_356, tj_136);    


	// now compute boundary range intertwertwiners
	dspin ti1_min = max(abs(ts_1[0][0]-ts_1[0][1]),abs(ts_1[0][2] - ts_1[0][3]));
	dspin ti1_max = min(ts_1[0][0] + ts_1[0][1], ts_1[0][2] + ts_1[0][3]);
	dspin ti2_min = max(abs(ts_1[1][0]-ts_1[1][1]),abs(ts_1[1][2] - ts_1[1][3]));
	dspin ti2_max = min(ts_1[1][0] + ts_1[1][1], ts_1[1][2] + ts_1[1][3]);
	dspin ti3_min = max(abs(ts_1[2][0]-ts_1[2][1]),abs(ts_1[2][2] - ts_1[2][3]));
	dspin ti3_max = min(ts_1[2][0] + ts_1[2][1], ts_1[2][2] + ts_1[2][3]);

	dspin ti4_min = max(abs(ts_2[0][0]-ts_2[0][1]),abs(ts_2[0][2] - ts_2[0][3]));
	dspin ti4_max = min(ts_2[0][0] + ts_2[0][1], ts_2[0][2] + ts_2[0][3]);
	dspin ti5_min = max(abs(ts_2[1][0]-ts_2[1][1]),abs(ts_2[1][2] - ts_2[1][3]));
	dspin ti5_max = min(ts_2[1][0] + ts_2[1][1], ts_2[1][2] + ts_2[1][3]);
	dspin ti6_min = max(abs(ts_2[2][0]-ts_2[2][1]),abs(ts_2[2][2] - ts_2[2][3]));
	dspin ti6_max = min(ts_2[2][0] + ts_2[2][1], ts_2[2][2] + ts_2[2][3]);

	dspin ti7_min = max(abs(ts_3[0][0]-ts_3[0][1]),abs(ts_3[0][2] - ts_3[0][3]));
	dspin ti7_max = min(ts_3[0][0] + ts_3[0][1], ts_3[0][2] + ts_3[0][3]);
	dspin ti8_min = max(abs(ts_3[1][0]-ts_3[1][1]),abs(ts_3[1][2] - ts_3[1][3]));
	dspin ti8_max = min(ts_3[1][0] + ts_3[1][1], ts_3[1][2] + ts_3[1][3]);
	dspin ti9_min = max(abs(ts_3[2][0]-ts_3[2][1]),abs(ts_3[2][2] - ts_3[2][3]));
	dspin ti9_max = min(ts_3[2][0] + ts_3[2][1], ts_3[2][2] + ts_3[2][3]);


	/////////////////////////////////////////////////////////////////////
    // compute coherent states - creation and hashing
    /////////////////////////////////////////////////////////////////////

	khash_t(HashTableCoherent) *hcoherent = kh_init(HashTableCoherent);


	coherent_hash (hcoherent, ts_1[0][0], ts_1[0][1], ts_1[0][2], ts_1[0][3], f, t, ti1_min, ti1_max, 0, -1, -1, -1, -1);
	coherent_hash (hcoherent, ts_1[1][0], ts_1[1][1], ts_1[1][2], ts_1[1][3], f, t, ti2_min, ti2_max, 1, -1, -1, -1, +1);
	coherent_hash (hcoherent, ts_1[2][0], ts_1[2][1], ts_1[2][2], ts_1[2][3], f, t, ti3_min, ti3_max, 2, -1, -1, +1, +1);

	coherent_hash (hcoherent, ts_2[0][0], ts_2[0][1], ts_2[0][2], ts_2[0][3], f, t, ti4_min, ti4_max, 3, -1, +1, +1, +1);
	coherent_hash (hcoherent, ts_2[1][0], ts_2[1][1], ts_2[1][2], ts_2[1][3], f, t, ti5_min, ti5_max, 4, +1, +1, +1, +1);
	coherent_hash (hcoherent, ts_2[2][0], ts_2[2][1], ts_2[2][2], ts_2[2][3], f, t, ti6_min, ti6_max, 5, +1, -1, -1, -1);

	coherent_hash (hcoherent, ts_3[0][0], ts_3[0][1], ts_3[0][2], ts_3[0][3], f, t, ti7_min, ti7_max, 6, -1, -1, +1, +1);
	coherent_hash (hcoherent, ts_3[1][0], ts_3[1][1], ts_3[1][2], ts_3[1][3], f, t, ti8_min, ti8_max, 7, -1, -1, +1, +1);
	coherent_hash (hcoherent, ts_3[2][0], ts_3[2][1], ts_3[2][2], ts_3[2][3], f, t, ti9_min, ti9_max, 8, -1, +1, +1, +1);	  

		
	if (verbose) {
		printf("coherent hash done\n");
	}

	// we have to perform the cut
	// over boundary gaussians

    dspin tr_i1[2], tr_i2[2], tr_i3[2];
	dspin tr_i4[2], tr_i5[2], tr_i6[2];
	dspin tr_i7[2], tr_i8[2], tr_i9[2]; 


	// and we store new ranges for
	// boundary intertwiners
    coherent_cut ( hcoherent, tr_i1,
                   ts_1[0][0], ts_1[0][1], ts_1[0][2], ts_1[0][3], cutoff_gaussian, 0);
    coherent_cut ( hcoherent, tr_i2,
                   ts_1[1][0], ts_1[1][1], ts_1[1][2], ts_1[1][3], cutoff_gaussian, 1);
    coherent_cut ( hcoherent, tr_i3,
                   ts_1[2][0], ts_1[2][1], ts_1[2][2], ts_1[2][3], cutoff_gaussian, 2);

	coherent_cut ( hcoherent, tr_i4,
				ts_2[0][0], ts_2[0][1], ts_2[0][2], ts_2[0][3], cutoff_gaussian, 3);
	coherent_cut ( hcoherent, tr_i5,
				ts_2[1][0], ts_2[1][1], ts_2[1][2], ts_2[1][3], cutoff_gaussian, 4);
	coherent_cut ( hcoherent, tr_i6,
				ts_2[2][0], ts_2[2][1], ts_2[2][2], ts_2[2][3], cutoff_gaussian, 5);				  

	coherent_cut ( hcoherent, tr_i7,
				ts_3[0][0], ts_3[0][1], ts_3[0][2], ts_3[0][3], cutoff_gaussian, 6);
	coherent_cut ( hcoherent, tr_i8,
				ts_3[1][0], ts_3[1][1], ts_3[1][2], ts_3[1][3], cutoff_gaussian, 7);
	coherent_cut ( hcoherent, tr_i9,
				ts_3[2][0], ts_3[2][1], ts_3[2][2], ts_3[2][3], cutoff_gaussian, 8);



	/////////////////////////////////////////////////////////////////////
    // hash 6j symbols
    /////////////////////////////////////////////////////////////////////		

	// for every internal spin
	// we hash all needed 15j symbols.
	#pragma omp parallel  
	{ 

		init_wigxjpf_thread();

		// we can parallelize, 
		// data are saved labeled by x	
		#pragma omp for	
		for (dspin tx = tx_min; tx <= tx_max; tx +=2) {

			for (dspin ti_1 = tr_i1[0]; ti_1 <= tr_i1[1]; ti_1 +=2) { 
				
				char filename[1024];
				char path_6j[1024];
				char path_ampl[1024];
				
				// build path for 6j symbols
				sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.6j", tj_245, tj_234, tj_345, tj_235, tx, ti_1);
				strcpy(path_6j, fd->foursimp_hashtable6j_6j);
				strcat(path_6j, filename);
				
				// build path for amp
				sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_245, tj_234, tj_345, tj_235, tx, ti_1);
				strcpy(path_ampl, fd->foursimp_bf_hashtableampl_ampl);
				strcat(path_ampl, filename);

				//////////////////// Hash Table initialization ///////////////////

				// check for previously computed tables.
				khash_t(HashTableJ6) *h = NULL;

				//if (file_exist(path_6j) != 0 || file_exist(path_ampl) != 0) continue;
				if (file_exist(path_6j) != 0) continue;

				if (h == NULL) {
					h = kh_init(HashTableJ6);
				}

			for (dspin ti_2 = tr_i2[0]; ti_2 <= tr_i2[1]; ti_2 +=2) { 
			for (dspin ti_3 = tr_i3[0]; ti_3 <= tr_i3[1]; ti_3 +=2) { 

				// first 4-simplex 12345
				// first internal intertwiner
				// j134, j145, j345

				j15_symbol_hash_delta3_2(h,
									  	 tj_124, tj_245, tj_234, tj_235, tj_345,
									     tx    , tj_125, tj_123, tj_134, tj_145, 
									     ti_1, ti_2, ti_3);
			} //i3
			} //i2

			// write hashtables to disk
    		kh_write(HashTableJ6, h, path_6j);
    		// free memory
    		kh_destroy(HashTableJ6, h);

			} //i1
		} // x

		if (do_sym == false) {

			#pragma omp for	
			for (dspin tx = tx_min; tx <= tx_max; tx +=2) {

				for (dspin ti_4 = tr_i4[0]; ti_4 <= tr_i4[1]; ti_4 +=2) { 
					
					char filename[1024];
					char path_6j[1024];
					char path_ampl[1024];
					
					// build path for 6j symbols
					sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.6j",tj_236, tj_126, tj_123, tj_136, tx, ti_4);
					strcpy(path_6j, fd->foursimp_hashtable6j_6j);
					strcat(path_6j, filename);

					// build path for amp
					sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_236, tj_126, tj_123, tj_136, tx, ti_4);
					strcpy(path_ampl, fd->foursimp_bf_hashtableampl_ampl);
					strcat(path_ampl, filename);	

					//////////////////// Hash Table initialization ///////////////////

					// check for previously computed tables.
					khash_t(HashTableJ6) *h = NULL;
				
					//if (file_exist(path_6j) != 0 || file_exist(path_ampl) != 0) continue;
					if (file_exist(path_6j) != 0 ) continue;
					
					if (h == NULL) {
						h = kh_init(HashTableJ6);
					}

				for (dspin ti_5 = tr_i5[0]; ti_5 <= tr_i5[1]; ti_5 +=2) { 
				for (dspin ti_6 = tr_i6[0]; ti_6 <= tr_i6[1]; ti_6 +=2) { 

					// second 4-simplex 12356
					// second internal intertwiner
					// j125, j235, j123

					j15_symbol_hash_delta3_2(h,
										     tj_256, tj_236, tj_126, tj_136, tj_123,
										     tx    , tj_356, tj_156, tj_125, tj_235,
										     ti_4, ti_5, ti_6);	
				} //i6
				} //i5

					// write hashtables to disk
					kh_write(HashTableJ6, h, path_6j);
					// free memory
					kh_destroy(HashTableJ6, h);

				} //i4
			} //x

			#pragma omp for	
			for (dspin tx = tx_min; tx <= tx_max; tx +=2) {
				for (dspin ti_7 = tr_i7[0]; ti_7 <= tr_i7[1]; ti_7 +=2) { 

					char filename[1024];
					char path_6j[1024];
					char path_ampl[1024];

					// build path for 6j symbols
					sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.6j",tj_146, tj_456, tj_156, tj_145, tx, ti_7);
					strcpy(path_6j, fd->foursimp_hashtable6j_6j);
					strcat(path_6j, filename);

					// build path for amp
					sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_146, tj_456, tj_156, tj_145, tx, ti_7);
					strcpy(path_ampl, fd->foursimp_bf_hashtableampl_ampl);
					strcat(path_ampl, filename);

					//////////////////// Hash Table initialization ///////////////////

					// check for previously computed tables.
					khash_t(HashTableJ6) *h = NULL;

					//if (file_exist(path_6j) != 0 || file_exist(path_ampl) != 0 ) continue;
					if (file_exist(path_6j) != 0 ) continue;
					
					if (h == NULL) {
						h = kh_init(HashTableJ6);
					}

				for (dspin ti_8 = tr_i8[0]; ti_8 <= tr_i8[1]; ti_8 +=2) { 
				for (dspin ti_9 = tr_i9[0]; ti_9 <= tr_i9[1]; ti_9 +=2) { 


					// third 4-simplex 13456
					// third internal intertwiner
					// j356, j136, j156

					j15_symbol_hash_delta3_2(h,
										     tj_346, tj_146, tj_456, tj_145, tj_156,
										     tx    , tj_134, tj_345, tj_356, tj_136,
										     ti_7, ti_8, ti_9);

				} // i9
				} // i8

				// write hashtables to disk
				kh_write(HashTableJ6, h, path_6j);
				// free memory
				kh_destroy(HashTableJ6, h);

				} // i7
			} // x
		}

	} // omp parallel
		
	if (verbose) {
		printf("15j hash done\n");
	}	

	/////////////////////////////////////////////////////////////////////
	// hash BF amplitudes
	/////////////////////////////////////////////////////////////////////
	#pragma omp parallel  
	{ 

		#pragma omp for	
		for (dspin tx = tx_min; tx <= tx_max; tx +=2) {
				
			for (dspin ti_1 = tr_i1[0]; ti_1 <= tr_i1[1]; ti_1 +=2) { 

				// build path for 6js
				char filename[1024];
				char path_6j[1024];
				sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.6j", tj_245, tj_234, tj_345, tj_235, tx, ti_1);
				strcpy(path_6j, fd->foursimp_hashtable6j_6j);
				strcat(path_6j, filename);

				// load hash tables for 6Js
				khash_t(HashTableJ6) *h = kh_load(HashTableJ6, path_6j);

				char path_ampl[1024];

				// build path for amplitude values
				sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_245, tj_234, tj_345, tj_235, tx, ti_1);
				strcpy(path_ampl, fd->foursimp_bf_hashtableampl_ampl);
				strcat(path_ampl, filename);

				// initialize or load hash table of amplitudes on disk
				kh_HashTableEPRL_t *h3 = NULL;
				
				if (file_exist(path_ampl) != 0 ) {
					continue;
					//h3 = kh_load(HashTableEPRL,path_ampl);
				} else {

				// table not found, initialize it
				h3 = kh_init(HashTableEPRL);
				}
					

			for (dspin ti_2 = tr_i2[0]; ti_2 <= tr_i2[1]; ti_2 +=2) { 
			for (dspin ti_3 = tr_i3[0]; ti_3 <= tr_i3[1]; ti_3 +=2) { 

				// first 4-simplex 12345
				sl2cfoam_hash_four_ampl_BF_delta3_2(h, h3,
													tj_124, tj_245, tj_234, tj_235, tj_345,
													tx, tj_125, tj_123, tj_134, tj_145, 
												  	ti_1,   ti_2, ti_3);

			} // i3
			} // i2

				// write hash table for this x,i1 to disk
        		kh_write(HashTableEPRL, h3, path_ampl);
        		// free memory
        		kh_destroy(HashTableJ6, h);
        		kh_destroy(HashTableEPRL, h3);

			} // i1
		} // x

		if (do_sym == false) {
		
			#pragma omp for	
			for (dspin tx = tx_min; tx <= tx_max; tx +=2) {

				for (dspin ti_4 = tr_i4[0]; ti_4 <= tr_i4[1]; ti_4 +=2) { 
					
					// build path for 6js
					char filename[1024];
					char path_6j[1024];
					sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.6j", tj_236, tj_126, tj_123, tj_136, tx, ti_4);
					strcpy(path_6j, fd->foursimp_hashtable6j_6j);
					strcat(path_6j, filename);

					// load hash tables for 6Js
					khash_t(HashTableJ6) *h = kh_load(HashTableJ6, path_6j);

					char path_ampl[1024];

					// build path for amplitude values
					sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_236, tj_126, tj_123, tj_136, tx, ti_4);
					strcpy(path_ampl, fd->foursimp_bf_hashtableampl_ampl);
					strcat(path_ampl, filename);

					// initialize or load hash table of amplitudes on disk
					kh_HashTableEPRL_t *h3 = NULL;


					if (file_exist(path_ampl) != 0 ) {
						//h3 = kh_load(HashTableEPRL,path_ampl);
						continue;
					} else {

					// table not found, initialize it
					h3 = kh_init(HashTableEPRL);
					}
						
				for (dspin ti_5 = tr_i5[0]; ti_5 <= tr_i5[1]; ti_5 +=2) { 
				for (dspin ti_6 = tr_i6[0]; ti_6 <= tr_i6[1]; ti_6 +=2) { 

					// second 4-simplex 12356
					sl2cfoam_hash_four_ampl_BF_delta3_2(h, h3,
													  	tj_256, tj_236, tj_126, tj_136, tj_123,
													  	tx, tj_356, tj_156, tj_125, tj_235,
													 	ti_4, ti_5, ti_6);

				} // i6
				} // i5
					// write hash table for this x,i1 to disk
					kh_write(HashTableEPRL, h3, path_ampl);
					// free memory
					kh_destroy(HashTableJ6, h);
					kh_destroy(HashTableEPRL, h3);
				} // i4	
			} //x

			#pragma omp for	
			for (dspin tx = tx_min; tx <= tx_max; tx +=2) {

				for (dspin ti_7 = tr_i7[0]; ti_7 <= tr_i7[1]; ti_7 +=2) { 

					// build path for 6js
					char filename[1024];
					char path_6j[1024];
					sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.6j", tj_146, tj_456, tj_156, tj_145, tx, ti_7);
					strcpy(path_6j, fd->foursimp_hashtable6j_6j);
					strcat(path_6j, filename);

					// load hash tables for 6Js
					khash_t(HashTableJ6) *h = kh_load(HashTableJ6, path_6j);

					char path_ampl[1024];

					// build path for amplitude values
					sprintf(filename, "%i.%i.%i.%i.%i_%i_d3.bf", tj_146, tj_456, tj_156, tj_145, tx, ti_7);
					strcpy(path_ampl, fd->foursimp_bf_hashtableampl_ampl);
					strcat(path_ampl, filename);

					// initialize or load hash table of amplitudes on disk
					kh_HashTableEPRL_t *h3 = NULL;
					
					if (file_exist(path_ampl) != 0 ) {
						continue;
						//h3 = kh_load(HashTableEPRL,path_ampl);
					} else {

					// table not found, initialize it
					h3 = kh_init(HashTableEPRL);
					}
						

				for (dspin ti_8 = tr_i8[0]; ti_8 <= tr_i8[1]; ti_8 +=2) { 
				for (dspin ti_9 = tr_i9[0]; ti_9 <= tr_i9[1]; ti_9 +=2) { 

					// third 4-simplex 13456
					sl2cfoam_hash_four_ampl_BF_delta3_2(h, h3,
													    tj_346, tj_146, tj_456, tj_145, tj_156,
													    tx, tj_134, tj_345, tj_356, tj_136,
													    ti_7, ti_8, ti_9);

				} // i9
				} // i8

					// write hash table for this x,i1 to disk
					kh_write(HashTableEPRL, h3, path_ampl);
					// free memory
					kh_destroy(HashTableJ6, h);
					kh_destroy(HashTableEPRL, h3);

				} // i7	
			} // x
		}

		clear_wigxjpf_thread();
	
	}// omp parallel

	if (verbose) {
		printf("amp hash done\n");
	}

	/////////////////////////////////////////////////////////////////////
	// hash Coherent BF amplitudes
	/////////////////////////////////////////////////////////////////////

	// initialize table
	khash_t(HashTableCoherentBF) *hbfcoherent = NULL;

	if (hbfcoherent == NULL) {
		hbfcoherent = kh_init(HashTableCoherentBF);
	}
	else{
		error("problem in initializing hash table");
	}

	double complex ***coherent_amp_value = (double complex***) malloc(COHERENT_ARRAY*sizeof(double complex**));
            
	// initialize double complex array
	int i,j,k;

	for (i = 0; i < COHERENT_ARRAY; i++) {
		coherent_amp_value[i] = (double complex**) malloc(COHERENT_ARRAY*sizeof(double complex*));
		for (j = 0; j < COHERENT_ARRAY; j++) {
			coherent_amp_value[i][j] = (double complex*) malloc(COHERENT_ARRAY*sizeof(double complex));
		}
	}

	for (i = 0; i < COHERENT_ARRAY; i++) {
		for (j = 0; j < COHERENT_ARRAY; j++) {
			for (k = 0; k < COHERENT_ARRAY; k++) {
				coherent_amp_value[i][j][k] = 0;

			}
		}
    }

	// start hashing the first 4 simplex
	// with 3 coherent states
	#pragma omp parallel  
	{ 
		
	#pragma omp for	
	for (dspin tx = tx_min; tx <= tx_max; tx +=2) {
		
		dspin tk1_min, tk1_max;

		tk1_min = (dspin) max(abs(tj_134-tj_145), abs(tx-tj_345));
		tk1_max = (dspin) min(tj_134 + tj_145, tx + tj_345);

		dspin tk2_min, tk2_max;

		tk2_min = (dspin) max(abs(tj_125-tj_235), abs(tx-tj_123));
		tk2_max = (dspin) min(tj_125 + tj_235, tx + tj_123);
	

		// for every combination of internal int
		// we resum the boundary intertwiners
		// weightened by coherent states
		for (dspin tk_1 = tk1_min; tk_1 <= tk1_max; tk_1 += 2) {
		for (dspin tk_2 = tk2_min; tk_2 <= tk2_max; tk_2 += 2) {

		// first 4-simplex 12345
		coherent_amp_value[tx][tk_1][tk_2] = delta3_simplex_coherent_bf4_hashing (fd, hcoherent, hbfcoherent,
																				  tj_124,  tj_245,  tj_234,  tj_235,  tj_345,
																				  tx,  tj_125,  tj_123,  tj_134,  tj_145, 
																				  tr_i1,  tr_i2, tr_i3, tk_1, tk_2) ;	

		} // k2
		} // k1
	} //x

	} // omp parallel

	// now we store the date on an hash table
	for (dspin tx = tx_min; tx <= tx_max; tx +=2) {
			dspin tk1_min, tk1_max;

		tk1_min = (dspin) max(abs(tj_134-tj_145), abs(tx-tj_345));
		tk1_max = (dspin) min(tj_134 + tj_145, tx + tj_345);

		dspin tk2_min, tk2_max;

		tk2_min = (dspin) max(abs(tj_125-tj_235), abs(tx-tj_123));
		tk2_max = (dspin) min(tj_125 + tj_235, tx + tj_123);
	
		for (dspin tk_1 = tk1_min; tk_1 <= tk1_max; tk_1 += 2) {
		for (dspin tk_2 = tk2_min; tk_2 <= tk2_max; tk_2 += 2) {
		
		// put the value in the hash table
		khint_t k;
		int ret;
		HashTableCoherentBF_key_t key = {1 , tx, tk_1, tk_2};
		k = kh_put(HashTableCoherentBF, hbfcoherent, key, &ret);
		//printf("A %i %i %i %17g %17g\n", tx, tk_1, tk_2, creal(coherent_simplex),cimag(coherent_simplex));

		if ( ret == 1 ){

			kh_value(hbfcoherent, k) = coherent_amp_value[tx][tk_1][tk_2];
			k = kh_get(HashTableCoherentBF, hbfcoherent, key);

			if( kh_val(hbfcoherent,k) == 0 || kh_val(hbfcoherent,k) != coherent_amp_value[tx][tk_1][tk_2]){
			kh_val(hbfcoherent,k) = coherent_amp_value[tx][tk_1][tk_2];

			}
		}
	
		} // k2
		} //k1

	} // x

	// reset array to zero
	for (i = 0; i < COHERENT_ARRAY; i++) {
		for (j = 0; j < COHERENT_ARRAY; j++) {
			for (k = 0; k < COHERENT_ARRAY; k++) {
				coherent_amp_value[i][j][k] = 0;

			}
		}
    }

	#pragma omp parallel  
	{ 

	#pragma omp for	
	for (dspin tx = tx_min; tx <= tx_max; tx +=2) {
		
		dspin tk2_min, tk2_max;

		tk2_min = (dspin) max(abs(tj_125-tj_235), abs(tx-tj_123));
		tk2_max = (dspin) min(tj_125 + tj_235, tx + tj_123);

		dspin tk3_min, tk3_max;

		tk3_min = (dspin) max(abs(tj_356-tj_136), abs(tx-tj_156));
		tk3_max = (dspin) min(tj_356 + tj_136, tx + tj_156); 


		// for every combination of internal int
		// we resum the boundary intertwiners
		// weightened by coherent states
		for (dspin tk_2 = tk2_min; tk_2 <= tk2_max; tk_2 += 2) {
		for (dspin tk_3 = tk3_min; tk_3 <= tk3_max; tk_3 += 2) {

			// second 4-simplex 12356
			coherent_amp_value[tx][tk_2][tk_3] = delta3_simplex_coherent_bf4_hashing(fd, hcoherent, hbfcoherent,
																					tj_256, tj_236, tj_126, tj_136, tj_123,
																					tx, tj_356, tj_156, tj_125, tj_235,
																					tr_i4,  tr_i5, tr_i6, tk_2, tk_3);
		
		} // k3
		} // k2

	}  // x

	} // omp parallel

	// now we store the date on an hash table
	for (dspin tx = tx_min; tx <= tx_max; tx +=2) {
		
		dspin tk2_min, tk2_max;

		tk2_min = (dspin) max(abs(tj_125-tj_235), abs(tx-tj_123));
		tk2_max = (dspin) min(tj_125 + tj_235, tx + tj_123);

		dspin tk3_min, tk3_max;

		tk3_min = (dspin) max(abs(tj_356-tj_136), abs(tx-tj_156));
		tk3_max = (dspin) min(tj_356 + tj_136, tx + tj_156); 

		for (dspin tk_2 = tk2_min; tk_2 <= tk2_max; tk_2 += 2) {
		for (dspin tk_3 = tk3_min; tk_3 <= tk3_max; tk_3 += 2) {
		
			// put the value in the hash table
			khint_t k;
			int ret;
			HashTableCoherentBF_key_t key = {2 , tx, tk_2, tk_3};
			k = kh_put(HashTableCoherentBF, hbfcoherent, key, &ret);

			if ( ret == 1 ){

				kh_value(hbfcoherent, k) = coherent_amp_value[tx][tk_2][tk_3];
				k = kh_get(HashTableCoherentBF, hbfcoherent, key);

				if( kh_val(hbfcoherent,k) == 0 || kh_val(hbfcoherent,k) != coherent_amp_value[tx][tk_2][tk_3]){
				kh_val(hbfcoherent,k) = coherent_amp_value[tx][tk_2][tk_3];

				}
			}
	
		} // k3
		} // k2

	} // x

	// reset array to zero
	for (i = 0; i < COHERENT_ARRAY; i++) {
		for (j = 0; j < COHERENT_ARRAY; j++) {
			for (k = 0; k < COHERENT_ARRAY; k++) {
				coherent_amp_value[i][j][k] = 0;
			}
		}
    }


	#pragma omp parallel  
	{ 

	#pragma omp for	
	for (dspin tx = tx_min; tx <= tx_max; tx +=2) {

		dspin tk3_min, tk3_max;

		tk3_min = (dspin) max(abs(tj_356-tj_136), abs(tx-tj_156));
		tk3_max = (dspin) min(tj_356 + tj_136, tx + tj_156); 

		dspin tk1_min, tk1_max;

		tk1_min = (dspin) max(abs(tj_134-tj_145), abs(tx-tj_345));
		tk1_max = (dspin) min(tj_134 + tj_145, tx + tj_345);

		// for every combination of internal int
		// we resum the boundary intertwiners
		// weightened by coherent states
		for (dspin tk_3 = tk3_min; tk_3 <= tk3_max; tk_3 += 2) {
		for (dspin tk_1 = tk1_min; tk_1 <= tk1_max; tk_1 += 2) {

			// third 4-simplex 13456
			coherent_amp_value[tx][tk_3][tk_1] = delta3_simplex_coherent_bf4_hashing(fd, hcoherent, hbfcoherent,
																					tj_346, tj_146, tj_456, tj_145, tj_156,
																					tx, tj_134, tj_345, tj_356, tj_136,
																					tr_i7, tr_i8, tr_i9, tk_3, tk_1);			
		} // k3
		} // k1								
	
	} // x

	} // omp parallel

	// now we store the date on an hash table
	for (dspin tx = tx_min; tx <= tx_max; tx +=2) {

		dspin tk3_min, tk3_max;

		tk3_min = (dspin) max(abs(tj_356-tj_136), abs(tx-tj_156));
		tk3_max = (dspin) min(tj_356 + tj_136, tx + tj_156); 

		dspin tk1_min, tk1_max;

		tk1_min = (dspin) max(abs(tj_134-tj_145), abs(tx-tj_345));
		tk1_max = (dspin) min(tj_134 + tj_145, tx + tj_345);

		for (dspin tk_3 = tk3_min; tk_3 <= tk3_max; tk_3 += 2) {
		for (dspin tk_1 = tk1_min; tk_1 <= tk1_max; tk_1 += 2) {
			
			// put the value in the hash table
			khint_t k;
			int ret;
			HashTableCoherentBF_key_t key = {3, tx, tk_3, tk_1};
			k = kh_put(HashTableCoherentBF, hbfcoherent, key, &ret);

			if ( ret == 1 ){

				kh_value(hbfcoherent, k) = coherent_amp_value[tx][tk_3][tk_1];
				k = kh_get(HashTableCoherentBF, hbfcoherent, key);

				if( kh_val(hbfcoherent,k) == 0 || kh_val(hbfcoherent,k) != coherent_amp_value[tx][tk_3][tk_1]){
				kh_val(hbfcoherent,k) = coherent_amp_value[tx][tk_3][tk_1];

				}
			}

		} // k3
		} // k1

	} // x


	// clear array
	for (i = 0; i < COHERENT_ARRAY; i++) {
		for (j = 0; j < COHERENT_ARRAY; j++) {
			free(coherent_amp_value[i][j]);
		}
	}

	for (i = 0; i < COHERENT_ARRAY; i++) {
		free(coherent_amp_value[i]);
	}

	free(coherent_amp_value);


	if (verbose) {
			printf("coherent amp hash done\n");
	}

	/////////////////////////////////////////////////////////////////////
	// constructing the amplitude
	/////////////////////////////////////////////////////////////////////

	// buffer array to safely parallelize
	size_t nx = (tx_max - tx_min) / 2 + 1;
	double* buf_re = calloc(nx, sizeof(double));
	double* buf_im = calloc(nx, sizeof(double));

	int nthreads;
	size_t x_per_thread;
	size_t x_last;
	int filled = 0;

	// threshold to avoid zeros
	double threshold = 1e-40;

	#pragma omp parallel reduction (+:filled)
	{

	sl2cfoam_init_thread();

	#pragma omp single
	{

		// get number of threads
		nthreads = omp_get_num_threads();
		// get the block of xs per each thread
		// there will be x missing to be filled at the end (division with remainder...)
		
		// if we have les terms then core
		// execute one computation per core
		if (nx <= nthreads) {

			x_per_thread = 1;
			x_last = 0;
			nthreads = nx ;

		} else {

			x_per_thread = nx / nthreads;		
			x_last = nx % nthreads;

		}
	} //omp single

	#pragma omp for
	for (int thread_id = 0; thread_id < nthreads; thread_id++) {

		dspin two_x_min_thread = max(tx_min + 2*thread_id*x_per_thread, tx_min);
		dspin two_x_max_thread = min(two_x_min_thread + 2*x_per_thread, tx_max + 2);

		int bufi = x_per_thread * thread_id;
		
		for (dspin tx = two_x_min_thread; tx < two_x_max_thread; tx += 2) {	
	
			dspin tk1_min, tk1_max;

			tk1_min = (dspin) max(abs(tj_134-tj_145), abs(tx-tj_345));
			tk1_max = (dspin) min(tj_134 + tj_145, tx + tj_345);

			dspin tk2_min, tk2_max;

			tk2_min = (dspin) max(abs(tj_125-tj_235), abs(tx-tj_123));
			tk2_max = (dspin) min(tj_125 + tj_235, tx + tj_123);

			dspin tk3_min, tk3_max;

			tk3_min = (dspin) max(abs(tj_356-tj_136), abs(tx-tj_156));
			tk3_max = (dspin) min(tj_356 + tj_136, tx + tj_156); 

			double complex delta3_x = 0.0 + 0.0*I;
			double complex err = 0.0 + 0.0*I;
			double complex simplex_1 = 0.0 + 0.0*I;
			double complex simplex_2 = 0.0 + 0.0*I;
			double complex simplex_3 = 0.0 + 0.0*I;

			for (dspin tk_1 = tk1_min; tk_1 <= tk1_max; tk_1 += 2) {
			for (dspin tk_2 = tk2_min; tk_2 <= tk2_max; tk_2 += 2) {
				
				// get first coherent 4 simplex from table
				HashTableCoherentBF_key_t keyCbf_1 = {1, tx, tk_1, tk_2};
				khint_t s_1 = kh_get(HashTableCoherentBF, hbfcoherent, keyCbf_1);

				if (s_1 != kh_end(hbfcoherent)) {		
					simplex_1 = kh_val(hbfcoherent, s_1);
				} else {
					error("value for coherent simplex not found");
				}

				if(cabs(simplex_1) < threshold) continue;

			for (dspin tk_3 = tk3_min; tk_3 <= tk3_max; tk_3 += 2) {

				// get second coherent 4 simplex from table

				HashTableCoherentBF_key_t keyCbf_2 = {2, tx, tk_2, tk_3};
				khint_t s_2 = kh_get(HashTableCoherentBF, hbfcoherent, keyCbf_2);

				if (s_2 != kh_end(hbfcoherent)) {		
					simplex_2 = kh_val(hbfcoherent, s_2);
				} else {
					error("value for coherent simplex not found");
				}

				if(cabs(simplex_2) < threshold) continue;		
			
				// get first coherent 4 simplex from table
				HashTableCoherentBF_key_t keyCbf_3 = {3, tx, tk_3, tk_1};
				khint_t s_3 = kh_get(HashTableCoherentBF, hbfcoherent, keyCbf_3);

				if (s_3 != kh_end(hbfcoherent)) {		
					simplex_3 = kh_val(hbfcoherent, s_3);
				} else {
					error("value for coherent simplex not found");
				}

				if(cabs(simplex_3) < threshold) continue;
				
				double complex d3value = 0.0 + 0.0*I;

				d3value = mpfr_complex_mul (simplex_1, simplex_2);
				d3value = mpfr_complex_mul (d3value, simplex_3);
				d3value = d3value*d(tx)*real_negpow(3*(tk_1 + tk_2 + tk_3) + 3*tx)/(d(tk_1)*d(tk_2)*d(tk_3));

				// assembly everything at arbitraty precision
				mpfr_t simplexMPFR_Re, simplexMPFR_Im;

				mpfr_init_set_d(simplexMPFR_Re, creal(d3value), MPFR_RNDN);
				mpfr_init_set_d(simplexMPFR_Im, cimag(d3value), MPFR_RNDN);

				compsum_mpfr_complex (&err, &delta3_x, simplexMPFR_Re, simplexMPFR_Im);	

				mpfr_clears( simplexMPFR_Re, simplexMPFR_Im, NULL);		
			}
			}
			}

			buf_re[bufi] = creal(delta3_x);
			buf_im[bufi] = cimag(delta3_x);
			bufi++;
			filled++;
		
		} // x

	} //threads

	int tid = omp_get_thread_num();
	if (tid != 0) {
		// not in master thread
		wig_temp_free();
	}


	
	} // omp parallel

	// last xs
	if (x_last > 0) {
		
		for (dspin tx = tx_max - 2*(x_last-1); tx <= tx_max; tx += 2) {	

			dspin tk1_min, tk1_max;

			tk1_min = (dspin) max(abs(tj_134-tj_145), abs(tx-tj_345));
			tk1_max = (dspin) min(tj_134 + tj_145, tx + tj_345);

			dspin tk2_min, tk2_max;

			tk2_min = (dspin) max(abs(tj_125-tj_235), abs(tx-tj_123));
			tk2_max = (dspin) min(tj_125 + tj_235, tx + tj_123);

			dspin tk3_min, tk3_max;

			tk3_min = (dspin) max(abs(tj_356-tj_136), abs(tx-tj_156));
			tk3_max = (dspin) min(tj_356 + tj_136, tx + tj_156); 

			double complex delta3_x = 0.0 + 0.0*I;
			double complex err = 0.0 + 0.0*I;
			double complex simplex_1 = 0.0 + 0.0*I;
			double complex simplex_2 = 0.0 + 0.0*I;
			double complex simplex_3 = 0.0 + 0.0*I;

			for (dspin tk_1 = tk1_min; tk_1 <= tk1_max; tk_1 += 2) {
			for (dspin tk_2 = tk2_min; tk_2 <= tk2_max; tk_2 += 2) {
				
				// get first coherent 4 simplex from table
				HashTableCoherentBF_key_t keyCbf_1 = {1, tx, tk_1, tk_2};
				khint_t s_1 = kh_get(HashTableCoherentBF, hbfcoherent, keyCbf_1);

				if (s_1 != kh_end(hbfcoherent)) {		
					simplex_1 = kh_val(hbfcoherent, s_1);
				} else {
					error("value for coherent simplex not found");
				}

				if(cabs(simplex_1) < threshold) continue;

			for (dspin tk_3 = tk3_min; tk_3 <= tk3_max; tk_3 += 2) {

				// get second coherent 4 simplex from table

				HashTableCoherentBF_key_t keyCbf_2 = {2, tx, tk_2, tk_3};
				khint_t s_2 = kh_get(HashTableCoherentBF, hbfcoherent, keyCbf_2);

				if (s_2 != kh_end(hbfcoherent)) {		
					simplex_2 = kh_val(hbfcoherent, s_2);
				} else {
					error("value for coherent simplex not found");
				}

				if(cabs(simplex_2) < threshold) continue;		
			
				// get first coherent 4 simplex from table
				HashTableCoherentBF_key_t keyCbf_3 = {3, tx, tk_3, tk_1};
				khint_t s_3 = kh_get(HashTableCoherentBF, hbfcoherent, keyCbf_3);

				if (s_3 != kh_end(hbfcoherent)) {		
					simplex_3 = kh_val(hbfcoherent, s_3);
				} else {
					error("value for coherent simplex not found");
				}

				if(cabs(simplex_3) < threshold) continue;
				
				double complex d3value = 0.0 + 0.0*I;

				d3value = mpfr_complex_mul (simplex_1, simplex_2);
				d3value = mpfr_complex_mul (d3value, simplex_3);
				d3value = d3value*d(tx)*real_negpow(3*(tk_1 + tk_2 + tk_3) + 3*tx)/(d(tk_1)*d(tk_2)*d(tk_3));

				// assembly everything at arbitraty precision
				mpfr_t simplexMPFR_Re, simplexMPFR_Im;

				mpfr_init_set_d(simplexMPFR_Re, creal(d3value), MPFR_RNDN);
				mpfr_init_set_d(simplexMPFR_Im, cimag(d3value), MPFR_RNDN);

				compsum_mpfr_complex (&err, &delta3_x, simplexMPFR_Re, simplexMPFR_Im);	

				mpfr_clears( simplexMPFR_Re, simplexMPFR_Im, NULL);		
			}
			}
			}

			buf_re[filled] = creal(delta3_x);
			buf_im[filled] = cimag(delta3_x);
			filled++;
	
		} //tx

	} //last xs

	// destroy table for coherent states
	kh_destroy(HashTableCoherent, hcoherent);
	kh_destroy(HashTableCoherentBF, hbfcoherent);

	// check that we computed all x required
	if (filled != nx) {
		error("computed %i values but %i were expected", filled, (int)nx);
	}

	// print loop
	dspin two_x = tx_min;
	for (int i = 0; i < nx; i++) {

		fprintf(file, "%i %.3e %.3e\n", two_x, buf_re[i], buf_im[i]);
		two_x += 2;

	}
	
	fclose(file);
	//free(buf_re);
	//free(buf_im);

	if (verbose) {
		printf(" done.\n");
	}

	sl2cfoam_free();

	return EXIT_SUCCESS;

}
