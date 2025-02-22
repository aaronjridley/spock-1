// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at

//   http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#include <mpi.h>
#include "propagator.h"
#include "options.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_sf_legendre.h"
#include "f2c.h"

// computes orbit average of parameters. so far: arg perigee.
// in the main output file, the average is reported as 9999.999999 unless eaxctly one full orbit has been travelled, in which case the average fo the parameter over this orbit is reported and the string "ORB" is indluced at the end of the output file
int orbit_ave(  SPACECRAFT_T *SC,
			    CONSTELLATION_T *CONSTELLATION,
			    int iProc,
			    int iDebugLevel,
		double previous_an_to_sc
		){


  // to compute orbit average orbital elements
  int jump360 = 0;
  if ( fabs(SC->OE.an_to_sc - previous_an_to_sc) > 350/RAD2DEG){ // this condition occurs when th sc crosses the AN
      jump360 = 1;
    }
  

  if ( ( ( (SC->OE.an_to_sc >= SC->OE.initial_an_to_sc) && ( previous_an_to_sc < SC->OE.initial_an_to_sc) ) || ((SC->OE.an_to_sc <= SC->OE.initial_an_to_sc) && ( previous_an_to_sc > SC->OE.initial_an_to_sc))) && (jump360 == 0) ) { // the sc has traveled a multiple of one orbit around the Earth since the start of the simu


    //      etprint(SC->et,"");
	SC->OE.w_ave = SC->OE.w_ave_temp / SC->OE.ave_increm;  
	SC->OE.w_ave = SC->OE.w_ave_temp / SC->OE.ave_increm;  

	SC->OE.sma_ave = SC->OE.sma_ave_temp / SC->OE.ave_increm;  
	SC->OE.sma_ave = SC->OE.sma_ave_temp / SC->OE.ave_increm;  

	SC->OE.ecc_ave = SC->OE.ecc_ave_temp / SC->OE.ave_increm;  
	SC->OE.ecc_ave = SC->OE.ecc_ave_temp / SC->OE.ave_increm;  

/*  printf("w %f\n", SC->OE.w_ave * 180/M_PI);  */
/*  printf("w %f\n", SC->OE.sma_ave );  */
      SC->OE.w_ave_temp = SC->OE.w;
      SC->OE.sma_ave_temp = SC->OE.sma;
      SC->OE.ecc_ave_temp = SC->OE.eccentricity;
    SC->OE.ave_increm = 1;
    SC->et_last_orbit = SC->et;
    SC->orbit_number = SC->orbit_number + 1;
  }
  else{
    //    SC->OE.w_ave = //9999.999999/RAD2DEG;
    SC->OE.w_ave_temp = SC->OE.w_ave_temp + SC->OE.w;

    //    SC->OE.sma_ave = 9999.999999;
    SC->OE.sma_ave_temp = SC->OE.sma_ave_temp + SC->OE.sma;
    SC->OE.ecc_ave_temp = SC->OE.ecc_ave_temp + SC->OE.eccentricity;


    SC->OE.ave_increm = SC->OE.ave_increm + 1;
  }
    // end of to cmopute orbita average elements


  return 0;
}

int eci2lla(double pos[3], double et,  double geodetic[3]  ){ // convert ECI xyz to lon/lat/alt. x,y,z in km. time in seconds past J2000.

  // convert time seconds past J2000 to julian date in ut1 days from 4713 bc   
  // // convert time seconds past J2000 to time string 
  char cu_time[290];
  et2utc_c(et, "ISOC" ,3 ,255 , cu_time); 
  // // convert time string to year month day hour minute seconds (all int, except seconds which is double)
  int year_current, month_current, day_current, hour_current, minute_current;
  double second_current;
  char year_current_str[10];	
  strcpy(year_current_str, ""); strncat(year_current_str, &cu_time[0], 4); year_current = atoi(year_current_str);
  char month_current_str[10];	
  strcpy(month_current_str, ""); strncat(month_current_str, &cu_time[5], 2); month_current = atoi(month_current_str);
  char day_current_str[10];	
  strcpy(day_current_str, ""); strncat(day_current_str, &cu_time[8], 2); day_current = atoi(day_current_str);
  char hour_current_str[10];	
  strcpy(hour_current_str, ""); strncat(hour_current_str, &cu_time[11], 2); hour_current = atoi(hour_current_str);
  char minute_current_str[10];	
  strcpy(minute_current_str, ""); strncat(minute_current_str, &cu_time[14], 2); minute_current = atoi(minute_current_str);
  char second_current_str[10];	
  int second_current_no_microsec;
  strcpy(second_current_str, ""); strncat(second_current_str, &cu_time[17], 2); second_current_no_microsec = atoi(second_current_str);
  char microsecond_current_str[10];	
  int microsecond;
  strcpy(microsecond_current_str, ""); strncat(microsecond_current_str, &cu_time[20], 3); microsecond = atoi(microsecond_current_str);
  second_current = second_current_no_microsec + microsecond/1000.;

  //  printf("%d %d %d %d %d %f\n", year_current, month_current, day_current, hour_current, minute_current, second_current);
  double time ;//= 2453101.82741;
  // // Convert year month day hhour min sec to Julian day (in ut1 days from 4713 bc)
  jday (year_current, month_current, day_current, hour_current, minute_current, second_current, &time );



  double f = 1.0/298.26;
  double xkmper = 6378.135;
  double theta = atan2(pos[1],pos[0]);
  double lon = fmod( theta - gstime( time ), 2*M_PI ) ; 
    double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double e2 = f*(2-f);
    double lat = atan2(pos[2],r);
    
    double phi = lat;
    double c = 1.0/sqrt( 1 - e2*sin(phi)*sin(phi) );
    lat = atan2(pos[2]+xkmper*c*e2*sin(phi),r);
    while (fabs(lat - phi) > 0.0000000001){ 
      phi = lat;
      c = 1.0/sqrt( 1 - e2*sin(phi)*sin(phi) );
      lat = atan2(pos[2]+xkmper*c*e2*sin(phi),r);
      //      printf("%e %d\n", fabs(lat - phi), lat == phi);
    }

    double   alt = r/cos(lat) - xkmper*c;
    geodetic[0] = lat;
      geodetic[1] = lon;
      geodetic[2] = alt;
      if (lon < -M_PI){
	printf("lon = %f\n", lon*180./M_PI);
      }
      

  return 0;
}

/* -----------------------------------------------------------------------------
*
*                           function gstime
*
*  this function finds the greenwich sidereal time.
*
*  author        : Charles Bussy-Virat, largely inspired from SGP4 by david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    jdut1       - julian date in ut1             days from 4713 bc
*
*  outputs       :
*    gstime      - greenwich sidereal time        0 to 2pi rad
*
*  locals        :
*    temp        - temporary variable for doubles   rad
*    tut1        - julian centuries from the
*                  jan 1, 2000 12 h epoch (ut1)
*
*  coupling      :
*    none
*
*  references    :
*    vallado       2004, 191, eq 3-45
* --------------------------------------------------------------------------- */

double  gstime (double jdut1 )
   {
     const double twopi = 2.0 * M_PI;
     const double deg2rad = M_PI / 180.0;
     double       temp, tut1;

     tut1 = (jdut1 - 2451545.0) / 36525.0;

     temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
             (876600.0*3600 + 8640184.812866) * tut1 + 67310.54841;  // sec
     temp = fmod(temp * deg2rad / 240.0, twopi); //360/86400 = 1/240, to deg, to rad

     // ------------------------ check quadrants ---------------------
     if (temp < 0.0)
         temp += twopi;

     return temp;
   }  // end gstime

/* -----------------------------------------------------------------------------
*
*                           procedure jday
*
*  this procedure finds the julian date given the year, month, day, and time.
*    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
*
*  algorithm     : calculate the answer in one step for efficiency
*
*  author        : Charles Bussy-Virat, largely inspired from SGP4 by david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    year        - year                           1900 .. 2100
*    mon         - month                          1 .. 12
*    day         - day                            1 .. 28,29,30,31
*    hr          - universal time hour            0 .. 23
*    min         - universal time min             0 .. 59
*    sec         - universal time sec             0.0 .. 59.999
*
*  outputs       :
*    jd          - julian date                    days from 4713 bc
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
*
*  references    :
*    vallado       2007, 189, alg 14, ex 3-14
*
* --------------------------------------------------------------------------- */

int jday (int year, int mon, int day, int hr, int minute, double sec, double *jd )
   {
     *jd = 367.0 * year -
          floor((7 * (year + floor((mon + 9) / 12.0))) * 0.25) +
          floor( 275 * mon / 9.0 ) +
          day + 1721013.5 +
          ((sec / 60.0 + minute) / 60.0 + hr) / 24.0;  // ut in days
          // - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
     return 0;
   }  // end jday





/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           coverage_ground_station
//  Purpose:        Computes the coverage of a satellite above all ground stations and calculate a few parameters: elevation, azimuth, range in different two coordinate systems (ECEF, ENU) and in two different frames of reference (satellite and ground station)
//  Assumptions:    if the sc is a GPS that was initialized in the second line of section #SPACECRAFT of the main input file, its attitude is assumed to be nadir pointing
//  References:     /
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 06/01/2016    |   ---     | Initial Implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int coverage_ground_station(   SPACECRAFT_T *SC, // in: spacecraft (position and time). out: elevation, azimuth, range 
			       GROUND_STATION_T *GROUND_STATION, // in: list of ground stations
			       PARAMS_T *PARAMS, // in: parameters for the propagation (Earth flattening and radius here)
			       int index_in_attitude_interpolated, // in: index of the current time step (take into account RK4)
			       INTEGRATOR_T    *INTEGRATOR, // in: paramters for the propagation (attitude of the sc here)
			       double          et_initial_epoch, // in: time of the inital epoch of the constellation
			       double et_sc_initial,
			       double sc_ecef_previous_time_step[3], // in: ECEF position of the spacecraft at the previous time step of the propagation (used to linear interpolate the position between the previous position and the current position)
			       double sc_eci_previous_time_step[3], //  in: ECI position of the spacecraft at the previous time step of the propagation (used to linear interpolate the position between the previous position and the current position)
			       double sc_v_eci_previous_time_step[3], //  in: ECI velocity of the spacecraft at the previous time step of the propagation (used to linear interpolate the position between the previous position and the current position)
			       double time_step_interpolation // in: time step of the linear interpolation of the sc position (in seconds)
			       )

{

  /* Declaration */


  //      char time_et[256];
  double norm_unit_ground_station_to_sc_in_enu_projected_on_east_north_local_plane[3];
  double norm_unit_sc_to_ground_station_body_projected_on_body_xy[3];
  double unit_sc_to_ground_station_body[3];
  double sc_to_ground_station_dot_minus_z_axis_body;	
  double unit_sc_to_ground_station_body_projected_on_body_xy_dot_y_axis_body;
  double x_axis_body_dot_sc_to_ground_station_body_projected_on_body_xy;
  double unit_sc_to_ground_station_body_projected_on_body_xy[3];
  double minus_z_axis_body[3];
  int step_end_interpolation;
  int iground;
  double ground_station_to_sc_in_ecef[3];
  double opposite_ground_station[3];
  double T_ECEF_to_J2000[3][3], T_inrtl_2_lvlh[3][3], T_sc_to_lvlh[3][3], T_lvlh_to_sc[3][3];
  double ground_station_to_sc_eci[3], ground_station_to_sc_lvlh[3],ground_station_to_sc_body[3];
  double v_angle[3];
  int order_rotation[3];
  double x_axis_body[3], unit_ground_station_to_sc_body[3];
  double ground_station_to_sc_in_enu[3];
  double ground_station_to_sc_dot_local_vertical_in_enu;
  double T_enu_to_ecef[3][3], T_ecef_to_enu[3][3];
  double local_vertical_in_enu[3];
  double local_north_in_enu[3];
  double ground_station_to_sc_in_enu_projected_on_east_north_local_plane[3];
  double ground_station_to_sc_in_enu_projected_on_east_north_local_plane_dot_local_north_in_enu;
  char time_interpolated[256];
  double et_interpolated, et_previous_time_step;
  double ecef_sc[3], eci_sc[3], eci_v_sc[3];
  int in_sight_of_ground_station;
  double elevation_wtr_to_ground_station_in_sc_refce, azimuth_wtr_to_ground_station_in_sc_refce, range_wtr_to_ground_station, elevation_wtr_to_ground_station_in_ground_station_refce, azimuth_wtr_to_ground_station_in_ground_station_refce;
  int istep;
  double unit_ground_station_to_sc_in_enu[3];
  double y_axis_body[3];
  double unit_ground_station_to_sc_in_enu_projected_on_east_north_local_plane[3];
  double ground_station_to_sc_in_enu_projected_on_east_north_local_plane_dot_local_east_in_enu;
  double local_east_in_enu[3];

  et_previous_time_step = SC->et - INTEGRATOR->dt;

  /* Algorithm */
  for (iground = 0; iground < GROUND_STATION->nb_ground_stations; iground++){
    // Write the header of the coverage file
    if ( (index_in_attitude_interpolated == INTEGRATOR->index_in_attitude_interpolated_first + 2) || (fabs(et_initial_epoch - SC->et) <= INTEGRATOR->dt/2.)  ){ // Write the header at the first time step of the propagation (the second condition is in case the initialization was made with TLEs, the fisrt time to write is when the SC->et is eqaual to initial_epohc (we dont want to write at the first time step right after the TLE initial epoch_
      fprintf(SC->fp_coverage_ground_station[iground], "// ---------------------------------------------------------------------------------- \n");
      fprintf(SC->fp_coverage_ground_station[iground], "// \n");
      fprintf(SC->fp_coverage_ground_station[iground], "// Please note this file is auto generated.\n");
      fprintf(SC->fp_coverage_ground_station[iground], "// Coverage of Spacecraft %s for the ground station %s.\n", SC->name_sat,GROUND_STATION->name_ground_station[iground]);
      fprintf(SC->fp_coverage_ground_station[iground], "// \n");
      fprintf(SC->fp_coverage_ground_station[iground], "// Version control is under Joel Getchius' Mac and Charles Bussy-Virat University of Michigan's Mac \n");
      fprintf(SC->fp_coverage_ground_station[iground], "// \n");
      fprintf(SC->fp_coverage_ground_station[iground], "// ---------------------------------------------------------------------------------- \n");
      fprintf(SC->fp_coverage_ground_station[iground], "// Ground station %s ECEF: %f %f %f\n", GROUND_STATION->name_ground_station[iground],GROUND_STATION->ecef_ground_station[iground][0], GROUND_STATION->ecef_ground_station[iground][1], GROUND_STATION->ecef_ground_station[iground][2]);
      fprintf(SC->fp_coverage_ground_station[iground], "// Ground station %s longitude(degree) latitude(degree) altitude(m): %f %f %f\n", GROUND_STATION->name_ground_station[iground],GROUND_STATION->longitude_ground_station[iground]*RAD2DEG, GROUND_STATION->latitude_ground_station[iground]*RAD2DEG, GROUND_STATION->altitude_ground_station[iground]);
      fprintf(SC->fp_coverage_ground_station[iground], "TIME ECEF_X_SC ECEF_Y_SC ECEF_Z_SC ELEV_REFCE_SC AZI_REFCE_SC ELEV_REFCE_GROUND_STATION AZI_REFCE_GROUND_STATION RANGE_SC_TO_GROUND_STATION\n");
      fprintf(SC->fp_coverage_ground_station[iground], "#START \n");
    }

    if ( ( index_in_attitude_interpolated == INTEGRATOR->index_in_attitude_interpolated_first + 2 ) ) {
      step_end_interpolation = (int)(nearbyint(INTEGRATOR->dt / time_step_interpolation + 1));
      //      printf("%s: %d | %d - %d\n", SC->name_sat,step_end_interpolation, GROUND_STATION->nb_ground_stations, iground);
      et_interpolated = et_previous_time_step - 1;
    }
    else{
      step_end_interpolation = (int)(nearbyint(INTEGRATOR->dt / time_step_interpolation));
      //      printf("%s: %d | %d - %d\n", SC->name_sat,step_end_interpolation, GROUND_STATION->nb_ground_stations, iground);
      et_interpolated = et_previous_time_step;
    }
    for (istep = 0; istep < step_end_interpolation; istep++){

      // Linear interpolation of the sc position between the previous time step and the current time step
      et_interpolated = et_interpolated + time_step_interpolation;
      /* printf("%d\n", istep); */
      /* etprint(et_interpolated, ""); */
      
      // // ECEF position interpolation
      ecef_sc[0] = sc_ecef_previous_time_step[0] + ( ( et_interpolated - et_previous_time_step ) / INTEGRATOR->dt  ) * ( SC->r_ecef2cg_ECEF[0] - sc_ecef_previous_time_step[0] ) ;
      ecef_sc[1] = sc_ecef_previous_time_step[1] + ( ( et_interpolated - et_previous_time_step ) / INTEGRATOR->dt ) * ( SC->r_ecef2cg_ECEF[1] - sc_ecef_previous_time_step[1] ) ;
      ecef_sc[2] = sc_ecef_previous_time_step[2] + ( ( et_interpolated - et_previous_time_step ) / INTEGRATOR->dt ) * ( SC->r_ecef2cg_ECEF[2] - sc_ecef_previous_time_step[2] ) ;
      // // ECI position interpolation
      eci_sc[0] = sc_eci_previous_time_step[0] + ( ( et_interpolated - et_previous_time_step ) / INTEGRATOR->dt  ) * ( SC->r_i2cg_INRTL[0] - sc_eci_previous_time_step[0] ) ;
      eci_sc[1] = sc_eci_previous_time_step[1] + ( ( et_interpolated - et_previous_time_step ) / INTEGRATOR->dt ) * ( SC->r_i2cg_INRTL[1] - sc_eci_previous_time_step[1] ) ;
      eci_sc[2] = sc_eci_previous_time_step[2] + ( ( et_interpolated - et_previous_time_step ) / INTEGRATOR->dt ) * ( SC->r_i2cg_INRTL[2] - sc_eci_previous_time_step[2] ) ;
      // // ECI velocity interpolation
      eci_v_sc[0] = sc_v_eci_previous_time_step[0] + ( ( et_interpolated - et_previous_time_step ) / INTEGRATOR->dt  ) * ( SC->v_i2cg_INRTL[0] - sc_v_eci_previous_time_step[0] ) ;
      eci_v_sc[1] = sc_v_eci_previous_time_step[1] + ( ( et_interpolated - et_previous_time_step ) / INTEGRATOR->dt ) * ( SC->v_i2cg_INRTL[1] - sc_v_eci_previous_time_step[1] ) ;
      eci_v_sc[2] = sc_v_eci_previous_time_step[2] + ( ( et_interpolated - et_previous_time_step ) / INTEGRATOR->dt ) * ( SC->v_i2cg_INRTL[2] - sc_v_eci_previous_time_step[2] ) ;

      // Figures out if the satellite is above the local horizon of the ground station
      // // Local vertical in the ENU frame attached to the ground station  (ENU: East-North-Up)
      local_vertical_in_enu[0] = 0; local_vertical_in_enu[1] = 0; local_vertical_in_enu[2] = 1;
      // // vector ground station to sc
      v_scale(opposite_ground_station, GROUND_STATION->ecef_ground_station[iground], -1);
      v_add(ground_station_to_sc_in_ecef, opposite_ground_station, ecef_sc);
      // // convert vector ground station to sc from ECEF to ENU coordinates of the ground station
      compute_T_enu_to_ecef( T_enu_to_ecef, GROUND_STATION->latitude_ground_station[iground], GROUND_STATION->longitude_ground_station[iground], PARAMS->EARTH.flattening);
      m_trans(T_ecef_to_enu, T_enu_to_ecef);
      m_x_v(ground_station_to_sc_in_enu, T_ecef_to_enu, ground_station_to_sc_in_ecef);
      // // Dot product between the local vertical and the vector ground station to sc in ENU reference system
      v_norm(unit_ground_station_to_sc_in_enu, ground_station_to_sc_in_enu);
      v_dot(&ground_station_to_sc_dot_local_vertical_in_enu, unit_ground_station_to_sc_in_enu, local_vertical_in_enu);
      // Elevation angle of the sc with respect to the ground station in the reference frame of the ground station
      elevation_wtr_to_ground_station_in_ground_station_refce = M_PI/2. - acos( ground_station_to_sc_dot_local_vertical_in_enu );

      if ( elevation_wtr_to_ground_station_in_ground_station_refce >= GROUND_STATION->min_elevation_angle_ground_station[iground] ){ // if the sc is above the local horizon + min elevation angle of the ground station
    	in_sight_of_ground_station = 1;

    	/************************************************************************************************************/
    	/* COMPUTE ELVATION AND AZIMUTH ANGLES IN THE REFERENCE FRAME OF THE SC */
    	// Elevation angle of the sc with respect to the ground station in the reference frame of the sc: angle between the sc xy body plane and the vector sc to ground station (expressed in the sc body reference system). Counted negative if the ground station is bwlo the body xy plane, which is always the case in a nadir configuration and when the sc is in sight.
    	// // convert the ground station to sc vector in the sc body reference system
    	// // // ECEF to ECI
    	pxform_c(PARAMS->EARTH.earth_fixed_frame, "J2000", et_interpolated, T_ECEF_to_J2000);
    	m_x_v(ground_station_to_sc_eci, T_ECEF_to_J2000, ground_station_to_sc_in_ecef);
    	// // // ECI to LVLH
    	compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, eci_sc, eci_v_sc);
    	m_x_v(ground_station_to_sc_lvlh, T_inrtl_2_lvlh, ground_station_to_sc_eci);
    	// // // LVLH to BODY
    	// // // // Find the current attitude of the sc
	if (INTEGRATOR->isGPS == 0){ // if the sc is not a GPS that was initialized in the second line of section #SPACECRAFT of the main input file
	  if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") == 0){
	    v_angle[0] = INTEGRATOR->attitude.pitch_angular_velocity_ensemble * ( et_interpolated - et_sc_initial ) + INTEGRATOR->attitude.pitch_ini_ensemble;
	    v_angle[1] = INTEGRATOR->attitude.roll_angular_velocity_ensemble * ( et_interpolated - et_sc_initial ) + INTEGRATOR->attitude.roll_ini_ensemble;
	    v_angle[2] = INTEGRATOR->attitude.yaw_angular_velocity_ensemble * ( et_interpolated - et_sc_initial ) + INTEGRATOR->attitude.yaw_ini_ensemble;
	    order_rotation[0] = 1; // !!!!!!!! we might want to change that in the future
	    order_rotation[1] = 2;
	    order_rotation[2] = 3;
	    INTEGRATOR->attitude.pitch_current = v_angle[0];
	    INTEGRATOR->attitude.roll_current = v_angle[1];
	    INTEGRATOR->attitude.yaw_current = v_angle[2];
	        INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


	  }

	  if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") == 0) {
	      v_angle[0] =  INTEGRATOR->attitude.pitch_for_attitude_ensemble  +  INTEGRATOR->attitude.pitch_angular_velocity_constant * ( et_interpolated - et_sc_initial );
	      v_angle[1] =  INTEGRATOR->attitude.roll_for_attitude_ensemble  +  INTEGRATOR->attitude.roll_angular_velocity_constant * ( et_interpolated - et_sc_initial );
	      v_angle[2] =  INTEGRATOR->attitude.yaw_for_attitude_ensemble  +  INTEGRATOR->attitude.yaw_angular_velocity_constant * ( et_interpolated - et_sc_initial );

	      order_rotation[0]  = 1; order_rotation[1]  = 2; order_rotation[2]  = 3;
	   
	      INTEGRATOR->attitude.pitch_current = v_angle[0];
	      INTEGRATOR->attitude.roll_current = v_angle[1];
	      INTEGRATOR->attitude.yaw_current = v_angle[2];
	          INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


	  }
   
	  if ( (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "sun_pointed") != 0) ){ // otherwise (atittude is nadir, sun_pointed or manual (from an input file))  
	    if (INTEGRATOR->file_is_quaternion == 0){
	    v_angle[0] = INTEGRATOR->attitude.pitch[index_in_attitude_interpolated];
	    v_angle[1] = INTEGRATOR->attitude.roll[index_in_attitude_interpolated];
	    v_angle[2] = INTEGRATOR->attitude.yaw[index_in_attitude_interpolated];
	    order_rotation[0] = INTEGRATOR->attitude.order_pitch[index_in_attitude_interpolated];
	    order_rotation[1] = INTEGRATOR->attitude.order_roll[index_in_attitude_interpolated];
	    order_rotation[2] = INTEGRATOR->attitude.order_yaw[index_in_attitude_interpolated];
	    INTEGRATOR->attitude.pitch_current = v_angle[0];
	    INTEGRATOR->attitude.roll_current = v_angle[1];
	    INTEGRATOR->attitude.yaw_current = v_angle[2];
	        INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


	    }
	    else{
	      q_copy( INTEGRATOR->attitude.quaternion_current, INTEGRATOR->attitude.quaternion[index_in_attitude_interpolated]);
	    }
	    

	  }
	}
	else{ // if the sc is a GPS that was initialized in the second line of section #SPACECRAFT of the main input file: its attitude is assumed to be nadir pointing
	  v_angle[0] = 0;
	  v_angle[1] = 0;
	  v_angle[2] = 0;
	  order_rotation[0] = 1;
	  order_rotation[1] = 2;
	  order_rotation[2] = 3;
	    INTEGRATOR->attitude.pitch_current = v_angle[0];
	    INTEGRATOR->attitude.roll_current = v_angle[1];
	    INTEGRATOR->attitude.yaw_current = v_angle[2];
    INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


	}

    	// // // // Now we know the current attitude of the sc, convert the ground station to sc vector from LVLH to SC body reference frame
	//

	compute_T_sc_to_lvlh( T_sc_to_lvlh, v_angle, order_rotation, INTEGRATOR->attitude.attitude_profile, &et_interpolated,  eci_sc, eci_v_sc, INTEGRATOR->file_is_quaternion, INTEGRATOR->attitude.quaternion_current, PARAMS);
    	m_trans(T_lvlh_to_sc, T_sc_to_lvlh);


    	m_x_v(ground_station_to_sc_body, T_lvlh_to_sc, ground_station_to_sc_lvlh);
	// // Calculate the elevation angle
    	v_norm( unit_ground_station_to_sc_body, ground_station_to_sc_body);
	v_scale(unit_sc_to_ground_station_body, unit_ground_station_to_sc_body, -1);
    	minus_z_axis_body[0] = 0, minus_z_axis_body[1] = 0, minus_z_axis_body[2] = -1; // coordinates of the -z_axis body in the sc body reference frame
    	v_dot( &sc_to_ground_station_dot_minus_z_axis_body, unit_sc_to_ground_station_body, minus_z_axis_body );
    	elevation_wtr_to_ground_station_in_sc_refce = -(M_PI/2- acos( sc_to_ground_station_dot_minus_z_axis_body )); // the elevation is counted from the sc xy body, negatively if the ground station is below the xy body plane
    	// Azimuth angle of the sc with respect to the ground station in the reference frame of the sc
    	// This is the angle between the x axis and the sc to ground station vector (expressed in the sc body reference system) projected on the body xy plane
    	// // Computes the angle, in the sc body system, between the body x axis and the sc to ground station vector projected on the body xy plane
    	x_axis_body[0] = 1, x_axis_body[1] = 0, x_axis_body[2] = 0; // coordinates of the x_axis body in the sc body reference frame
    	unit_sc_to_ground_station_body_projected_on_body_xy[0] = unit_sc_to_ground_station_body[0]; // project the unit_sc_to_ground_station_body on the body xy plane
    	unit_sc_to_ground_station_body_projected_on_body_xy[1] = unit_sc_to_ground_station_body[1]; // project the unit_sc_to_ground_station_body on the body xy plane
    	unit_sc_to_ground_station_body_projected_on_body_xy[2] = 0; // project the unit_sc_to_ground_station_body on the body xy plane
	v_norm(norm_unit_sc_to_ground_station_body_projected_on_body_xy, unit_sc_to_ground_station_body_projected_on_body_xy);
    	v_dot( &x_axis_body_dot_sc_to_ground_station_body_projected_on_body_xy, x_axis_body, norm_unit_sc_to_ground_station_body_projected_on_body_xy);
    	azimuth_wtr_to_ground_station_in_sc_refce = acos( x_axis_body_dot_sc_to_ground_station_body_projected_on_body_xy ) ;
    	// // Azimuth angle goes from 0 to 360 in the body -y axis direction (so to the right of the +x axis vector if sc seen from above) (ex: if azimuth = 90 then the ground station is in body -y axis direction of the sc, so starboard). This is to be conistent with the azimuth angle in the ground station reference frame where it is calculated in the +East direction (so the right of the North direction)
    	y_axis_body[0] = 0; y_axis_body[1] = 1; y_axis_body[2] = 0;
    	v_dot(&unit_sc_to_ground_station_body_projected_on_body_xy_dot_y_axis_body, y_axis_body, norm_unit_sc_to_ground_station_body_projected_on_body_xy); 
    	if ( unit_sc_to_ground_station_body_projected_on_body_xy_dot_y_axis_body > 0){ // if the ground station is in the positive body y plane (so port of the sc if sc nadir)
    	  azimuth_wtr_to_ground_station_in_sc_refce = 2*M_PI - azimuth_wtr_to_ground_station_in_sc_refce;
    	}
    	/************************************************************************************************************/
    	/* COMPUTE AZIMUTH ANGLES IN THE REFERENCE FRAME OF THE GROUND STATION */
    	// Azimuth angle of the sc with respect to the ground station in the reference frame of the ground station
    	// This is the angle between the local North of the ground station and the ground station to sc vector (expressed in the ENU reference system of the ground station)  projected on the East-North local plane (which is the xy plane of the ENU)
    	// // Local North in ENU reference frame
    	local_north_in_enu[0] = 0; local_north_in_enu[1] = 1; local_north_in_enu[2] = 0;
    	// // Projection of the ground station to sc vector (expressed in the ENU reference system of the ground station) on the East-North local plane
    	ground_station_to_sc_in_enu_projected_on_east_north_local_plane[0] = ground_station_to_sc_in_enu[0];
    	ground_station_to_sc_in_enu_projected_on_east_north_local_plane[1] = ground_station_to_sc_in_enu[1];
    	ground_station_to_sc_in_enu_projected_on_east_north_local_plane[2] = 0;
    	// // Dot product between the local North and the ground station to sc vector (expressed in the ENU reference system of the ground station) projected on the East-North local plane
	v_norm(unit_ground_station_to_sc_in_enu_projected_on_east_north_local_plane, ground_station_to_sc_in_enu_projected_on_east_north_local_plane);
	v_norm(norm_unit_ground_station_to_sc_in_enu_projected_on_east_north_local_plane, unit_ground_station_to_sc_in_enu_projected_on_east_north_local_plane);
    	v_dot(&ground_station_to_sc_in_enu_projected_on_east_north_local_plane_dot_local_north_in_enu, norm_unit_ground_station_to_sc_in_enu_projected_on_east_north_local_plane, local_north_in_enu);
    	azimuth_wtr_to_ground_station_in_ground_station_refce = acos( ground_station_to_sc_in_enu_projected_on_east_north_local_plane_dot_local_north_in_enu );
    	// // Azimuth angle goes from 0 to 360 in the East direction (ex: if azimuth = 90 then the sc is in the direction of the East wtr to the ground station)
    	local_east_in_enu[0] = 1; local_east_in_enu[1] = 0; local_east_in_enu[2] = 0;
    	v_dot(&ground_station_to_sc_in_enu_projected_on_east_north_local_plane_dot_local_east_in_enu, norm_unit_ground_station_to_sc_in_enu_projected_on_east_north_local_plane, local_east_in_enu);
    	if ( ground_station_to_sc_in_enu_projected_on_east_north_local_plane_dot_local_east_in_enu < 0){ // if the sc is in the West direction
    	  azimuth_wtr_to_ground_station_in_ground_station_refce = 2*M_PI - azimuth_wtr_to_ground_station_in_ground_station_refce;
    	}

    	/************************************************************************************************************/
    	// Range of the sc with respect to the ground station
    	v_mag(&range_wtr_to_ground_station,ground_station_to_sc_in_ecef);
      // Write the results in a file: Time ECEF_sc elevation_refce_sc azimuth_refce_sc elevation_refce_ground_station azimuth_refce_ground_station (in the header I indicate ECEF_ground_station and lla)
      et2utc_c(et_interpolated, "ISOC" ,0 ,255 , time_interpolated);
      //     printf("e_sc %f | a_sc %f | e_g %f | a_g %f\n",elevation_wtr_to_ground_station_in_sc_refce*RAD2DEG, azimuth_wtr_to_ground_station_in_sc_refce*RAD2DEG, elevation_wtr_to_ground_station_in_ground_station_refce*RAD2DEG, azimuth_wtr_to_ground_station_in_ground_station_refce*RAD2DEG);
      if (istep == step_end_interpolation-1){ // only print every dt_output
	//	      etprint(et_interpolated, "");
	//fprintf(SC->fp_coverage_ground_station[iground], "%s %f %f %f %f %f %f %f %f\n",time_interpolated, ecef_sc[0], ecef_sc[1], ecef_sc[2], elevation_wtr_to_ground_station_in_sc_refce*RAD2DEG, azimuth_wtr_to_ground_station_in_sc_refce*RAD2DEG, elevation_wtr_to_ground_station_in_ground_station_refce*RAD2DEG, azimuth_wtr_to_ground_station_in_ground_station_refce*RAD2DEG, range_wtr_to_ground_station); // !!!!!!!!!!!!! uncomment this line and comment line below
	//	fprintf(SC->fp_coverage_ground_station[iground], "%s %f\n",time_interpolated, elevation_wtr_to_ground_station_in_ground_station_refce*RAD2DEG); // !!!!!!!!!!!!! comment this line and uncommnet line above
		fprintf(SC->fp_coverage_ground_station[iground], "%s %f\n",time_interpolated,90 - (-elevation_wtr_to_ground_station_in_sc_refce*RAD2DEG)); // !!!!!!!! this is to see if the gs is in the field of view of the sc. 0 means the gs is nadir of the sc. CRAP


      }

      }
      else{
    	in_sight_of_ground_station = 0;
      et2utc_c(et_interpolated, "ISOC" ,0 ,255 , time_interpolated);
	elevation_wtr_to_ground_station_in_sc_refce = -999; azimuth_wtr_to_ground_station_in_sc_refce = -999; elevation_wtr_to_ground_station_in_ground_station_refce = -999; azimuth_wtr_to_ground_station_in_ground_station_refce = -999; range_wtr_to_ground_station = -999;
	//	fprintf(SC->fp_coverage_ground_station[iground], "%s %f %f %f %f %f %f %f %f\n",time_interpolated, ecef_sc[0], ecef_sc[1], ecef_sc[2], elevation_wtr_to_ground_station_in_sc_refce, azimuth_wtr_to_ground_station_in_sc_refce, elevation_wtr_to_ground_station_in_ground_station_refce, azimuth_wtr_to_ground_station_in_ground_station_refce, range_wtr_to_ground_station);
      }


    }
  }
  /* if ( index_in_attitude_interpolated == 16 ){ */
  /*   exit(0); */
  /* } */

  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           geodetic_to_geocentric
//  Purpose:        Computes planet fixed state based on lat/long/altitude FOR LOCATIONS ON THE EARTH'S SURFACE ONLY (NOT FOR SATELLITES) (although I don't know, seems it's working ok...)
//  Assumptions:    Works only FOR LOCATIONS ON THE EARTH'S SURFACE ONLY (NOT FOR SATELLITES) (although I don't know, seems it's working ok...). This assumes the Earth is an ellipsoid. To use a spherical Eart, set the flattening parameter to 0 
//  References      Valado 3rd edition section 3.2
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 10/01/2015    |   ---     | Initial Implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int geodetic_to_geocentric(double flattening,            /* IN:     flattening parameter     */
			   double h,                     /* IN:  M  height above ellipsoid  (km) */
			   double lat,                   /* IN:  r  geodetic latitude (rad)      */
			   double longitude,              /* IN:  r  longitude  (rad)             */
			   double equatorial_radius,       /* IN: equatorial radius  (km) */
			   double  R_ecef_2_cg_ECEF[3])   /* OUT:     vector in ECEF           */

{
  /* Declarations */
  double c_earth, s_earth, earth_eccentricity, deno, r_xy_mag;

  /* Algorithm */
  earth_eccentricity = sqrt( 1 - (1 - flattening)*(1-flattening) );
  deno = sqrt( 1-earth_eccentricity*earth_eccentricity*sin(lat)*sin(lat) );
  c_earth = equatorial_radius / deno;
  s_earth =  equatorial_radius * (1 - earth_eccentricity*earth_eccentricity) / deno;

  r_xy_mag = ( c_earth + h ) * cos(lat);
  R_ecef_2_cg_ECEF[0] = r_xy_mag * cos(longitude);
  R_ecef_2_cg_ECEF[1] = r_xy_mag * sin(longitude);

  R_ecef_2_cg_ECEF[2] = ( s_earth + h ) * sin(lat);

  // OHTER SOLUTION (SOURCE: http://www.mathworks.com/help/aeroblks/llatoecefposition.html)
  /* /\* Declaration *\/ */
  /* double geocentric_lat_mean_sea_level; */
  /* double radius_at_surface_point; */
  /* double term1, term2; */

  /* /\* Algorithm *\/ */
  /* geocentric_lat_mean_sea_level = atan( ( 1 - flattening ) * ( 1 - flattening ) * tan( lat ) ); */
  
  /* term1 = ( 1 - flattening ) * ( 1 - flattening ); */
  /* term2 = ( 1 / term1 ) - 1 ; */
  /* radius_at_surface_point = sqrt( equatorial_radius * equatorial_radius / ( 1 + term2 *  sin(geocentric_lat_mean_sea_level)  *  sin(geocentric_lat_mean_sea_level ) ) ); */

  /* R_ecef_2_cg_ECEF[0] = radius_at_surface_point * cos(geocentric_lat_mean_sea_level) * cos(longitude) + h * cos(lat) * cos(longitude); */
  /* R_ecef_2_cg_ECEF[1] = radius_at_surface_point * cos(geocentric_lat_mean_sea_level) * sin(longitude) + h * cos(lat) * sin(longitude); */
  /* R_ecef_2_cg_ECEF[2] = radius_at_surface_point * sin(geocentric_lat_mean_sea_level) + h * sin(lat); */

  /* v_print(R_ecef_2_cg_ECEF, "r_ecef"); */
  /* exit(0); */
  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           Geocentric_to_geodetic
//  Purpose:        Computes lat/long/altitude based on planet fixed state
//  Assumptions:    None.
//  References      Larson & Wertz
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/09/2015    |   ---     | Initial Implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int geocentric_to_geodetic(
			   double  R_ecef_2_cg_ECEF[3],   /* IN:     vector in ECEF           */
			   double *semimajor_axis,        /* IN:     planetary semimajor axis */
			   double *flattening,            /* IN:     flattening parameter     */
			   double *h,                     /* OUT:  M  height above ellipsoid  */
			   double *lat,                   /* OUT:  r  geodetic latitude       */
			   double *longitude   )          /* OUT:  r  longitude               */
{

  /* Variable declarations */
  double b;
  double R;
  double E;
  double F;
  double P;
  double Q;
  double D;
  double v;
  double G;
  double t;
  double inv3;
  double D2;
  double Domain;
  double D_minus_Q;
  double D_plus_Q;
  double D_minus_Q_pow;
  double D_plus_Q_pow;


  /* Algorithm */
  inv3 = 1.0 / 3.0;

  b = -((*flattening * *semimajor_axis) - *semimajor_axis);

  if (R_ecef_2_cg_ECEF[2] < 0.0) {
    b = (-b);
  }

  R = sqrt( R_ecef_2_cg_ECEF[0]* R_ecef_2_cg_ECEF[0] +  R_ecef_2_cg_ECEF[1]* R_ecef_2_cg_ECEF[1] );

  if ((*semimajor_axis * R ) > SMALL_NUM) {
    E = ( (b*R_ecef_2_cg_ECEF[2] - (*semimajor_axis * *semimajor_axis - b*b))) / (*semimajor_axis *R);
    F = ( (b*R_ecef_2_cg_ECEF[2] + (*semimajor_axis * *semimajor_axis - b*b))) / (*semimajor_axis *R);
  } else {
    E = 0.0;
    F = 0.0;
  }

  P = 4*(E*F + 1.0 )/ 3.0;
  Q = 2*(E*E - F*F);

  /* Sqrt protection */
  D2 = P*P*P + Q*Q;

  if (D2 > 0.0 ) {
    D = sqrt(D2);
  } else {
    D = 0.0;
  }

  /* Pow protection */
  D_minus_Q = D - Q;
  D_plus_Q = D + Q;

  if (D_minus_Q < 0.0) {
    D_minus_Q_pow = pow(fabs(D_minus_Q ), inv3 );
    D_minus_Q_pow = (-D_minus_Q_pow);
  } else {
    D_minus_Q_pow = pow( D_minus_Q , inv3);
  }

  if (D_plus_Q < 0.0) {
    D_plus_Q_pow = pow(fabs(D_plus_Q ), inv3 );
    D_plus_Q_pow = (-D_plus_Q_pow);
  } else {
    D_plus_Q_pow = pow( D_plus_Q , inv3);
  }

  v = D_minus_Q_pow - D_plus_Q_pow;

  Domain = E*E + v;

  /* Sqrt protection */
  if ( Domain > 0.0 ) {
    G = 0.5*(sqrt( Domain ) + E);
  }  else {
    G = 0.0;
  }

  if (fabs(2*G - E) > SMALL_NUM) {
    t = sqrt(G*G + (F - v*G) / (2*G - E)) - G;
  } else {
    t = 0.0 ;
  }


  if (fabs(2.0*b*t) > SMALL_NUM ) {
    *lat = atan( *semimajor_axis * (1.0 - t*t) / (2.0*b*t));
  } else {
    *lat = 0.0;
  }

  if (fabs( R_ecef_2_cg_ECEF[0] ) > SMALL_NUM ) {
    *longitude = atan2( R_ecef_2_cg_ECEF[1] , R_ecef_2_cg_ECEF[0] ) ;

  } else {
    *longitude = M_PI / 2;
  }

  if ( *longitude < 0.0 ) {
    *longitude += (2.0 * M_PI);
  }

  *h = (R - *semimajor_axis * t)*cos( *lat ) + (R_ecef_2_cg_ECEF[2] - b)*sin( *lat );

  return 0;

} /* ---------- end of function geocentric_to_geodetic ----------*/

/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           cart2kep
//  Purpose:        Computes the Keplerian elements based on the cartesian inputs
//  Assumptions:    Oscullating elements only.
//  References      BMW  // Note by CBV: all the following equations can also be found in Vallado3 p. 104 to 108.
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/09/2015    |   ---     | Initial Implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int cart2kep( ORBITAL_ELEMENTS_T   *OE,
              double                r[3],
              double                v[3],
              double                time,
              double                u)

{
  // Declarations
  double rrss;
  double vrss;
  double unit_r[3];
  double rdotv;
  double h[3];
  double mag_h;
  double K[3] = {0.0};
  double node[3];
  double mag_node;
  double unit_node[3];
  double specific_energy;
  double tempv1[3];
  double tempv2[3];
  double coeff;
  double e_vector[3];
  double unit_e[3];
  double e_vector_dot_r;
  double n;
    
  // Begin algorithm
  v_mag( &rrss, r);
  v_mag( &vrss, v);
  v_norm( unit_r, r);
  v_dot(&rdotv, r, v);

  //
  //  Solve for angular momentum (r x v)
  //
  v_cross(h,r,v);
  v_mag( &mag_h, h);

  //
  // Solve for node vector
  //
  K[2] = 1.0;
  v_cross(node, K,h);
  v_mag(&mag_node, node);
  v_norm( unit_node, node);

  //
  //  Solve for semi-major axis, sma
  //
  OE->sma = 1.0 / ( (2.0/rrss) - ( (vrss*vrss)/u ) );
    
  //
  //  Solve for ecentricity, e and the e vector
  //
  specific_energy = -1.0*(u/(2.0*OE->sma));
  v_scale(tempv1, v, rdotv);
    
  coeff = vrss*vrss - (u/rrss);
    
  v_scale(tempv2, r, coeff);
  v_sub( e_vector , tempv2, tempv1);
    
  coeff = 1.0/u;
  v_scale(e_vector, e_vector, (coeff));

  v_mag( &OE->eccentricity, e_vector );
  v_norm( unit_e, e_vector);

  //  Solve for inclination, i
  OE->inclination = acos(h[2]/mag_h);

  //  Solve for longitude of ascending node
  if (mag_node == 0.0) {
    
    OE->long_an = 0.0;
    
  } else if (node[1] >= 0){
    
    // TODO: Check this
    OE->long_an = acos(node[0]/mag_node);  // was checked by CBV (Vallado3 eq (2.82))
        
  }

  else if (node[1] < 0){
    OE->long_an = 2*M_PI - acos(node[0]/mag_node);
  }

  //
  //  Solve for argument of periapse
  //
  if (mag_node != 0.0) {
    v_dot(&coeff, unit_node, unit_e);
    OE->w = acos( coeff );
    if (e_vector[2] < 0.0) {
      OE->w = (2.0*M_PI - OE->w);
    }/*  else { */
    /*   OE->w = 0; */
    /* } */
  }

  //  Solve for true anomaly
  v_dot(&e_vector_dot_r, e_vector,r);
  if (OE->eccentricity != 0) {
    OE->f = acos(e_vector_dot_r/(OE->eccentricity*rrss));
    if (rdotv < 0) {
      OE->f = (2*M_PI) - fabs(OE->f);
      //OE->f = OE->f + M_PI;
    }
  } else {
    OE->f = 0;
  }
  //
  //  Solve for time of periapsis
  //
  OE->E = 2*atan(sqrt((1-OE->eccentricity)/(1+OE->eccentricity))*tan(OE->f/2));
  if (OE->E < 0) {
    OE->E = OE->E + (2*M_PI);
  }

  n = sqrt(u/(OE->sma*OE->sma*OE->sma)); // mean motion - Vallado3 eq (2.76) (CBV)
  OE->tp = -1 * ((OE->E/n) - ((OE->eccentricity*sin(OE->E))/n) - time);

  // Right ascension
  OE->ra =  atan2(r[1], r[0]);
  if ( OE->ra < 0){
    OE->ra = 2*M_PI + OE->ra ;
  }


  // orbital period
  OE->period = pow( OE->sma, 3.0);
  OE->period = OE->period / u;
  OE->period = 2.0 * M_PI * sqrt( OE->period );

  // AN to sc
  OE->an_to_sc = fmod(OE->w + OE->f, 2*M_PI); // angle ascending node to s/c


  // Solar zenith
  double x[6];
  double lt;
  double r_earth2sun_J2000[3];
  double r_cg2sun_J2000[3];
  double r_cg2sun_J2000_normalized[3];
    spkez_c(10, time, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.

    r_earth2sun_J2000[0] = x[0];
    r_earth2sun_J2000[1] = x[1];
    r_earth2sun_J2000[2] = x[2];
    double r_i2cg_INRTL_normalized[3];
    v_sub(r_cg2sun_J2000, r_earth2sun_J2000, r);
    v_norm(r_cg2sun_J2000_normalized, r_cg2sun_J2000);
    v_norm(r_i2cg_INRTL_normalized, r);
    double cos_earthtosc_sctosun;
    v_dot(&cos_earthtosc_sctosun, r_i2cg_INRTL_normalized, r_cg2sun_J2000_normalized);

    double cos_v_sctosun;
    double v_i2cg_INRTL_normalized[3];
    v_norm(v_i2cg_INRTL_normalized, v);
    v_dot(&cos_v_sctosun, v_i2cg_INRTL_normalized, r_cg2sun_J2000_normalized);
    if (cos_v_sctosun > 0){
    OE->zenith = 2*M_PI - acos(cos_earthtosc_sctosun);
    }
    else{
      OE->zenith = acos(cos_earthtosc_sctosun);
    }
  return 0;
    
} /* ---------- end of function cart2kep ----------*/

/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           kep2cart
//  Purpose:        Computes the ECI coordinates based on the Keplerian inputs
//  Assumptions:    None.
//  References      BMW
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/09/2015    |   ---     | Initial Implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int kep2cart(   double              r_i2cg_INRTL[3],
                double              v_i2cg_INRTL[3],
                double             *mu,
                ORBITAL_ELEMENTS_T *OE)

{

  /* Declarations */
  double slr;
  double rm;
  double arglat;
  double sarglat;
  double carglat;
  double c4;
  double c5;
  double c6;
  double sinc;
  double cinc;
  double sraan;
  double craan;
    
  /* Algorithm */
  slr = OE->sma * (1 - OE->eccentricity * OE->eccentricity);

  rm = slr / (1 + OE->eccentricity * cos(OE->f));
   
  arglat = OE->w + OE->f;

  sarglat = sin(arglat);
  carglat = cos(arglat);
   
  c4 = sqrt(mu[0] / slr);
  c5 = OE->eccentricity * cos(OE->w) + carglat;
  c6 = OE->eccentricity * sin(OE->w) + sarglat;

  sinc = sin(OE->inclination);
  cinc = cos(OE->inclination);

  sraan = sin(OE->long_an);
  craan = cos(OE->long_an);

  // position vector
  r_i2cg_INRTL[0] = rm * (craan * carglat - sraan * cinc * sarglat);
  r_i2cg_INRTL[1] = rm * (sraan * carglat + cinc * sarglat * craan);
  r_i2cg_INRTL[2] = rm * sinc * sarglat;

  // velocity vector
  v_i2cg_INRTL[0] = -c4 * (craan * c6 + sraan * cinc * c5);
  v_i2cg_INRTL[1] = -c4 * (sraan * c6 - craan * cinc * c5);
  v_i2cg_INRTL[2] = c4 * c5 * sinc;

  return 0;
} /* ---------- end of function kep2cart ----------*/

/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           propagate_spacecraft
//  Purpose:        RK4
//  Assumptions:    None.
//  References      BMW
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/09/2015    |   ---     | Initial Implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int propagate_spacecraft(   SPACECRAFT_T *SC,
                            PARAMS_T     *PARAMS,
			    double et_initial_epoch,
			    double et_sc_initial,
			    double *density,
			    GROUND_STATION_T *GROUND_STATION,
			    OPTIONS_T *OPTIONS,
			    CONSTELLATION_T *CONSTELLATION,
			    int iProc,
			    int iDebugLevel,
			    int *start_ensemble,
			    int *array_sc)

{
  //    etprint(SC->et, "before propagate");


  
  //      printf("\n");
  // Declarations
  //  int start_ensemble;
  double geodetic[3];
  SpiceDouble       xform[6][6];
  double estate[6], jstate[6];
  double k1r[3];
  double k1v[3];
  double k2r[3];
  double k2v[3];
  double k3r[3];
  double k3v[3];
  double k4r[3];
  double k4v[3];
  double dr[3];
  double dv[3];
  double etk;
  double rk[3];
  double vk[3];
  double sc_ecef_previous_time_step[3];  double sc_eci_previous_time_step[3];  double sc_v_eci_previous_time_step[3];
  double previous_an_to_sc = SC->OE.an_to_sc; // for orbit average computation
/*   /\*************************************************************************************\/ */
/*   /\*************************************************************************************\/ */
/*   /\*********************** SET THINGS UP FOR PARALLEL PROGRAMMING **********************\/ */
/*   /\*************************************************************************************\/ */
/*   /\*************************************************************************************\/ */
/*   /\*   Things we set up here: *\/ */
/*   /\*   - which iProc runs which main sc / which sc is run by which iProc *\/ */


/*   // For each iProc, set up the first main sc (iStart_save[iProcf]) and the last main sc (iEnd_save[iProcf]) that it's going to run. iStart_save and iEnd_save are two arrays that have the same values for all procs (so if you're iProc 0 or iProc 1, you have value recorded for iStart_save[0] and iEnd_save[0] and the same value recorded for iStart_save[1] and iEnd_save[1]) -> they are not "iProc-dependent" */
/*   int *iStart_save, *iEnd_save; */
/*   int nscEachPe, nscLeft; */
/*   int iProcf; */
/*   int i; */
/*   nscEachPe = (OPTIONS->n_satellites)/nProcs; */
/*   nscLeft = (OPTIONS->n_satellites) - (nscEachPe * nProcs); */

/*   iStart_save = malloc( nProcs * sizeof(int)); */
/*   iEnd_save = malloc( nProcs  * sizeof(int)); */
/*   for (iProcf = 0; iProcf < nProcs; iProcf++){ */
/*     iStart_save[iProcf] = 0; */
/*     iEnd_save[iProcf] = 0; */
/*   } */
/*   for (iProcf = 0; iProcf < nProcs; iProcf++){ */
/*     for (i=0; i<iProcf; i++) { */
/*       iStart_save[iProcf] += nscEachPe; */
/*       if (i < nscLeft && iProcf > 0) iStart_save[iProcf]++; */
/*     } */
/*     iEnd_save[iProcf] = iStart_save[iProcf]+nscEachPe; */
/*     if (iProcf  < nscLeft) iEnd_save[iProcf]++; */
/*     iStart_save[iProcf] = iStart_save[iProcf]; */
/*     iEnd_save[iProcf] = iEnd_save[iProcf]; */
/*   } */
    
/*   //    if (iProc == 0){ */
/*   /\*     for (iProcf = 0; iProcf < nProcs; iProcf++){ *\/ */
/*   /\*       printf("%d - %d\n", iStart_save[iProcf], iEnd_save[iProcf] - 1) ; *\/ */
/*   /\*     } *\/ */
/*   //} */

/*   // For each main sc, start_ensemble is 0 if the iProc runs this main sc. start_ensemble is 1 is the iProc does not run this main sc. (so each iProc has a different array start_ensemble -> start_ensemble is "iProc-dependent")  */
/*   int *start_ensemble; */
/*   start_ensemble = malloc(OPTIONS->n_satellites * sizeof(int)); */
/*   for (ii = 0; ii < OPTIONS->n_satellites; ii++){ */
/*     if ( (ii >= iStart_save[iProc]) & ( ii < iEnd_save[iProc]) ){ */
/*       start_ensemble[ii] = 0; */
/*     } */
/*     else{ */
/*       start_ensemble[ii] = 1; */
/*     } */
/*     //    printf("iProc %d | start_ensemble[%d] %d\n", iProc, ii, start_ensemble[ii]); */
/*   } */

/* /\*   if ( (iProc == 0) && (OPTIONS->first_run_to_find_tca_before_collision_assessment == 1)){ *\/ */
/* /\*     start_ensemble = 0; *\/ */
/* /\*   } *\/ */
/* /\*   else{ *\/ */
/* /\*     start_ensemble = 1; *\/ */
/* /\*   } *\/ */

/*   /\*************************************************************************************\/ */
/*   /\*************************************************************************************\/ */
/*   /\****************** end of SET THINGS UP FOR PARALLEL PROGRAMMING ********************\/ */
/*   /\*************************************************************************************\/ */
/*   /\*************************************************************************************\/ */

//  SC->INTEGRATOR.write_given_output = 0; commented on 073119


  if  (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) != 0 ){ // if sgp4 equations are not used for the orbit propagation
  // Compute k
  v_copy(sc_ecef_previous_time_step, SC->r_ecef2cg_ECEF);
  v_copy(sc_eci_previous_time_step, SC->r_i2cg_INRTL);
  v_copy(sc_v_eci_previous_time_step, SC->v_i2cg_INRTL);

  // the only case where the attitute needs to be set is the only case where it has not been set initially in initialize_constellation, which is if there is no ensemble at all run on the attitude
  if ( ( strcmp(OPTIONS->attitude_profile, "ensemble_angular_velocity") != 0 ) && ( strcmp(OPTIONS->attitude_profile, "ensemble_initial_attitude") != 0 ) ) { // if we do not run ensembles on the initial angular velocity

    if (OPTIONS->nb_ensembles_attitude <= 0){ // if we dont' run ensembles of any kind on the attitude

      if ( (SC->INTEGRATOR.solar_cell_efficiency != -1) || (GROUND_STATION->nb_ground_stations > 0) || (SC->INTEGRATOR.include_solar_pressure == 1) || (SC->INTEGRATOR.include_drag == 1) ){ // these are the case where SpOCK uses the attitude

	if (SC->INTEGRATOR.isGPS == 0){ // no attitude to set up for GPS because we dont compute the drag, solar rariatin pressu,re, power, and gorund station coverage (attitude is needed for coverage because we calculate azimuth and levation angles form the spaceecraft reference system)
	  
	  set_attitude(SC->INTEGRATOR.attitude, SC->INTEGRATOR.index_in_attitude_interpolated, OPTIONS, SC->INTEGRATOR.file_is_quaternion);

	}
      }
    }
  }
  // end of if the only case where the attitute needs to be set is the only case where it has not been set initially in initialize_constellation, which is if there is no ensemble at all run on the attitude
  SC->INTEGRATOR.last_compute_dxdt = 0; // dont worry about it


  compute_dxdt( k1r, k1v, &SC->et, SC->r_i2cg_INRTL, SC->v_i2cg_INRTL, PARAMS, &SC->INTEGRATOR, et_initial_epoch, et_sc_initial, density, SC->INTEGRATOR.index_in_attitude_interpolated, SC->INTEGRATOR.index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc, iDebugLevel, SC);
  v_scale(k1r, k1r, SC->INTEGRATOR.dt);
  v_scale(k1v, k1v, SC->INTEGRATOR.dt);

  // Compute k2
  etk = SC->et + SC->INTEGRATOR.dt_pos_neg / 2.0;
    
  v_copy( dr, k1r);
  v_scale(dr, dr, 0.5);    
  v_copy( dv, k1v);
  v_scale(dv, dv, 0.5);
  /* v_print(k1r, "k1r"); */
  /* v_print(k1v, "k1v"); */
  /* exitf(); */
  if (SC->INTEGRATOR.dt_pos_neg >=0){
    v_add(rk, SC->r_i2cg_INRTL, dr);
    v_add(vk, SC->v_i2cg_INRTL, dv);
  }
  else{// backward propagation
    v_sub(rk, SC->r_i2cg_INRTL, dr);
    v_sub(vk, SC->v_i2cg_INRTL, dv);
  }
  
  if ( ( SC->INTEGRATOR.include_drag != 0 ) || ( SC->INTEGRATOR.include_solar_pressure != 0 ) || ( SC->INTEGRATOR.solar_cell_efficiency != -1 ) || (GROUND_STATION->nb_ground_stations > 0) ){
    //    if (SC->et >= et_initial_epoch){ //  that's because if the initializtion of the orbit was made with a TLE, we don't want index_in_attitude_interpolated to be incremented for times between the TLE epoch and the inital epoch 
      SC->INTEGRATOR.index_in_attitude_interpolated =  SC->INTEGRATOR.index_in_attitude_interpolated + 1;
    
      if ( SC->INTEGRATOR.include_drag != 0 ){
	SC->INTEGRATOR.index_in_driver_interpolated =  SC->INTEGRATOR.index_in_driver_interpolated + 1;
      }
      //    }
  }

  if ( (OPTIONS->nb_ensembles_density > 0) && (OPTIONS->swpc_need_predictions) && ( SC->et >= OPTIONS->swpc_et_first_prediction ) ) { 
	if( ( (array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] >= 0) ) && (SC->INTEGRATOR.sc_ensemble_nb == array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] )){ // used to be	if( ( (array_sc[start_ensemble[SC->INTEGRATOR.sc_ensemble_nb]] >= 0) ) && (SC->INTEGRATOR.sc_ensemble_nb == array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] )){  but doesn't make senese because start_ensemble is defined only for the number of main sc so start_ensemble[xxx] doesnt exist for xxx = SC->INTEGRATOR.sc_ensemble_nb (imagine there's 2 main sc and 10 ensbles then start_ensemble is defined for xxx = 0 or 1, not 2, 3, ..., 10. The explanation right below was for this if that used to be here so not sure it's all correct either
      // if this iProc is running main sc SC->INTEGRATOR.sc_main_nb then start_ensemble[SC->INTEGRATOR.sc_main_nb] = 0. If this iProc is not running main sc SC->INTEGRATOR.sc_main_nb but is running an ensemble corresponding to main sc SC->INTEGRATOR.sc_main_nb, then start_ensemble[SC->INTEGRATOR.sc_main_nb] = 1. If none of these two, then start_ensemble[SC->INTEGRATOR.sc_main_nb] = -1. So here we don't want to get in this if the iProc does not run main sc or an ensemble corresponding to this main sc (if (array_sc[start_ensemble[SC->INTEGRATOR.sc_ensemble_nb]] >= 0)). And we want to get in only once, for the first sc that this iProc is running

      CONSTELLATION->aaa_sigma[SC->INTEGRATOR.sc_main_nb] = CONSTELLATION->aaa_sigma[SC->INTEGRATOR.sc_main_nb] + 1;
    }
  }
  if ( (strcmp(OPTIONS->test_omniweb_or_external_file, "swpc_mod") == 0) && (OPTIONS->swpc_need_predictions) && ( SC->et >= OPTIONS->swpc_et_first_prediction ) ){
	if( ( (array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] >= 0) ) && (SC->INTEGRATOR.sc_ensemble_nb == array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] )){ // used to be	if( ( (array_sc[start_ensemble[SC->INTEGRATOR.sc_ensemble_nb]] >= 0) ) && (SC->INTEGRATOR.sc_ensemble_nb == array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] )){  but doesn't make senese because start_ensemble is defined only for the number of main sc so start_ensemble[xxx] doesnt exist for xxx = SC->INTEGRATOR.sc_ensemble_nb (imagine there's 2 main sc and 10 ensbles then start_ensemble is defined for xxx = 0 or 1, not 2, 3, ..., 10. The explanation right below was for this if that used to be here so not sure it's all correct either
      // if this iProc is running main sc SC->INTEGRATOR.sc_main_nb then start_ensemble[SC->INTEGRATOR.sc_main_nb] = 0. If this iProc is not running main sc SC->INTEGRATOR.sc_main_nb but is running an ensemble corresponding to main sc SC->INTEGRATOR.sc_main_nb, then start_ensemble[SC->INTEGRATOR.sc_main_nb] = 1. If none of these two, then start_ensemble[SC->INTEGRATOR.sc_main_nb] = -1. So here we don't want to get in this if the iProc does not run main sc or an ensemble corresponding to this main sc (if (array_sc[start_ensemble[SC->INTEGRATOR.sc_ensemble_nb]] >= 0)). And we want to get in only once, for the first sc that this iProc is running

      CONSTELLATION->aaa_mod[SC->INTEGRATOR.sc_main_nb] = CONSTELLATION->aaa_mod[SC->INTEGRATOR.sc_main_nb] + 1;

    }
  }



  // the only case where the attitute needs to be set is the only case where it has not been set initially in initialize_constellation, which is if there is no ensemble at all run on the attitude
  if ( ( strcmp(OPTIONS->attitude_profile, "ensemble_angular_velocity") != 0 ) && ( strcmp(OPTIONS->attitude_profile, "ensemble_initial_attitude") != 0 ) ) { // if we do not run ensembles on the initial angular velocity
    if (OPTIONS->nb_ensembles_attitude <= 0){
      if ( (SC->INTEGRATOR.solar_cell_efficiency != -1) || (GROUND_STATION->nb_ground_stations > 0) || (SC->INTEGRATOR.include_solar_pressure == 1) || (SC->INTEGRATOR.include_drag == 1) ){ // these are the case where SpOCK uses the attitude
	if (SC->INTEGRATOR.isGPS == 0){ // no attitude to set up for GPS because we dont compute the drag, solar rariatin pressu,re, power, and gorund station coverage (attitude is needed for coverage because we calculate azimuth and levation angles form the spaceecraft reference system)
	  set_attitude(SC->INTEGRATOR.attitude, SC->INTEGRATOR.index_in_attitude_interpolated, OPTIONS, SC->INTEGRATOR.file_is_quaternion);
	}
      }
    }
  }
  // end of if the only case where the attitute needs to be set is the only case where it has not been set initially in initialize_constellation, which is if there is no ensemble at all run on the attitude

  compute_dxdt( k2r, k2v, &etk, rk, vk, PARAMS, &SC->INTEGRATOR, et_initial_epoch, et_sc_initial,density, SC->INTEGRATOR.index_in_attitude_interpolated, SC->INTEGRATOR.index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc,iDebugLevel, SC);


  v_scale(k2r, k2r, SC->INTEGRATOR.dt);
  v_scale(k2v, k2v, SC->INTEGRATOR.dt);

  // Compute k3
  etk = SC->et + SC->INTEGRATOR.dt_pos_neg / 2.0;
    
  v_copy( dr, k2r);
  v_scale(dr, dr, 0.5);

  v_copy( dv, k2v);
  v_scale(dv, dv, 0.5);
  if (SC->INTEGRATOR.dt_pos_neg >=0){
      v_add(rk, SC->r_i2cg_INRTL, dr);
      v_add(vk, SC->v_i2cg_INRTL, dv);
  }
  else{// backward propagation
    v_sub(rk, SC->r_i2cg_INRTL, dr);
    v_sub(vk, SC->v_i2cg_INRTL, dv);
  }

  
  compute_dxdt( k3r, k3v, &etk, rk, vk, PARAMS, &SC->INTEGRATOR, et_initial_epoch, et_sc_initial, density, SC->INTEGRATOR.index_in_attitude_interpolated, SC->INTEGRATOR.index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc, iDebugLevel, SC);

  v_scale(k3r, k3r, SC->INTEGRATOR.dt);
  v_scale(k3v, k3v, SC->INTEGRATOR.dt);
    
  // Compute k4
  etk = SC->et + SC->INTEGRATOR.dt_pos_neg;

  v_copy( dr, k3r);
  v_copy( dv, k3v);


  if (SC->INTEGRATOR.dt_pos_neg >=0){
      v_add(rk, SC->r_i2cg_INRTL, dr);
      v_add(vk, SC->v_i2cg_INRTL, dv);
  }
  else{// backward propagation
      v_sub(rk, SC->r_i2cg_INRTL, dr);
      v_sub(vk, SC->v_i2cg_INRTL, dv);
  }

  
  if ( ( SC->INTEGRATOR.include_drag != 0 ) || ( SC->INTEGRATOR.include_solar_pressure != 0 ) || ( SC->INTEGRATOR.solar_cell_efficiency != -1 )  || (GROUND_STATION->nb_ground_stations > 0) ){
    //    if (SC->et >= et_initial_epoch){ //  that's because if the initializtion of the orbit was made with a TLE, we don't want index_in_attitude_interpolated to be incremented for times between the TLE epoch and the inital epoch 
      SC->INTEGRATOR.index_in_attitude_interpolated =  SC->INTEGRATOR.index_in_attitude_interpolated + 1;
      if ( SC->INTEGRATOR.include_drag != 0 ){
	SC->INTEGRATOR.index_in_driver_interpolated =  SC->INTEGRATOR.index_in_driver_interpolated + 1;
      }
      //    }
  }
  if ( (OPTIONS->nb_ensembles_density > 0) && (OPTIONS->swpc_need_predictions) && ( SC->et >= OPTIONS->swpc_et_first_prediction ) ) { 
	if( ( (array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] >= 0) ) && (SC->INTEGRATOR.sc_ensemble_nb == array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] )){ // used to be	if( ( (array_sc[start_ensemble[SC->INTEGRATOR.sc_ensemble_nb]] >= 0) ) && (SC->INTEGRATOR.sc_ensemble_nb == array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] )){  but doesn't make senese because start_ensemble is defined only for the number of main sc so start_ensemble[xxx] doesnt exist for xxx = SC->INTEGRATOR.sc_ensemble_nb (imagine there's 2 main sc and 10 ensbles then start_ensemble is defined for xxx = 0 or 1, not 2, 3, ..., 10. The explanation right below was for this if that used to be here so not sure it's all correct either
      // if this iProc is running main sc SC->INTEGRATOR.sc_main_nb then start_ensemble[SC->INTEGRATOR.sc_main_nb] = 0. If this iProc is not running main sc SC->INTEGRATOR.sc_main_nb but is running an ensemble corresponding to main sc SC->INTEGRATOR.sc_main_nb, then start_ensemble[SC->INTEGRATOR.sc_main_nb] = 1. If none of these two, then start_ensemble[SC->INTEGRATOR.sc_main_nb] = -1. So here we don't want to get in this if the iProc does not run main sc or an ensemble corresponding to this main sc (if (array_sc[start_ensemble[SC->INTEGRATOR.sc_ensemble_nb]] >= 0)). And we want to get in only once, for the first sc that this iProc is running

      CONSTELLATION->aaa_sigma[SC->INTEGRATOR.sc_main_nb] = CONSTELLATION->aaa_sigma[SC->INTEGRATOR.sc_main_nb] + 1;
    }
  }
  if ( (strcmp(OPTIONS->test_omniweb_or_external_file, "swpc_mod") == 0) && (OPTIONS->swpc_need_predictions) && ( SC->et >= OPTIONS->swpc_et_first_prediction ) ){
	if( ( (array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] >= 0) ) && (SC->INTEGRATOR.sc_ensemble_nb == array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] )){ // used to be	if( ( (array_sc[start_ensemble[SC->INTEGRATOR.sc_ensemble_nb]] >= 0) ) && (SC->INTEGRATOR.sc_ensemble_nb == array_sc[start_ensemble[SC->INTEGRATOR.sc_main_nb]] )){  but doesn't make senese because start_ensemble is defined only for the number of main sc so start_ensemble[xxx] doesnt exist for xxx = SC->INTEGRATOR.sc_ensemble_nb (imagine there's 2 main sc and 10 ensbles then start_ensemble is defined for xxx = 0 or 1, not 2, 3, ..., 10. The explanation right below was for this if that used to be here so not sure it's all correct either
      // if this iProc is running main sc SC->INTEGRATOR.sc_main_nb then start_ensemble[SC->INTEGRATOR.sc_main_nb] = 0. If this iProc is not running main sc SC->INTEGRATOR.sc_main_nb but is running an ensemble corresponding to main sc SC->INTEGRATOR.sc_main_nb, then start_ensemble[SC->INTEGRATOR.sc_main_nb] = 1. If none of these two, then start_ensemble[SC->INTEGRATOR.sc_main_nb] = -1. So here we don't want to get in this if the iProc does not run main sc or an ensemble corresponding to this main sc (if (array_sc[start_ensemble[SC->INTEGRATOR.sc_ensemble_nb]] >= 0)). And we want to get in only once, for the first sc that this iProc is running

      CONSTELLATION->aaa_mod[SC->INTEGRATOR.sc_main_nb] = CONSTELLATION->aaa_mod[SC->INTEGRATOR.sc_main_nb] + 1;
    }
  }

  // the only case where the attitute needs to be set is the only case where it has not been set initially in initialize_constellation, which is if there is no ensemble at all run on the attitude
  if ( ( strcmp(OPTIONS->attitude_profile, "ensemble_angular_velocity") != 0 ) && ( strcmp(OPTIONS->attitude_profile, "ensemble_initial_attitude") != 0 ) ) { // if we do not run ensembles on the initial angular velocity
    if (OPTIONS->nb_ensembles_attitude <= 0){
      if ( (SC->INTEGRATOR.solar_cell_efficiency != -1) || (GROUND_STATION->nb_ground_stations > 0) || (SC->INTEGRATOR.include_solar_pressure == 1) || (SC->INTEGRATOR.include_drag == 1) ){ // these are the case where SpOCK uses the attitude
	if (SC->INTEGRATOR.isGPS == 0){ // no attitude to set up for GPS because we dont compute the drag, solar rariatin pressu,re, power, and gorund station coverage (attitude is needed for coverage because we calculate azimuth and levation angles form the spaceecraft reference system)
	  set_attitude(SC->INTEGRATOR.attitude, SC->INTEGRATOR.index_in_attitude_interpolated, OPTIONS, SC->INTEGRATOR.file_is_quaternion);
	}
      }
    }
  }
  // end of if the only case where the attitute needs to be set is the only case where it has not been set initially in initialize_constellation, which is if there is no ensemble at all run on the attitude
  //  printf("%f %f %f %d %d %d\n", SC->INTEGRATOR.attitude.pitch[SC->INTEGRATOR.index_in_attitude_interpolated], SC->INTEGRATOR.attitude.roll[SC->INTEGRATOR.index_in_attitude_interpolated], SC->INTEGRATOR.attitude.yaw[SC->INTEGRATOR.index_in_attitude_interpolated], SC->INTEGRATOR.attitude.order_pitch[SC->INTEGRATOR.index_in_attitude_interpolated], SC->INTEGRATOR.attitude.order_roll[SC->INTEGRATOR.index_in_attitude_interpolated], SC->INTEGRATOR.attitude.order_yaw[SC->INTEGRATOR.index_in_attitude_interpolated]);
  compute_dxdt( k4r, k4v, &etk, rk, vk, PARAMS, &SC->INTEGRATOR, et_initial_epoch, et_sc_initial, density, SC->INTEGRATOR.index_in_attitude_interpolated, SC->INTEGRATOR.index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc, iDebugLevel, SC);

  v_scale(k4r, k4r, SC->INTEGRATOR.dt);
  v_scale(k4v, k4v, SC->INTEGRATOR.dt);
  
  //  Compute the delta r, delta v
  v_copy( dr, k4r);
  v_scale( k3r, k3r, 2.0);
  if (SC->INTEGRATOR.dt_pos_neg >=0){
    v_add(dr, dr, k3r);    
  }
  else{// backward propagation 
    v_add(dr, dr, k3r);
  }

  v_scale( k2r, k2r, 2.0);
  if (SC->INTEGRATOR.dt_pos_neg >=0){
      v_add(dr, dr, k2r);
      v_add(dr, dr, k1r);

  }
  else{// backward propagation
      v_add(dr, dr, k2r);
      v_add(dr, dr, k1r);
  }

  v_scale( dr, dr, (1.0/6.0));
    
  v_copy( dv, k4v);
  v_scale( k3v, k3v, 2.0);
  if (SC->INTEGRATOR.dt_pos_neg >=0){
    v_add(dv, dv, k3v);
  }
  else{// backward propagation
    v_add(dv, dv, k3v);
  }
  v_scale( k2v, k2v, 2.0);

    if (SC->INTEGRATOR.dt_pos_neg >=0){
      v_add(dv, dv, k2v);
  v_add(dv, dv, k1v);

  }
  else{// backward propagation
    v_add(dv, dv, k2v);
  v_add(dv, dv, k1v);
  }
  v_scale( dv, dv, (1.0/6.0));

  
  // Update Inertial State
  SC->et = SC->et + SC->INTEGRATOR.dt_pos_neg;

    if (SC->INTEGRATOR.dt_pos_neg >=0){
        v_add( SC->r_i2cg_INRTL, SC->r_i2cg_INRTL, dr);
  v_add( SC->v_i2cg_INRTL, SC->v_i2cg_INRTL, dv);

  }
  else{// backward propagation
    v_sub( SC->r_i2cg_INRTL, SC->r_i2cg_INRTL, dr);// the dr is the sum of k1, k2, k3, k4 for either the fowrward or backward propagation so for backward propagation we need to remove this entire dr
  v_sub( SC->v_i2cg_INRTL, SC->v_i2cg_INRTL, dv);
  }



  double starttime;
  str2et_c(OPTIONS->initial_epoch, &starttime);
  double min_end_time;
  str2et_c(OPTIONS->final_epoch, &min_end_time);

  double dvdt[3], drdt[3];
	    if ( ( fmod( SC->et - starttime, OPTIONS->dt_output ) < OPTIONS->dt / 2.) || ( fabs( fmod( SC->et - starttime, OPTIONS->dt_output ) - OPTIONS->dt_output ) < OPTIONS->dt / 2.) || ( SC->et > min_end_time - 0.01)  )  {

	      //  SC->INTEGRATOR.write_given_output = 0; commented on 073119  
	    }
	    //	    printf("\n\n");
	    SC->INTEGRATOR.last_compute_dxdt = 1;

	    compute_dxdt( drdt, dvdt, &SC->et,  SC->r_i2cg_INRTL, SC->v_i2cg_INRTL, PARAMS, &SC->INTEGRATOR, et_initial_epoch, et_sc_initial, density, SC->INTEGRATOR.index_in_attitude_interpolated, SC->INTEGRATOR.index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc, iDebugLevel, SC);
	    SC->INTEGRATOR.last_compute_dxdt = 0;
  v_copy( SC->a_i2cg_INRTL,  dvdt);

  
  // For Kalman: convert acceleration from inertial to lvlh
  double T_inrtl_2_lvlh[3][3];
  compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, SC->r_i2cg_INRTL, SC->v_i2cg_INRTL);
  m_x_v(SC->a_i2cg_LVLH, T_inrtl_2_lvlh, SC->a_i2cg_INRTL);
  m_x_v(SC->a_i2cg_LVLH_gravity, T_inrtl_2_lvlh, SC->a_i2cg_INRTL_gravity);
  m_x_v(SC->a_i2cg_LVLH_drag, T_inrtl_2_lvlh, SC->a_i2cg_INRTL_drag);

} // end of if sgp4 eqautions are not used for the orbit propagation
  else{ // if sgp4 eqautions are not used for the orbit propagation
    	    // !!!!!!!!!!!!! TEMP remove 080619
    SpiceDouble state[6];
    SC->et = SC->et + SC->INTEGRATOR.dt_pos_neg;
    	    extern /* Subroutine */ int ev2lin_(SpiceDouble *, SpiceDouble *,
	    					SpiceDouble *, SpiceDouble *);

    ev2lin_( &SC->et, PARAMS->geophs, SC->INTEGRATOR.elems, state );
    /* etprint(SC->et, ""); */
    /* 			    int lll; */
    /* 			    for (lll = 0; lll < 6; lll++){ */
    /* 			      printf("state[%d]: %.8f\n", lll, state[lll]); */
    /* 			    } */
    /* 			    //			    exitf(); */
    /* 			    // end of TEMP remove 080619 */
	    static doublereal  precm[36];
	        static doublereal invprc[36]	/* was [6][6] */;
	    			    extern /* Subroutine */ int zzteme_(doublereal *, doublereal *);
	    	        extern /* Subroutine */ int invstm_(doublereal *, doublereal *);
	    		    extern /* Subroutine */ int  mxvg_(
	    doublereal *, doublereal *, integer *, integer *, doublereal *);

	        zzteme_(&SC->et, precm);

/*     ...now convert STATE to J2000. Invert the state transformation */
/*     operator (important to correctly do this). */

    invstm_(precm, invprc);
    static integer c__6 = 6;
        static doublereal tmpsta[6];
	int pp;
	    mxvg_(invprc, state, &c__6, &c__6, tmpsta);
    moved_(tmpsta, &c__6, state);
	    for (pp = 0; pp<3; pp++){
	      //	      	      	      printf("%f ", state[pp]);
	      SC->r_i2cg_INRTL[pp] = state[pp];
	      SC->v_i2cg_INRTL[pp] = state[pp+3];
	    }

  } // enf of if sgp4 eqautions are not used for the orbit propagation
 
  // Update ECEF state
  estate[0] = SC->r_i2cg_INRTL[0];estate[1] = SC->r_i2cg_INRTL[1];estate[2] = SC->r_i2cg_INRTL[2];
  estate[3] = SC->v_i2cg_INRTL[0];estate[4] = SC->v_i2cg_INRTL[1];estate[5] = SC->v_i2cg_INRTL[2];
  sxform_c (  "J2000", PARAMS->EARTH.earth_fixed_frame,  SC->et,    xform  );

  mxvg_c   (  xform,       estate,   6,  6, jstate );
  SC->r_ecef2cg_ECEF[0] = jstate[0]; SC->r_ecef2cg_ECEF[1] = jstate[1]; SC->r_ecef2cg_ECEF[2] = jstate[2];
  SC->v_ecef2cg_ECEF[0] = jstate[3]; SC->v_ecef2cg_ECEF[1] = jstate[4]; SC->v_ecef2cg_ECEF[2] = jstate[5];

   // Update Geodetic State
  geocentric_to_geodetic( SC->r_ecef2cg_ECEF,
			  &PARAMS->EARTH.radius,
			  &PARAMS->EARTH.flattening,
			  &SC->GEODETIC.altitude,
			  &SC->GEODETIC.latitude,
			  &SC->GEODETIC.longitude);



/* /\*    // Update Geodetic State *\/ */
/* 	eci2lla(SC->r_i2cg_INRTL , SC->et, geodetic ); */
	
/* 	SC->GEODETIC.altitude = geodetic[2]; */
/* 	SC->GEODETIC.latitude = geodetic[0]; */
/* 	SC->GEODETIC.longitude = geodetic[1];  */



/* 	// update the planet fixed state		 */
/* 	  geodetic_to_geocentric(PARAMS->EARTH.flattening,             */
/* 				 SC->GEODETIC.altitude, */
/* 				 SC->GEODETIC.latitude, */
/* 				 SC->GEODETIC.longitude, */
/* 				 PARAMS->EARTH.radius,        */
/* 				 SC->r_ecef2cg_ECEF) ; */

  // Update Keplerian State
  cart2kep(&SC->OE, SC->r_i2cg_INRTL, SC->v_i2cg_INRTL, SC->et, PARAMS->EARTH.GRAVITY.mu);

  orbit_ave(SC, CONSTELLATION, iProc, iDebugLevel, previous_an_to_sc);


  // Returns if the satellite is or not in the shadow of the Earth
  shadow_light( SC->INTEGRATOR.shadow, SC->r_i2cg_INRTL, SC->et, PARAMS);

  // Returns if the satellite is or not in the shadow of the Moon
  shadow_light_moon( SC->INTEGRATOR.shadow_moon, SC->r_i2cg_INRTL, SC->et, PARAMS);


  // Compute power from solar arrays

  if (SC->INTEGRATOR.solar_cell_efficiency != -1){
    compute_power(&SC->INTEGRATOR, CONSTELLATION, SC->r_i2cg_INRTL, SC->v_i2cg_INRTL, &SC->et, PARAMS, et_initial_epoch, et_sc_initial,SC->INTEGRATOR.index_in_attitude_interpolated);
  }  

  // printf("%d\n", SC->INTEGRATOR.index_in_attitude_interpolated);
  if (SC->INTEGRATOR.isGPS == 0){
    if (GROUND_STATION->nb_ground_stations > 0){
      if ((SC->INTEGRATOR.index_in_attitude_interpolated > 0) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 )){
	if (SC->et >= et_initial_epoch){
	coverage_ground_station( SC, GROUND_STATION, PARAMS, SC->INTEGRATOR.index_in_attitude_interpolated, &SC->INTEGRATOR, et_initial_epoch, et_sc_initial,sc_ecef_previous_time_step, sc_eci_previous_time_step, sc_v_eci_previous_time_step, 1 ); // last argument is the time step in seconds for the linear interpolation of the sc position when computing the coverage of the ground stations by the sc
	}
      }
    }
  }

  //    etprint(SC->et, "after propagate");
  //  printf("//\n\n");
  return 0;
}  /* ---------- end of function propagate_spacecraft ----------*/

/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           compute_dxdt
//  Purpose:        Computes the dxdt for position / velocity in RK4 integrator
//  Assumptions:    None
//  References      BMW
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/09/2015    |   ---     | Initial Implementation
//      | J. Getchius   | 05/20/2015    |   ---     | Update with MSIS and EGM96 gravity
//      | C. Bussy-Virat| 08/02/2015    |   ---     | Sun and Moon perturbations: correction of the direct effect (Vallado3 page 35 (eq (1-35)))
//
/////////////////////////////////////////////////////////////////////////////////////////
int compute_dxdt(   double          drdt[3],
                    double          dvdt[3],
                    double          *et,
                    double          r_i2cg_INRTL[3],
                    double          v_i2cg_INRTL[3],
                    PARAMS_T        *PARAMS,
                    INTEGRATOR_T    *INTEGRATOR,
		    double          et_initial_epoch, 
		    double          et_sc_initial, 
		    double          *density,
		    int             index_in_attitude_interpolated, 
		    int             index_in_driver_interpolated,
		    CONSTELLATION_T *CONSTELLATION,
		    OPTIONS_T       *OPTIONS,
		    int iProc,
		    int iDebugLevel,
		    SPACECRAFT_T *SC)




{

  // Declarations
  double r;
  double r3;
  double r_earth2moon_J2000[3];
  double x[6];
  double lt;
  double pert_direct[3];
  double pert_indirect[3];
  double pert_total[3];
  double coeff;
  double r_cg2moon_J2000[3];
  double r_cg2sun_J2000[3];
  double r_earth2sun_J2000[3];
  double a_solar_pressure_INRTL[3];
  double adrag_i2cg_INRTL[3];
      double athrust_i2cg_INRTL[3];
  double a_earth_pressure_INRTL[3];


  /* // drdt */
  v_copy( drdt, v_i2cg_INRTL);

  // dvdt

  // Gravity
  compute_gravity( dvdt, r_i2cg_INRTL, et[0], &PARAMS->EARTH.GRAVITY, INTEGRATOR->degree, PARAMS->EARTH.earth_fixed_frame, PARAMS->EARTH.flattening, PARAMS->EARTH.radius, SC, OPTIONS->gravity_map, CONSTELLATION);

  /*   if (dvdt[0] != aa){ */
  /*     print_test(); */
  /*     exitf(); */
  /*   } */
  //  Moon Perturbations
  if (INTEGRATOR->include_moon == 1){

    spkez_c(301, et[0], "J2000", "NONE", 399, x, &lt);
    r_earth2moon_J2000[0] = x[0];
    r_earth2moon_J2000[1] = x[1];
    r_earth2moon_J2000[2] = x[2];
    
    v_mag( &r, r_earth2moon_J2000);
    r3 = pow( r, 3.0 );
    coeff = PARAMS->MOON.GRAVITY.mu / r3;
    v_copy(pert_indirect, r_earth2moon_J2000);
    v_scale(pert_indirect, pert_indirect, coeff);
    
    v_sub(r_cg2moon_J2000, r_earth2moon_J2000, r_i2cg_INRTL);
    v_mag( &r, r_cg2moon_J2000);
    r3 = pow( r, 3.0 );
    coeff = PARAMS->MOON.GRAVITY.mu / r3;
    v_copy(pert_direct, r_cg2moon_J2000); // used to be: v_copy(pert_direct, r_earth2moon_J2000) but corrected based on Vallado3 page 35 equation (1-35)
    v_scale(pert_direct, pert_direct, coeff);
    v_sub(pert_total, pert_direct , pert_indirect);
    v_add(dvdt, dvdt , pert_total);
  }

  /* Sun Perturbations */
  if (INTEGRATOR->include_sun == 1){

    spkez_c(10, et[0], "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.

    r_earth2sun_J2000[0] = x[0];
    r_earth2sun_J2000[1] = x[1];
    r_earth2sun_J2000[2] = x[2];
    
    v_mag( &r, r_earth2sun_J2000);
    r3 = pow( r, 3.0 );
    coeff = PARAMS->SUN.GRAVITY.mu / r3;
    v_copy(pert_indirect, r_earth2sun_J2000);
    v_scale(pert_indirect, pert_indirect, coeff);

    v_sub(r_cg2sun_J2000, r_earth2sun_J2000, r_i2cg_INRTL);
    v_mag( &r, r_cg2sun_J2000);
    r3 = pow( r, 3.0 );
    coeff = PARAMS->SUN.GRAVITY.mu / r3;
    v_copy(pert_direct, r_cg2sun_J2000); // used to be: v_copy(pert_direct, r_earth2sun_J2000) but corrected based on Vallado3 page 35 equation (1-35)
    v_scale(pert_direct, pert_direct, coeff);
    v_sub(pert_total, pert_direct , pert_indirect);
    v_add(dvdt, dvdt , pert_total);
  }


  /* Solar radiation pressure */
  if (INTEGRATOR->include_solar_pressure == 1){
    compute_solar_pressure(a_solar_pressure_INRTL, r_i2cg_INRTL, v_i2cg_INRTL, et[0], PARAMS, INTEGRATOR, CONSTELLATION, et_initial_epoch, et_sc_initial, index_in_attitude_interpolated);
    v_add(dvdt, dvdt , a_solar_pressure_INRTL);
    //      v_norm_print(a_solar_pressure_INRTL, "\nSun");
  }

  /* Earth radiation pressure */
  if (INTEGRATOR->include_earth_pressure == 1){
    if (OPTIONS->opengl != 1){
    compute_earth_pressure(a_earth_pressure_INRTL, r_i2cg_INRTL, v_i2cg_INRTL, et[0], PARAMS, INTEGRATOR, CONSTELLATION, et_initial_epoch, et_sc_initial, index_in_attitude_interpolated);
    //           v_norm_print(a_earth_pressure_INRTL, "Earth");

    v_add(dvdt, dvdt , a_earth_pressure_INRTL);
    }
    else{
      printf("!***!\nFor now, SpOCK can't compute the Earth radiation pressure with a 3D model. It will be ignored.\n!***!\n");
    }
  }

  
  // Drag

  if (INTEGRATOR->include_drag == 1){

    compute_drag( adrag_i2cg_INRTL, r_i2cg_INRTL, v_i2cg_INRTL, et[0], PARAMS, INTEGRATOR, et_initial_epoch, et_sc_initial, density, index_in_attitude_interpolated, index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc, iDebugLevel, SC);
    v_add( dvdt, dvdt, adrag_i2cg_INRTL);
  }



  if (INTEGRATOR->thrust == 1){

    if ((et[0] >= OPTIONS->et_thrust_start) && (et[0] <= OPTIONS->et_thrust_stop)){
      compute_thrust(athrust_i2cg_INRTL, r_i2cg_INRTL, v_i2cg_INRTL, OPTIONS);
      //      etprint(et[0], "et");
      v_add( dvdt, dvdt, athrust_i2cg_INRTL);
    }
    
  }

  //  printf("\n");
  return 0;
}

/* /\* // STILL NEED TO CHECK THAT FUNCTION IS CORRECT *\/ */

/*  ///////////////////////////////////////////////////////////////////////////////////////// */
/* // */
/* //  Name:           load_gravity */
/* //  Purpose:        Loads EGM96 gravity coefficients */
/* //  Assumptions:    This isn't agnositic to file format, don't change current file */
/* //  References      BMW */
/* // */
/* //  Change Log: */
/* //      |   Developer   |       Date    |   SCR     |   Notes */
/* //      | --------------|---------------|-----------|------------------------------- */
/* //      | J. Getchius   | 05/17/2015    |   ---     | Initial Implementation */
/* //      | C. Bussy-Virat| 10/14/2015    |   ---     | Change the inputs */
/* // */
/* ///////////////////////////////////////////////////////////////////////////////////////// */
/* //newstructure */
/* //int load_gravity( GRAVITY_T *GRAVITY, char main_directory_location[256]) */
/* int load_gravity( GRAVITY_T *GRAVITY,  char path_to_spice[256]) */
/* //newstructure */
/* { */
/*   FILE *fp = NULL; */
/*   char *line = NULL; */
/*   size_t len = 0; */
/*   ssize_t read; */
/*   int l = 0; */
/*   int m = 0; */
/*   float C; */
/*   float S; */
/*   float Csig; */
/*   float Ssig; */
/*   char text_location[256]; */
/*       double fac, facnum, facden, kfac; */

/*   //  strcpy(text_location, main_directory_location); */
/*   //newstructure */
/*   //  egm96_to360_not_norm.txt is put with the SPICE files path_to_spice. ASSUMPTION: the path of spice specified in the main inoput file in the sectioN #SPICE must be the same as the path of SPICE installation in the Makefile, with /data at the end. So if in the Makefile the SPICE directory is hello/hi/ then the path of spice in the main input file in section #SPICE must be hello/hi/data/ */
/*   strcpy(text_location,path_to_spice); */
/*   strcat(text_location, "EGM2008_to2190_TideFree_stop85");//EGM2008_to2190_TideFree_stop85");//"egm96_to360.txt");//_not_norm.txt"); */
/*     //  strcpy(text_location, "input/egm96_to360_not_norm.txt"); */
/*   //newstructure */
/*   // !!!!!!!!!! to read normalized coefficients then put ./code/egm96_to360.txt */
/*   fp = fopen(text_location, "r"); */
/*   if (fp) { */
/*     while ( (read = getline(&line, &len, fp)) != -1 ) { */

/*       sscanf(line, "%i %i %e %e %e %e" , &l, &m, &C, &S, &Csig, &Ssig); */
/*       GRAVITY->Clm[l][m] = (double)C; // l is the degree (CBV) */
/*       GRAVITY->Slm[l][m] = (double)S; // m is the order (CBV) */


/*       if (m == 0){ */
/* 	kfac = 1; */
/*       } */
/*       else if (m > 0){ */
/* 	kfac = 2; */
/*       } */
/*       else */
/* 	{ */
/* 	  print_test(); */
/* 	  exitf(); */
/* 	} */

/*       facnum = kfac * (2*l+1)*(double)factorial(l-m); */
/*       facden = (double)factorial(l+m); */
/*       //printf("%d %d %e %e\n", l,m,facnum, facden); */
/*       fac = sqrt(facnum/facden); */
/* /\*            printf("%d %d %e %e %e %e\n" , l, m, GRAVITY->Clm[l][m], GRAVITY->Slm[l][m]); *\/ */
/* /\*            printf("%d %d %e %e\n" , l, m, (double)factorial(l-m), (double)factorial(l+m)); *\/ */

/*       GRAVITY->Clm[l][m]  = GRAVITY->Clm[l][m] * fac; */
/*       GRAVITY->Slm[l][m] = GRAVITY->Slm[l][m] * fac; */

/*     } */
    
/*     free(line); */
/*     printf("emg08\n"); */
/*     fclose(fp); */

/*   } else { */

/*     printf("The file input/egm96_to360_not_norm.txt has not been opened. The program will stop.\n"); */
/*     MPI_Finalize(); */
/*     exit(0); */

/*   } */

/*   return 0; */

/* } */


// FUNCTION IS CORRECT
/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           load_gravity
//  Purpose:        Loads EGM96 gravity coefficients
//  Assumptions:    This isn't agnositic to file format, don't change current file
//  References      BMW
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/17/2015    |   ---     | Initial Implementation
//      | C. Bussy-Virat| 10/14/2015    |   ---     | Change the inputs
//
/////////////////////////////////////////////////////////////////////////////////////////
//newstructure
//int load_gravity( GRAVITY_T *GRAVITY, char main_directory_location[256])
int load_gravity( GRAVITY_T *GRAVITY,  char path_to_spice[256])
//newstructure
{
  FILE *fp = NULL;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  int l = 0;
  int m = 0;
  float C;
  float S;
  float Csig;
  float Ssig;
  char text_location[256];


  //  strcpy(text_location, main_directory_location);
  //newstructure
  //  egm96_to360_not_norm.txt is put with the SPICE files path_to_spice. ASSUMPTION: the path of spice specified in the main inoput file in the sectioN #SPICE must be the same as the path of SPICE installation in the Makefile, with /data at the end. So if in the Makefile the SPICE directory is hello/hi/ then the path of spice in the main input file in section #SPICE must be hello/hi/data/
  strcpy(text_location,path_to_spice);
  strcat(text_location, "egm96_to360_not_norm.txt");
    //  strcpy(text_location, "input/egm96_to360_not_norm.txt");
  //newstructure
  // !!!!!!!!!! to read normalized coefficients then put ./code/egm96_to360.txt
  fp = fopen(text_location, "r");
  if (fp) {
    while ( (read = getline(&line, &len, fp)) != -1 ) {

      sscanf(line, "%i %i %e %e %e %e" , &l, &m, &C, &S, &Csig, &Ssig);
      GRAVITY->Clm[l][m] = (double)C; // l is the degree (CBV)
      GRAVITY->Slm[l][m] = (double)S; // m is the order (CBV)

    }
    
    free(line);
    
    fclose(fp);

  } else {

    printf("The file input/egm96_to360_not_norm.txt has not been opened. The program will stop.\n");
    MPI_Finalize();
    exit(0);

  }

  return 0;
    
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           compute_gravity
//  Purpose:        Computes gravity acceleration using EGM96
//  Assumptions:    Limited to 360x360
//  References      Vallado
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|---------------------------------------------
//      | J. Getchius   | 05/09/2015    |   ---     | Initial Implementatio
//      | C. Bussy-Virat| 07/19/2015    |   ---     | Corrections of term2 (from Vallado's theory (Vallado3 p. 548))
//
/////////////////////////////////////////////////////////////////////////////////////////
int compute_gravity(    double      a_i2cg_INRTL[3],
                        double      r_i2cg_INRTL[3],
                        double      et,
                        GRAVITY_T   *Gravity,
                        int         degree,
			char earth_fixed_frame[100],
			double earth_flattening, double earth_radius,
			SPACECRAFT_T *SC,
			int gravity_map,
			CONSTELLATION_T *CONSTELLATION)
			//                        int         order)
{

  // Declarations
  int iradius_arr[4], ilat_arr[4], ilon_arr[4];
  double *xinter_radius, *xinter_lon, *xinter_lat;
    double *yinter;
    int order_interpo_map = 4;
    xinter_radius = malloc(4 * sizeof(double));
    yinter = malloc(4 * sizeof(double));
    xinter_lon = malloc(4 * sizeof(double));
    xinter_lat = malloc(4 * sizeof(double));
    double y_radius1_lat1, y_radius1_lat2, y_radius2_lat1, y_radius2_lat2;
    int iradius1, ilat1, ilon1, iradius2, ilat2, ilon2, ilon3, ilat3, iradius3, ilon0, ilat0, iradius0;
        double long_gc_corr;
	   double xradius, xlat, xlon;
        int ii;
    double y_radius1_lat1_lon1, y_radius1_lat1_lon2, y_radius1_lat2_lon1, y_radius1_lat2_lon2, y_radius1;
        double y_radius2_lat1_lon1, y_radius2_lat1_lon2, y_radius2_lat2_lon1, y_radius2_lat2_lon2, y_radius2;

	double y_radius1_lat1_lon3, y_radius1_lat2_lon3;
        double y_radius1_lat3_lon1, y_radius1_lat3_lon2, y_radius1_lat3_lon3;
	    double y_radius2_lat1_lon3, y_radius2_lat2_lon3;

        double y_radius3_lat1_lon1, y_radius3_lat1_lon2, y_radius3_lat2_lon1, y_radius3_lat2_lon2, y_radius3;
	    double y_radius3_lat1_lon3, y_radius3_lat2_lon3;

	    double y_radius1_lat3;  double y_radius2_lat3_lon3; double y_radius2_lat3_lon1; double y_radius2_lat3_lon2; double y_radius2_lat3;

	    double y_radius3_lat3_lon3; double y_radius3_lat3_lon1; double y_radius3_lat3_lon2; double y_radius3_lat3;
	    	    double y_radius3_lat1; double y_radius3_lat2; 

	double y_radius1_lat1_lon0, y_radius1_lat2_lon0;
        double y_radius1_lat0_lon1, y_radius1_lat0_lon2, y_radius1_lat0_lon0;
	    double y_radius2_lat1_lon0, y_radius2_lat2_lon0;

        double y_radius0_lat1_lon1, y_radius0_lat1_lon2, y_radius0_lat2_lon1, y_radius0_lat2_lon2, y_radius0;
	    double y_radius0_lat1_lon0, y_radius0_lat2_lon0;

	    double y_radius1_lat0;  double y_radius2_lat0_lon0; double y_radius2_lat0_lon1; double y_radius2_lat0_lon2; double y_radius2_lat0;

	    double y_radius0_lat0_lon0; double y_radius0_lat0_lon1; double y_radius0_lat0_lon2; double y_radius0_lat3;

	    double y_radius1_lat0_lon3, y_radius1_lat3_lon0;
	    double y_radius0_lat0_lon3; double y_radius0_lat1_lon3; double y_radius0_lat2_lon3; double y_radius0_lat3_lon0; double y_radius0_lat3_lon1; double y_radius0_lat3_lon2; double y_radius0_lat3_lon3; double y_radius0_lat0; double y_radius0_lat1; double y_radius0_lat2;
	    double y_radius2_lat0_lon3; double y_radius2_lat3_lon0; double y_radius3_lat0_lon0; double y_radius3_lat0_lon1; double y_radius3_lat0_lon2; double y_radius3_lat0_lon3; double y_radius3_lat1_lon0; double y_radius3_lat2_lon0; double y_radius3_lat3_lon0; double y_radius3_lat0;
		    
	    double r_ecef2cg_ECEF[3];
  double a_ecef2cg_ECEF[3];
  double T_J2000_to_ECEF[3][3];
  double T_ECEF_to_J2000[3][3];
  double rmag;
  double rmag2;
  double rmag3;
  double lat_gc; // phi
  double slat_gc;
  double long_gc;
  double dUdr     = 0.0;
  double dUdlat   = 0.0;
  double dUdlong  = 0.0;
  double rad_over_r;
  double Plm;
  double Plm_plus1;
  double term;
  double term2;
  double coeff;
  double r_xy;
  int l;
  int m;
    
  // Begin algorithm
/*    // Update Geodetic State */
/*   double geodetic[3]; */
/* 	eci2lla(r_i2cg_INRTL , et, geodetic ); */
/* 	double altitude, latitude, longitude; */
/* 	altitude = geodetic[2]; */
/* 	latitude = geodetic[0]; */
/* 	longitude = geodetic[1];  */


/* 	// update the planet fixed state		 */
/* 	  geodetic_to_geocentric(earth_flattening,             */
/* 				 altitude, */
/* 				 latitude, */
/* 				 longitude, */
/* 				 earth_radius,        */
/* 				 r_ecef2cg_ECEF) ; */


  // Get Earth Fixed state
  pxform_c("J2000", earth_fixed_frame, et, T_J2000_to_ECEF);
  m_x_v(r_ecef2cg_ECEF, T_J2000_to_ECEF, r_i2cg_INRTL);

  // Compute geocentric latitude and longitude
  v_mag( &rmag, r_ecef2cg_ECEF);
  lat_gc  = asin( r_ecef2cg_ECEF[2] / rmag );
  long_gc = atan2(r_ecef2cg_ECEF[1], r_ecef2cg_ECEF[0]);
    
  // Some intermediate
  rad_over_r  = Gravity->radius / rmag;
  slat_gc     = sin( lat_gc );
  r_xy        = sqrt( r_ecef2cg_ECEF[0]*r_ecef2cg_ECEF[0] + r_ecef2cg_ECEF[1]*r_ecef2cg_ECEF[1]);
  rmag2       = pow( rmag, 2);
  rmag3       = pow( rmag, 3);

  //   Compute Partial of potential function wrt range, lat, long
  //  printf("degree = %d\n", degree);


  if (gravity_map != 1){ // if the user doesn't want to use the 3d gravity map option
  for (l = 2; l <= degree; l++) {

    for ( m = 0; m <= l; m++) { // !!!!!!!!!!!!! replace 0 with l

      Plm         = pow(-1,m) * gsl_sf_legendre_Plm( l, m, slat_gc ); // in ~/gsl-1.16/specfunc/legendre_poly.c // pow(-1,m) has been added to agree with Vallado's theory (the C library uses a pow(-1,m) that the theory does not))

      if ((m+1) > l) {
            
  	Plm_plus1 = 0.0;
            
      } else {
                
  	Plm_plus1   = pow(-1,m+1) * gsl_sf_legendre_Plm( l, (m+1), slat_gc ); // !!!!!!!!!!  pow(-1,m+1) has been added to agree with Vallado's theory (the C library uses a pow(-1,m+1) that the theory does not))
            
      }


      coeff       = pow(rad_over_r, l);
      //      printf("%d %d %e %e\n", l,m,Gravity->Clm[l][m], Gravity->Slm[l][m]);
      term        = (Gravity->Clm[l][m] * cos( m*long_gc ) + Gravity->Slm[l][m] * sin( m*long_gc ));
      term2       = (-Gravity->Clm[l][m] * sin( m*long_gc ) + Gravity->Slm[l][m] * cos( m*long_gc )); // Modif by CBV 07-19-2015 from Vallado3 p. 548
            
      dUdr    += coeff * (l + 1.0) * Plm * term;
      dUdlat  += coeff * ( Plm_plus1  - m * tan( lat_gc) * Plm ) * term;
      dUdlong += coeff * m * Plm * term2;
            
    }
  }


  //  exit(0);

  dUdr    = -dUdr     * Gravity->mu / rmag2;
  dUdlat  = dUdlat    * Gravity->mu / rmag;
  dUdlong = dUdlong   * Gravity->mu / rmag;

  } // end of if the user doesn't want to use the 3d gravity map option
  else{ // if the user wants to use the 3d gravity map option
    // determine the bin for the radius

    // determine the bin for the radius
    compute_iradius_gravity_map(iradius_arr, Gravity, rmag);
    iradius0 = iradius_arr[0]; iradius1 = iradius_arr[1]; iradius2 = iradius_arr[2]; iradius3 = iradius_arr[3];

    // determine the bin for the lat
    compute_ilat_gravity_map(ilat_arr, Gravity, lat_gc);
    ilat0 = ilat_arr[0]; ilat1 = ilat_arr[1]; ilat2 = ilat_arr[2]; ilat3 = ilat_arr[3];

    // determine the bin for the lon
    if (long_gc >= 0){ // long_gc_corr varies from 0 to 2*M_PI
      long_gc_corr = long_gc;
    }
    else{
      long_gc_corr = 2*M_PI + long_gc;
    }    
    compute_ilon_gravity_map(ilon_arr, Gravity, long_gc_corr);
    ilon0 = ilon_arr[0]; ilon1 = ilon_arr[1]; ilon2 = ilon_arr[2]; ilon3 = ilon_arr[3];

    // Determine the x bins (longitude, latitude, radius) for the interpolations
    gravity_map_xinter(xinter_lon, xinter_lat, xinter_radius,   long_gc_corr,  lat_gc,  rmag, Gravity, ilon_arr,  ilat_arr,  iradius_arr);

    // Interpolate over longitude, latitude, and radius      
    for (ii = 0; ii < 3; ii++){
      
      // RADIUS0
      y_radius0_lat0_lon0 = Gravity->gravity_map[iradius0][ilat0][ilon0][ii];
      y_radius0_lat0_lon1 = Gravity->gravity_map[iradius0][ilat0][ilon1][ii];
      y_radius0_lat0_lon2 = Gravity->gravity_map[iradius0][ilat0][ilon2][ii];
      y_radius0_lat0_lon3 = Gravity->gravity_map[iradius0][ilat0][ilon3][ii];
      y_radius0_lat1_lon0 = Gravity->gravity_map[iradius0][ilat1][ilon0][ii];
      y_radius0_lat1_lon1 = Gravity->gravity_map[iradius0][ilat1][ilon1][ii];
      y_radius0_lat1_lon2 = Gravity->gravity_map[iradius0][ilat1][ilon2][ii];
      y_radius0_lat1_lon3 = Gravity->gravity_map[iradius0][ilat1][ilon3][ii];
      y_radius0_lat2_lon0 = Gravity->gravity_map[iradius0][ilat2][ilon0][ii];
      y_radius0_lat2_lon1 = Gravity->gravity_map[iradius0][ilat2][ilon1][ii];
      y_radius0_lat2_lon2 = Gravity->gravity_map[iradius0][ilat2][ilon2][ii];
      y_radius0_lat2_lon3 = Gravity->gravity_map[iradius0][ilat2][ilon3][ii];
      y_radius0_lat3_lon0 = Gravity->gravity_map[iradius0][ilat3][ilon0][ii];
      y_radius0_lat3_lon1 = Gravity->gravity_map[iradius0][ilat3][ilon1][ii];
      y_radius0_lat3_lon2 = Gravity->gravity_map[iradius0][ilat3][ilon2][ii];
      y_radius0_lat3_lon3 = Gravity->gravity_map[iradius0][ilat3][ilon3][ii];

      // longitude
      order_interpo_map = 4;
      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius0_lat0_lon0, y_radius0_lat0_lon1, y_radius0_lat0_lon2, y_radius0_lat0_lon3);
      polynomial_interpo(&y_radius0_lat0, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius0_lat1_lon0, y_radius0_lat1_lon1, y_radius0_lat1_lon2, y_radius0_lat1_lon3);
      polynomial_interpo(&y_radius0_lat1, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius0_lat2_lon0, y_radius0_lat2_lon1, y_radius0_lat2_lon2, y_radius0_lat2_lon3);
      polynomial_interpo(&y_radius0_lat2, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius0_lat3_lon0, y_radius0_lat3_lon1, y_radius0_lat3_lon2, y_radius0_lat3_lon3);      
      polynomial_interpo(&y_radius0_lat3, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);      

      // latitude
      gravity_map_yinter_lat(yinter, &order_interpo_map, lat_gc, Gravity, y_radius0_lat0, y_radius0_lat1, y_radius0_lat2, y_radius0_lat3);
      polynomial_interpo(&y_radius0, order_interpo_map, xinter_lat, yinter, lat_gc*180/M_PI);  

      
      // RADIUS1
      y_radius1_lat0_lon0 = Gravity->gravity_map[iradius1][ilat0][ilon0][ii];
      y_radius1_lat0_lon1 = Gravity->gravity_map[iradius1][ilat0][ilon1][ii];
      y_radius1_lat0_lon2 = Gravity->gravity_map[iradius1][ilat0][ilon2][ii];
      y_radius1_lat0_lon3 = Gravity->gravity_map[iradius1][ilat0][ilon3][ii];
      y_radius1_lat1_lon0 = Gravity->gravity_map[iradius1][ilat1][ilon0][ii];
      y_radius1_lat1_lon1 = Gravity->gravity_map[iradius1][ilat1][ilon1][ii];
      y_radius1_lat1_lon2 = Gravity->gravity_map[iradius1][ilat1][ilon2][ii];
      y_radius1_lat1_lon3 = Gravity->gravity_map[iradius1][ilat1][ilon3][ii];
      y_radius1_lat2_lon0 = Gravity->gravity_map[iradius1][ilat2][ilon0][ii];
      y_radius1_lat2_lon1 = Gravity->gravity_map[iradius1][ilat2][ilon1][ii];
      y_radius1_lat2_lon2 = Gravity->gravity_map[iradius1][ilat2][ilon2][ii];
      y_radius1_lat2_lon3 = Gravity->gravity_map[iradius1][ilat2][ilon3][ii];
      y_radius1_lat3_lon0 = Gravity->gravity_map[iradius1][ilat3][ilon0][ii];
      y_radius1_lat3_lon1 = Gravity->gravity_map[iradius1][ilat3][ilon1][ii];
      y_radius1_lat3_lon2 = Gravity->gravity_map[iradius1][ilat3][ilon2][ii];
      y_radius1_lat3_lon3 = Gravity->gravity_map[iradius1][ilat3][ilon3][ii];

      // longitude
      order_interpo_map = 4;
      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius1_lat0_lon0, y_radius1_lat0_lon1, y_radius1_lat0_lon2, y_radius1_lat0_lon3);
      polynomial_interpo(&y_radius1_lat0, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius1_lat1_lon0, y_radius1_lat1_lon1, y_radius1_lat1_lon2, y_radius1_lat1_lon3);
      polynomial_interpo(&y_radius1_lat1, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius1_lat2_lon0, y_radius1_lat2_lon1, y_radius1_lat2_lon2, y_radius1_lat2_lon3);
      polynomial_interpo(&y_radius1_lat2, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius1_lat3_lon0, y_radius1_lat3_lon1, y_radius1_lat3_lon2, y_radius1_lat3_lon3);      
      polynomial_interpo(&y_radius1_lat3, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);      

      // Latitude
      gravity_map_yinter_lat(yinter, &order_interpo_map, lat_gc, Gravity, y_radius1_lat0, y_radius1_lat1, y_radius1_lat2, y_radius1_lat3);
      polynomial_interpo(&y_radius1, order_interpo_map, xinter_lat, yinter, lat_gc*180/M_PI);  
      	  /* printf("%e\n", y_radius1              ); */
	  /* 	  exitf(); */

      // RADIUS2
      y_radius2_lat0_lon0 = Gravity->gravity_map[iradius2][ilat0][ilon0][ii];
      y_radius2_lat0_lon1 = Gravity->gravity_map[iradius2][ilat0][ilon1][ii];
      y_radius2_lat0_lon2 = Gravity->gravity_map[iradius2][ilat0][ilon2][ii];
      y_radius2_lat0_lon3 = Gravity->gravity_map[iradius2][ilat0][ilon3][ii];
      y_radius2_lat1_lon0 = Gravity->gravity_map[iradius2][ilat1][ilon0][ii];
      y_radius2_lat1_lon1 = Gravity->gravity_map[iradius2][ilat1][ilon1][ii];
      y_radius2_lat1_lon2 = Gravity->gravity_map[iradius2][ilat1][ilon2][ii];
      y_radius2_lat1_lon3 = Gravity->gravity_map[iradius2][ilat1][ilon3][ii];
      y_radius2_lat2_lon0 = Gravity->gravity_map[iradius2][ilat2][ilon0][ii];
      y_radius2_lat2_lon1 = Gravity->gravity_map[iradius2][ilat2][ilon1][ii];
      y_radius2_lat2_lon2 = Gravity->gravity_map[iradius2][ilat2][ilon2][ii];
      y_radius2_lat2_lon3 = Gravity->gravity_map[iradius2][ilat2][ilon3][ii];
      y_radius2_lat3_lon0 = Gravity->gravity_map[iradius2][ilat3][ilon0][ii];
      y_radius2_lat3_lon1 = Gravity->gravity_map[iradius2][ilat3][ilon1][ii];
      y_radius2_lat3_lon2 = Gravity->gravity_map[iradius2][ilat3][ilon2][ii];
      y_radius2_lat3_lon3 = Gravity->gravity_map[iradius2][ilat3][ilon3][ii];
    //y_radius2_lat1 = xlon*y_radius2_lat1_lon2 + (1-xlon)*y_radius2_lat1_lon1;
    //        y_radius2_lat2 = xlon*y_radius2_lat2_lon2 + (1-xlon)*y_radius2_lat2_lon1;


      // longitude
      order_interpo_map = 4;
      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius2_lat0_lon0, y_radius2_lat0_lon1, y_radius2_lat0_lon2, y_radius2_lat0_lon3);
      polynomial_interpo(&y_radius2_lat0, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius2_lat1_lon0, y_radius2_lat1_lon1, y_radius2_lat1_lon2, y_radius2_lat1_lon3);
      polynomial_interpo(&y_radius2_lat1, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius2_lat2_lon0, y_radius2_lat2_lon1, y_radius2_lat2_lon2, y_radius2_lat2_lon3);
      polynomial_interpo(&y_radius2_lat2, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius2_lat3_lon0, y_radius2_lat3_lon1, y_radius2_lat3_lon2, y_radius2_lat3_lon3);      
      polynomial_interpo(&y_radius2_lat3, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);      

      // latitude
      gravity_map_yinter_lat(yinter, &order_interpo_map, lat_gc, Gravity, y_radius2_lat0, y_radius2_lat1, y_radius2_lat2, y_radius2_lat3);
      polynomial_interpo(&y_radius2, order_interpo_map, xinter_lat, yinter, lat_gc*180/M_PI);  
      	  /* printf("%e\n", y_radius2); */
	  /* 	  exitf(); */


      // RADIUS3
      y_radius3_lat0_lon0 = Gravity->gravity_map[iradius3][ilat0][ilon0][ii];
      y_radius3_lat0_lon1 = Gravity->gravity_map[iradius3][ilat0][ilon1][ii];
      y_radius3_lat0_lon2 = Gravity->gravity_map[iradius3][ilat0][ilon2][ii];
      y_radius3_lat0_lon3 = Gravity->gravity_map[iradius3][ilat0][ilon3][ii];
      y_radius3_lat1_lon0 = Gravity->gravity_map[iradius3][ilat1][ilon0][ii];
      y_radius3_lat1_lon1 = Gravity->gravity_map[iradius3][ilat1][ilon1][ii];
      y_radius3_lat1_lon2 = Gravity->gravity_map[iradius3][ilat1][ilon2][ii];
      y_radius3_lat1_lon3 = Gravity->gravity_map[iradius3][ilat1][ilon3][ii];
      y_radius3_lat2_lon0 = Gravity->gravity_map[iradius3][ilat2][ilon0][ii];
      y_radius3_lat2_lon1 = Gravity->gravity_map[iradius3][ilat2][ilon1][ii];
      y_radius3_lat2_lon2 = Gravity->gravity_map[iradius3][ilat2][ilon2][ii];
      y_radius3_lat2_lon3 = Gravity->gravity_map[iradius3][ilat2][ilon3][ii];
      y_radius3_lat3_lon0 = Gravity->gravity_map[iradius3][ilat3][ilon0][ii];
      y_radius3_lat3_lon1 = Gravity->gravity_map[iradius3][ilat3][ilon1][ii];
      y_radius3_lat3_lon2 = Gravity->gravity_map[iradius3][ilat3][ilon2][ii];
      y_radius3_lat3_lon3 = Gravity->gravity_map[iradius3][ilat3][ilon3][ii];
    //y_radius3_lat1 = xlon*y_radius3_lat1_lon2 + (1-xlon)*y_radius3_lat1_lon1;
    //        y_radius3_lat2 = xlon*y_radius3_lat2_lon2 + (1-xlon)*y_radius3_lat2_lon1;

      // longitude
      order_interpo_map = 4;
      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius3_lat0_lon0, y_radius3_lat0_lon1, y_radius3_lat0_lon2, y_radius3_lat0_lon3);
      polynomial_interpo(&y_radius3_lat0, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius3_lat1_lon0, y_radius3_lat1_lon1, y_radius3_lat1_lon2, y_radius3_lat1_lon3);
      polynomial_interpo(&y_radius3_lat1, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius3_lat2_lon0, y_radius3_lat2_lon1, y_radius3_lat2_lon2, y_radius3_lat2_lon3);
      polynomial_interpo(&y_radius3_lat2, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);

      gravity_map_yinter_lon(yinter, &order_interpo_map, long_gc_corr, Gravity, y_radius3_lat3_lon0, y_radius3_lat3_lon1, y_radius3_lat3_lon2, y_radius3_lat3_lon3);      
      polynomial_interpo(&y_radius3_lat3, order_interpo_map, xinter_lon, yinter, long_gc_corr*180/M_PI);      

      // latitude
      gravity_map_yinter_lat(yinter, &order_interpo_map, lat_gc, Gravity, y_radius3_lat0, y_radius3_lat1, y_radius3_lat2, y_radius3_lat3);
      polynomial_interpo(&y_radius3, order_interpo_map, xinter_lat, yinter, lat_gc*180/M_PI);  
      // Interpo over radius
      gravity_map_yinter_radius( yinter, &order_interpo_map,  rmag, Gravity,  y_radius0,  y_radius1,  y_radius2,  y_radius3);
	// Interpolate dUdr, Dudlat, DUdlong
    if (ii == 0){
      // dUdr = y_radius2*xradius + y_radius1*(1-xradius);
      polynomial_interpo(&dUdr, order_interpo_map, xinter_radius, yinter, rmag);
    }
    else if (ii == 1){
      //      dUdlat = y_radius2*xradius + y_radius1*(1-xradius);
      polynomial_interpo(&dUdlat, order_interpo_map, xinter_radius, yinter, rmag);    
    }
    else{
      //      dUdlong = y_radius2*xradius + y_radius1*(1-xradius);

      polynomial_interpo(&dUdlong, order_interpo_map, xinter_radius, yinter, rmag);    
    }


    }
      //      printf("%f %d, %f %d, %f %d\n", rmag,iradius, lat_gc*180/M_PI, ilat,long_gc_corr*180/M_PI,ilon);

      //      exitf();
  }

  //  printf("%f (%f %d) %f (%d) %f (%d) -> %e %e %e\n", rmag, rmag - Gravity->radius, iradius1, lat_gc*180/M_PI, ilat1, long_gc*180/M_PI, ilon1, dUdr, dUdlat, dUdlong);

  
  //  // Compute the Earth fixed accels
  term                = (1/rmag) * dUdr - (r_ecef2cg_ECEF[2] / (rmag*rmag * r_xy)) * dUdlat;
  term2               = (1/(r_xy * r_xy)) * dUdlong;
  a_ecef2cg_ECEF[0]   = term * r_ecef2cg_ECEF[0] - term2 * r_ecef2cg_ECEF[1] - Gravity->mu * r_ecef2cg_ECEF[0] / (rmag3);
  a_ecef2cg_ECEF[1]   = term * r_ecef2cg_ECEF[1] + term2 * r_ecef2cg_ECEF[0] - Gravity->mu * r_ecef2cg_ECEF[1] / (rmag3);
  a_ecef2cg_ECEF[2]   = (1/rmag) * dUdr * r_ecef2cg_ECEF[2] + r_xy / ( rmag2) * dUdlat - Gravity->mu  * r_ecef2cg_ECEF[2] / (rmag3);



  // Rotate to Inertial accels
  pxform_c(earth_fixed_frame, "J2000", et, T_ECEF_to_J2000);
  m_x_v( a_i2cg_INRTL, T_ECEF_to_J2000, a_ecef2cg_ECEF);
  v_copy(SC->a_i2cg_INRTL_gravity, a_i2cg_INRTL);
  //    v_print(a_ecef2cg_ECEF, "a_ecef2cg_ECEF"); 


  /* // J2 only  */
  /* double R0 = 6378.1370; */
  /* double J2 = 1082.65e-6; */

  /* double j2pert = J2*(3.0/2.0)*(R0/rmag)*(R0/rmag); */
  /* double j2sub = 5.0*(r_i2cg_INRTL[2]*r_i2cg_INRTL[2]/R0/R0); */

  /* a_i2cg_INRTL[0]   = -Gravity->mu*r_i2cg_INRTL[0]/(rmag*rmag*rmag) * (1.0 - j2pert*(j2sub-1)); */
  /* a_i2cg_INRTL[1]   = -Gravity->mu*r_i2cg_INRTL[1]/(rmag*rmag*rmag) * (1.0 - j2pert*(j2sub-1)); */
  /* a_i2cg_INRTL[2]   = -Gravity->mu*r_i2cg_INRTL[2]/(rmag*rmag*rmag) * (1.0 - j2pert*(j2sub-3)); */
    

  return 0;
}






/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           compute_power
//  Purpose:        Computes power from the solar arrays
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|---------------------------------------------------------
//      | C. Bussy-Virat| 08/02/2015    |   ---     | Initial implementation
//
/////////////////////////////////////////////////////////////////////////////////////////

int compute_power(INTEGRATOR_T    *INTEGRATOR,
		  CONSTELLATION_T *CONSTELLATION,
		  double          r_i2cg_INRTL[3],
		  double          v_i2cg_INRTL[3],
		  double          *et,
		  PARAMS_T        *PARAMS,
		  double          et_initial_epoch,
		  double          et_sc_initial,
		  int             index_in_attitude_interpolated)	
		
{

  // Declarations

  //  int index_in_attitude_interpolated;
  double cos_sun_elevation;
  //  double angle_normal_to_sun;
  double z_body[3];
  double solar_flux = 1358.;//1367.0; // solar flux is in W/m^2 (SI), which is in kg/s^3. So the fact that we express distances in km and not in m does not change the value of the solar_flux here.
  double x[6];
  double lt;
  double r_cg2sun_J2000[3];
  double r_cg2sun_J2000_normalized[3];
  double r_cg2sun_LVLH_normalized[3];
  double r_cg2sun_SC_normalized[3];
  double r_earth2sun_J2000[3];
  double T_inrtl_2_lvlh[3][3];
  double T_sc_to_lvlh[3][3];
  double T_lvlh_to_sc[3][3];
  double cos_phi; // phi is the angle between the satellite-Sun direction and the normal to the surface
  int sss;
  char shadow[256];
  double v_angle[3];
  int order_rotation[3];

  shadow_light( shadow, r_i2cg_INRTL, et[0], PARAMS);
  if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") == 0){  // if we run ensembles on the attitude. This block was moved out of the condition "if ( strcmp(shadow, "light") == 0)" to compute the (pitch, roll, yaw) angles even in the shadow, so that we can write these in an output file to compare with STK. If you don't want to compare with STK, then you can put it in the the condition "if ( strcmp(shadow, "light") == 0)"
    v_angle[0] = INTEGRATOR->attitude.pitch_angular_velocity_ensemble * ( et[0] - et_sc_initial ) + INTEGRATOR->attitude.pitch_ini_ensemble;
    v_angle[1] = INTEGRATOR->attitude.roll_angular_velocity_ensemble * ( et[0] - et_sc_initial ) + INTEGRATOR->attitude.roll_ini_ensemble;
    v_angle[2] = INTEGRATOR->attitude.yaw_angular_velocity_ensemble * ( et[0] - et_sc_initial ) + INTEGRATOR->attitude.yaw_ini_ensemble;
    order_rotation[0] = 1; // !!!!!!!! we might want to change that in the future
    order_rotation[1] = 2;
    order_rotation[2] = 3;
    INTEGRATOR->attitude.pitch_current = v_angle[0];
    INTEGRATOR->attitude.roll_current = v_angle[1];
    INTEGRATOR->attitude.yaw_current = v_angle[2];
        INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


  }

	  if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") == 0) {
	      v_angle[0] =  INTEGRATOR->attitude.pitch_for_attitude_ensemble  +  INTEGRATOR->attitude.pitch_angular_velocity_constant * ( et[0] - et_sc_initial );
	      v_angle[1] =  INTEGRATOR->attitude.roll_for_attitude_ensemble  +  INTEGRATOR->attitude.roll_angular_velocity_constant * ( et[0] - et_sc_initial );
	      v_angle[2] =  INTEGRATOR->attitude.yaw_for_attitude_ensemble  +  INTEGRATOR->attitude.yaw_angular_velocity_constant * ( et[0] - et_sc_initial );

	      order_rotation[0]  = 1; order_rotation[1]  = 2; order_rotation[2]  = 3;
	   
	      INTEGRATOR->attitude.pitch_current = v_angle[0];
	      INTEGRATOR->attitude.roll_current = v_angle[1];
	      INTEGRATOR->attitude.yaw_current = v_angle[2];
	          INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


	      
	  }

  //  char cu_time[256];
	  if ( (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "sun_pointed") != 0) ){ // otherwise (atittude is nadir, sun_pointed or manual (from an input file))  
    //    et2utc_c(et[0], "ISOC" ,0 ,255 , cu_time);
    //    index_in_attitude_interpolated = floor( ( et[0] - et_initial_epoch ) / ( INTEGRATOR->dt / 2.0) ) ; 
	    if (INTEGRATOR->file_is_quaternion == 0){
    v_angle[0] = INTEGRATOR->attitude.pitch[index_in_attitude_interpolated];
    v_angle[1] = INTEGRATOR->attitude.roll[index_in_attitude_interpolated];
    v_angle[2] = INTEGRATOR->attitude.yaw[index_in_attitude_interpolated];
    order_rotation[0] = INTEGRATOR->attitude.order_pitch[index_in_attitude_interpolated];
    order_rotation[1] = INTEGRATOR->attitude.order_roll[index_in_attitude_interpolated];
    order_rotation[2] = INTEGRATOR->attitude.order_yaw[index_in_attitude_interpolated];
    INTEGRATOR->attitude.pitch_current = v_angle[0];
    INTEGRATOR->attitude.roll_current = v_angle[1];
    INTEGRATOR->attitude.yaw_current = v_angle[2];
        INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


	    }
	    else{
	      q_copy( INTEGRATOR->attitude.quaternion_current, INTEGRATOR->attitude.quaternion[index_in_attitude_interpolated]);
	    }


  }



  if ( strcmp(shadow, "light") == 0) {

    spkez_c(10, et[0], "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.

    r_earth2sun_J2000[0] = x[0];
    r_earth2sun_J2000[1] = x[1];
    r_earth2sun_J2000[2] = x[2];
     
    v_sub(r_cg2sun_J2000, r_earth2sun_J2000, r_i2cg_INRTL);
    v_norm(r_cg2sun_J2000_normalized, r_cg2sun_J2000);


      /* r_cg2sun_J2000_normalized inertial to LVLH */
      compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, r_i2cg_INRTL, v_i2cg_INRTL);
      m_x_v(r_cg2sun_LVLH_normalized, T_inrtl_2_lvlh, r_cg2sun_J2000_normalized);
      /* r_cg2sun_J2000_normalized LVLH to body */
      compute_T_sc_to_lvlh( T_sc_to_lvlh, v_angle, order_rotation, INTEGRATOR->attitude.attitude_profile, et,  r_i2cg_INRTL, v_i2cg_INRTL, INTEGRATOR->file_is_quaternion, INTEGRATOR->attitude.quaternion_current, PARAMS);
      m_trans(T_lvlh_to_sc, T_sc_to_lvlh);
      m_x_v(r_cg2sun_SC_normalized, T_lvlh_to_sc, r_cg2sun_LVLH_normalized ); 


    z_body[0] = 0; z_body[1] = 0; z_body[2] = 1;
    v_dot(&cos_sun_elevation, r_cg2sun_SC_normalized, z_body);
    INTEGRATOR->sun_elevation =  M_PI/2 - acos(cos_sun_elevation)  ;



  // openGL computation of cross-section area
  if (INTEGRATOR->opengl_power == 1){ 
    // // compute the azimuth and elevation of the sun vector in body frame. 
    // // Actually the spherical coordinates are used. 
    // // theta is the angle between z_body and the sun vector (from 0 to -180 deg). 
    // // phi is the angle between x_body and the projection of the sun vector on the (x_body, y_body) plane (from 0 to 360 deg).
    double x_body[3], y_body[3], z_body[3];
    x_body[0] = 1; x_body[1] = 0; x_body[2] = 0;
    y_body[0] = 0; y_body[1] = 1; y_body[2] = 0;
    double theta, phi;
    // // // theta
    theta = acos(r_cg2sun_SC_normalized[2]); // r_cg2sun_SC_normalized[2] is euqal to the projection of r_cg2sun_SC_normalized on z_body
    // // // phi
    phi = acos(r_cg2sun_SC_normalized[0]);
    if (r_cg2sun_SC_normalized[1] < 0){
      phi = 2*M_PI - phi;
    }


    double theta_deg, phi_deg;
    theta_deg = theta * RAD2DEG;
    phi_deg = phi * RAD2DEG;
    // // get the cross-section area from the file generated by openGL
    // // // compute the indices itheta and iphi
    int itheta, iphi, iface;
    double A_ref = 0;
    itheta = (int)((theta_deg - CONSTELLATION->area_attitude_opengl_theta0) / CONSTELLATION->area_attitude_opengl_dtheta);
    iphi = (int)((phi_deg - CONSTELLATION->area_attitude_opengl_phi0) / CONSTELLATION->area_attitude_opengl_dphi);

    A_ref = CONSTELLATION->area_solar_panel_attitude_opengl[itheta][iphi];
	
	
/*       printf("%f %f\n", theta_deg, phi_deg); */
/*       printf("theta %.2f, phi %.2f\n", itheta*CONSTELLATION->area_attitude_opengl_dtheta + CONSTELLATION->area_attitude_opengl_theta0, iphi*CONSTELLATION->area_attitude_opengl_dphi + CONSTELLATION->area_attitude_opengl_phi0); */
/*       printf("A_ref %f\n" , A_ref*1e10); */
      // we assume all surfaces have the same solar cell efficiency. In generate_ephemerides.c the total power is computed at the sum of the power of each surface. This was done when the 3d opengl modeling version was not yet a thing. With opengl, the total power is output in opengl_filename_solar_power and we don't know the power of each surface of the satellite from the geometry file (first name of 5th line in section #SPACECRAFT of main input file). So we assign this toal power to surface 0 snd leave all other surfaces to 0. This can be paradoxical with the fact that surafec 0 is not necessarily assigned a solar panel. But it works
      for (sss = 0; sss < INTEGRATOR->nb_surfaces; sss++){
	INTEGRATOR->surface[sss].power_per_surface = 0;
      }
      INTEGRATOR->surface[0].power_per_surface = A_ref * solar_flux * INTEGRATOR->solar_cell_efficiency * 1e6; // "1e6" is here to convert the area of the solar panel from km^2 to m^2
      //      printf("power %f\n", INTEGRATOR->surface[0].power_per_surface);
      //      printf("A_ref = %f\n", A_ref*1e10);
  }


  else{
    for (sss = 0; sss < INTEGRATOR->nb_surfaces; sss++){
      INTEGRATOR->surface[sss].power_per_surface = 0.0;
      if (INTEGRATOR->surface[sss].area_solar_panel > 0){
      v_dot(&cos_phi, r_cg2sun_SC_normalized, INTEGRATOR->surface[sss].normal);
      // // BLOCK BELOW USELESS, SO ERASE IT (OR COMMENT IT)

      /* angle_normal_to_sun = acos(cos_phi) * 180 / M_PI; */
      
      /* //      printf("angle = %f\n", angle_normal_to_sun); */
      // // END OF BLOCK BELOW USELESS, SO ERASE IT (OR COMMENT IT)
      //      printf("%d: %f\n", sss,cos_phi);
      if (cos_phi > 0){
	
	//printf("surface %d | A_ref = %f\n", sss, INTEGRATOR->surface[sss].area_solar_panel * cos_phi * 1e10);
	INTEGRATOR->surface[sss].power_per_surface = INTEGRATOR->surface[sss].area_solar_panel * solar_flux * cos_phi * INTEGRATOR->solar_cell_efficiency * 1e6; // "1e6" is here to convert the area of the solar panel from km^2 to m^2
	//      printf("power for surface %d: %f\n", sss, INTEGRATOR->surface[sss].power_per_surface);
      }
      else{
	INTEGRATOR->surface[sss].power_per_surface = 0.0;
      }
      }
      else{
	INTEGRATOR->surface[sss].power_per_surface = 0.0;
      }
    }
  }

  }
  else{
    INTEGRATOR->sun_elevation = -999*DEG2RAD;
    for (sss = 0; sss < INTEGRATOR->nb_surfaces; sss++){
      INTEGRATOR->surface[sss].power_per_surface = 0.0;   
    } 
  }
  //  exit(0);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           compute_solar_pressure
//  Purpose:        Computes solar pressure acceleration
//  Assumptions:    VALLADO OR STK THEORY CAN BE USED, CHOOSE THE ONE YOU WANT (SEE COMMENTS IN THE CODE)
//  References      Vallado version 3 (section 8.6.4) AND STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm)
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|---------------------------------------------------------
//      | C. Bussy-Virat| 08/02/2015    |   ---     | Initial implementation
//
/////////////////////////////////////////////////////////////////////////////////////////

int compute_solar_pressure(double          a_solar_pressure_INRTL[3],
			   double          r_i2cg_INRTL[3],
			   double          v_i2cg_INRTL[3],
			   double          et,
			   PARAMS_T        *PARAMS,
			   INTEGRATOR_T    *INTEGRATOR,
			   CONSTELLATION_T *CONSTELLATION,
			   double          et_initial_epoch,
			   double          et_sc_initial,
			   int             index_in_attitude_interpolated)
		
{


  // Declarations
	  double k = 1; // in STK 'fraction of the solar disk visible at satellite location', set to 1 here
	  double solar_luminosity = 3.823e26; // in Watts

    double dist_sat_to_sun;
    //  double solar_flux = 1358.0; // solar flux is in W/m^2 (SI), which is in kg/s^3. So the fact that we express distances in km and not in m does not change the value of the solar_flux here.
  double a_solar_pressure_in_body[3]; // solar pressure acceleration in the SC reference system
  double x[6];
  double lt;
  double r_cg2sun_J2000[3];
  double r_cg2sun_J2000_normalized[3];
  double r_cg2sun_LVLH_normalized[3];
  double r_cg2sun_SC_normalized[3];
  double r_earth2sun_J2000[3];
  double T_inrtl_2_lvlh[3][3];
  double T_lvlh_2_inrtl[3][3];
  double T_sc_to_lvlh[3][3];
  double T_lvlh_to_sc[3][3];
  double cos_phi; // phi is the angle between the satellite-Sun direction and the normal to the surface
  double light_speed = 299792.458; // spped of light in km/s
  double a_solar_pressure_in_LVLH[3];
  //  double term1, term2;
  int sss;
  char shadow[256];
  //  int index_in_attitude_interpolated;
  double v_angle[3];
  int order_rotation[3];

  shadow_light( shadow, r_i2cg_INRTL, et, PARAMS);

  if (INTEGRATOR->coll_vcm != 1){
  if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") == 0){
    v_angle[0] = INTEGRATOR->attitude.pitch_angular_velocity_ensemble * ( et - et_sc_initial ) + INTEGRATOR->attitude.pitch_ini_ensemble;
    v_angle[1] = INTEGRATOR->attitude.roll_angular_velocity_ensemble * ( et - et_sc_initial ) + INTEGRATOR->attitude.roll_ini_ensemble;
    v_angle[2] = INTEGRATOR->attitude.yaw_angular_velocity_ensemble * ( et - et_sc_initial ) + INTEGRATOR->attitude.yaw_ini_ensemble;
    order_rotation[0] = 1; // !!!!!!!! we might want to change that in the future                                                                     
    order_rotation[1] = 2;
    order_rotation[2] = 3;
    INTEGRATOR->attitude.pitch_current = v_angle[0];
    INTEGRATOR->attitude.roll_current = v_angle[1];
    INTEGRATOR->attitude.yaw_current = v_angle[2];

  }

	  if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") == 0) {
	      v_angle[0] =  INTEGRATOR->attitude.pitch_for_attitude_ensemble  +  INTEGRATOR->attitude.pitch_angular_velocity_constant * ( et - et_sc_initial );
	      v_angle[1] =  INTEGRATOR->attitude.roll_for_attitude_ensemble  +  INTEGRATOR->attitude.roll_angular_velocity_constant * ( et - et_sc_initial );
	      v_angle[2] =  INTEGRATOR->attitude.yaw_for_attitude_ensemble  +  INTEGRATOR->attitude.yaw_angular_velocity_constant * ( et - et_sc_initial );

	      order_rotation[0]  = 1; order_rotation[1]  = 2; order_rotation[2]  = 3;
	   
	      INTEGRATOR->attitude.pitch_current = v_angle[0];
	      INTEGRATOR->attitude.roll_current = v_angle[1];
	      INTEGRATOR->attitude.yaw_current = v_angle[2];
	          INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


	      
	  }


	  if ( (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "sun_pointed") != 0) ){ // otherwise (atittude is nadir, sun_pointed or manual (from an input file))    
    //    index_in_attitude_interpolated = floor( ( et - et_sc_initial ) / ( INTEGRATOR->dt / 2.0) ) ; 
	    if (INTEGRATOR->file_is_quaternion == 0){
    v_angle[0] = INTEGRATOR->attitude.pitch[index_in_attitude_interpolated];
    v_angle[1] = INTEGRATOR->attitude.roll[index_in_attitude_interpolated];
    v_angle[2] = INTEGRATOR->attitude.yaw[index_in_attitude_interpolated];
    order_rotation[0] = INTEGRATOR->attitude.order_pitch[index_in_attitude_interpolated];
    order_rotation[1] = INTEGRATOR->attitude.order_roll[index_in_attitude_interpolated];
    order_rotation[2] = INTEGRATOR->attitude.order_yaw[index_in_attitude_interpolated];
    INTEGRATOR->attitude.pitch_current = v_angle[0];
    INTEGRATOR->attitude.roll_current = v_angle[1];
    INTEGRATOR->attitude.yaw_current = v_angle[2];
        INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


	    }
	    else{
	      q_copy( INTEGRATOR->attitude.quaternion_current, INTEGRATOR->attitude.quaternion[index_in_attitude_interpolated]);
	    }
  }
  } // end of no collision with VCM as colllision input file

    // !!! to comment beta angle 
    spkez_c(10, et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.

    r_earth2sun_J2000[0] = x[0];
    r_earth2sun_J2000[1] = x[1];
    r_earth2sun_J2000[2] = x[2];

    double r_earth2sun_J2000_norm[3];
    v_norm(r_earth2sun_J2000_norm, r_earth2sun_J2000);
    double r_cross_v[3];
    v_cross(r_cross_v,r_i2cg_INRTL,v_i2cg_INRTL);
    double r_cross_v_norm[3];
    v_norm(r_cross_v_norm, r_cross_v);
    double r_cross_v_norm_dot_r_earth2sun_J2000_norm;
    v_dot(&r_cross_v_norm_dot_r_earth2sun_J2000_norm, r_cross_v_norm, r_earth2sun_J2000_norm);
    double angle_h_to_earth_sun;
    angle_h_to_earth_sun = acos(r_cross_v_norm_dot_r_earth2sun_J2000_norm);
    
    INTEGRATOR->beta_angle = M_PI/2. - angle_h_to_earth_sun;
    //printf("beta: %f\n", INTEGRATOR->beta_angle*180./M_PI);
    // !!! end of to comment beta angle

  if ( strcmp(shadow, "light") == 0) {


    spkez_c(10, et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.

    r_earth2sun_J2000[0] = x[0];
    r_earth2sun_J2000[1] = x[1];
    r_earth2sun_J2000[2] = x[2];


     
    v_sub(r_cg2sun_J2000, r_earth2sun_J2000, r_i2cg_INRTL);
    v_mag( &dist_sat_to_sun, r_cg2sun_J2000 );
    v_norm(r_cg2sun_J2000_normalized, r_cg2sun_J2000);


  if (INTEGRATOR->coll_vcm != 1){
    /* r_cg2sun_J2000_normalized inertial to LVLH */
    compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, r_i2cg_INRTL, v_i2cg_INRTL);
    m_x_v(r_cg2sun_LVLH_normalized, T_inrtl_2_lvlh, r_cg2sun_J2000_normalized);

    /* r_cg2sun_J2000_normalized LVLH to body */
    compute_T_sc_to_lvlh( T_sc_to_lvlh, v_angle, order_rotation, INTEGRATOR->attitude.attitude_profile, &et,  r_i2cg_INRTL, v_i2cg_INRTL, INTEGRATOR->file_is_quaternion, INTEGRATOR->attitude.quaternion_current, PARAMS);
    //  compute_T_sc_to_lvlh(T_sc_to_lvlh, INTEGRATOR->attitude.lvlh_alongtrack_in_body_cartesian, INTEGRATOR->attitude.lvlh_crosstrack_in_body_cartesian, &et, r_i2cg_INRTL, v_i2cg_INRTL, INTEGRATOR); 
    m_trans(T_lvlh_to_sc, T_sc_to_lvlh);
    m_x_v(r_cg2sun_SC_normalized, T_lvlh_to_sc, r_cg2sun_LVLH_normalized ); 
    






  // openGL computation of cross-section area
  if (INTEGRATOR->opengl == 1){ 
    // // compute the azimuth and elevation of the sun vector in body frame. 
    // // Actually the spherical coordinates are used. 
    // // theta is the angle between z_body and the sun vector (from 0 to -180 deg). 
    // // phi is the angle between x_body and the projection of the sun vector on the (x_body, y_body) plane (from 0 to 360 deg).
    double x_body[3], y_body[3], z_body[3];
    x_body[0] = 1; x_body[1] = 0; x_body[2] = 0;
    y_body[0] = 0; y_body[1] = 1; y_body[2] = 0;
    z_body[0] = 0; z_body[1] = 0; z_body[2] = 1;
    double theta, phi;
    // // // theta
    theta = acos(r_cg2sun_SC_normalized[2]); // r_cg2sun_SC_normalized[2] is euqal to the projection of r_cg2sun_SC_normalized on z_body
    // // // phi
    phi = acos(r_cg2sun_SC_normalized[0]);
    if (r_cg2sun_SC_normalized[1] < 0){
      phi = 2*M_PI - phi;
    }


    double theta_deg, phi_deg;
    theta_deg = theta * RAD2DEG;
    phi_deg = phi * RAD2DEG;
    // // get the cross-section area from the file generated by openGL
    // // // compute the indices itheta and iphi
    int itheta, iphi, iface;
    double A_ref = 0;
    itheta = (int)((theta_deg - CONSTELLATION->area_attitude_opengl_theta0) / CONSTELLATION->area_attitude_opengl_dtheta);
    iphi = (int)((phi_deg - CONSTELLATION->area_attitude_opengl_phi0) / CONSTELLATION->area_attitude_opengl_dphi);

    A_ref = CONSTELLATION->area_attitude_opengl_total[itheta][iphi];
	
	
/*       printf("pres %f %f\n", theta_deg, phi_deg); */
/*       printf("pres theta %.2f, phi %.2f\n", itheta*CONSTELLATION->area_attitude_opengl_dtheta + CONSTELLATION->area_attitude_opengl_theta0, iphi*CONSTELLATION->area_attitude_opengl_dphi + CONSTELLATION->area_attitude_opengl_phi0); */
/*       printf("pres A_ref %f\n" , A_ref*1e10); */
      // we assume all surfaces have the same solar radiation coefficient so take the one of surface 0
	  // !!!!!!!!!!! THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES). Also (finally), if want to use Vallado, need to do some work to get the nomal of each face from the file OPTIONS->filename_area_attitude_opengl. Basically the same was you do for the drag coefficient, for wchih you need to know the normal of each triangle. I just didn't do it for the solar radiation pressure because wehn I'm writing these lines, i'm susing the equation from STK, wchih does not require to inkow the normal of each triangle (just the project area). 

	  a_solar_pressure_in_body[0] = - k * INTEGRATOR->surface[0].solar_radiation_coefficient * A_ref / INTEGRATOR->mass * solar_luminosity / (4 * M_PI * light_speed * dist_sat_to_sun * dist_sat_to_sun * 1000000.) * r_cg2sun_SC_normalized[0];
	  a_solar_pressure_in_body[1] = - k * INTEGRATOR->surface[0].solar_radiation_coefficient * A_ref / INTEGRATOR->mass * solar_luminosity / (4 * M_PI * light_speed * dist_sat_to_sun * dist_sat_to_sun * 1000000.) * r_cg2sun_SC_normalized[1];
	  a_solar_pressure_in_body[2] = - k * INTEGRATOR->surface[0].solar_radiation_coefficient * A_ref / INTEGRATOR->mass * solar_luminosity / (4 * M_PI * light_speed * dist_sat_to_sun * dist_sat_to_sun * 1000000.) * r_cg2sun_SC_normalized[2];

	  // !!!!!!!!!!! end of THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES). Also (finally) need to do somw rok to get the nomal of each face from the file OPTIONS->filename_area_attitude_opengl. Basically the same wasy you do for the drag coefficient, for wchih you need to know the normal of each triangle. I just didn't do it for the solar radiation pressure because wehn I'm writing these lines, i'm susing the equation from STK, wchih does not require to inkow the normal of each triangle (just the project area). 
    
	    // !!!!!!!!!! THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES). Also (finally) need to do somw rok to get the nomal of each face from the file OPTIONS->filename_area_attitude_opengl. Basically the same wasy you do for the drag coefficient, for wchih you need to know the normal of each triangle. I just didn't do it for the solar radiation pressure because wehn I'm writing these lines, i'm susing the equation from STK, wchih does not require to inkow the normal of each triangle (just the project area). So basically, don't just uncmment this block if ou want to use Vallado's equations.
	  /* term1 = INTEGRATOR->surface[0].diffuse_reflectivity / 3.0 + INTEGRATOR->surface[0].specular_reflectivity * cos_phi; */
	  /* term2 = 1 - INTEGRATOR->surface[0].specular_reflectivity; */

/* 	  a_solar_pressure_in_body[0] = - solar_flux * INTEGRATOR->surface[0].area * ( 2 * term1 * INTEGRATOR->surface[0].normal[0] + term2 * r_cg2sun_SC_normalized[0] ) / (light_speed * INTEGRATOR->mass); */
/* 	  a_solar_pressure_in_body[1] = - solar_flux * INTEGRATOR->surface[0].area * cos_phi * ( 2 * term1 * INTEGRATOR->surface[0].normal[1] + term2 * r_cg2sun_SC_normalized[1] ) / (light_speed * INTEGRATOR->mass); */
/* 	  a_solar_pressure_in_body[2] = - solar_flux * INTEGRATOR->surface[0].area * cos_phi * ( 2 * term1 * INTEGRATOR->surface[0].normal[2] + term2 * r_cg2sun_SC_normalized[2] ) / (light_speed * INTEGRATOR->mass); */
	    // !!!!!!!!!! END OF THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES).Also (finally) need to do somw rok to get the nomal of each face from the file OPTIONS->filename_area_attitude_opengl. Basically the same wasy you do for the drag coefficient, for wchih you need to know the normal of each triangle. I just didn't do it for the solar radiation pressure because wehn I'm writing these lines, i'm susing the equation from STK, wchih does not require to inkow the normal of each triangle (just the project area). So basically, don't just uncmment this block if ou want to use Vallado's equations.


  } // end of if (OPTIONS->opengl == 1)

  
  else{



    for (sss = 0; sss < INTEGRATOR->nb_surfaces; sss++){
      if (sss == 0){
	
	v_dot(&cos_phi, r_cg2sun_SC_normalized, INTEGRATOR->surface[0].normal);
	if (cos_phi > 0){

	  // !!!!!!!!!!! THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  a_solar_pressure_in_body[0] = - k * INTEGRATOR->surface[0].solar_radiation_coefficient * INTEGRATOR->surface[0].area * cos_phi / INTEGRATOR->mass * solar_luminosity / (4 * M_PI * light_speed * dist_sat_to_sun * dist_sat_to_sun * 1000000.) * r_cg2sun_SC_normalized[0];
	  a_solar_pressure_in_body[1] = - k * INTEGRATOR->surface[0].solar_radiation_coefficient * INTEGRATOR->surface[0].area * cos_phi / INTEGRATOR->mass * solar_luminosity / (4 * M_PI * light_speed * dist_sat_to_sun * dist_sat_to_sun * 1000000.) * r_cg2sun_SC_normalized[1];
	  a_solar_pressure_in_body[2] = - k * INTEGRATOR->surface[0].solar_radiation_coefficient * INTEGRATOR->surface[0].area * cos_phi / INTEGRATOR->mass * solar_luminosity / (4 * M_PI * light_speed * dist_sat_to_sun * dist_sat_to_sun * 1000000.) * r_cg2sun_SC_normalized[2];

	  // !!!!!!!!!!! END OF THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm).  UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)

	    // !!!!!!!!!! THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  /* term1 = INTEGRATOR->surface[0].diffuse_reflectivity / 3.0 + INTEGRATOR->surface[0].specular_reflectivity * cos_phi; */
	  /* term2 = 1 - INTEGRATOR->surface[0].specular_reflectivity; */

	  /* a_solar_pressure_in_body[0] = - solar_flux * INTEGRATOR->surface[0].area * cos_phi * ( 2 * term1 * INTEGRATOR->surface[0].normal[0] + term2 * r_cg2sun_SC_normalized[0] ) / (light_speed * INTEGRATOR->mass); */
	  /* a_solar_pressure_in_body[1] = - solar_flux * INTEGRATOR->surface[0].area * cos_phi * ( 2 * term1 * INTEGRATOR->surface[0].normal[1] + term2 * r_cg2sun_SC_normalized[1] ) / (light_speed * INTEGRATOR->mass); */
	  /* a_solar_pressure_in_body[2] = - solar_flux * INTEGRATOR->surface[0].area * cos_phi * ( 2 * term1 * INTEGRATOR->surface[0].normal[2] + term2 * r_cg2sun_SC_normalized[2] ) / (light_speed * INTEGRATOR->mass); */
	    // !!!!!!!!!! END OF THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	}
	else {
	  a_solar_pressure_in_body[0] = 0.0;
	  a_solar_pressure_in_body[1] = 0.0;
	  a_solar_pressure_in_body[2] = 0.0;
	}
      }
      else{
      
	v_dot(&cos_phi, r_cg2sun_SC_normalized, INTEGRATOR->surface[sss].normal);
      
	if (cos_phi > 0){
	  // !!!!!!!!!!! THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  a_solar_pressure_in_body[0] = a_solar_pressure_in_body[0]  - k * INTEGRATOR->surface[sss].solar_radiation_coefficient * INTEGRATOR->surface[sss].area * cos_phi / INTEGRATOR->mass * solar_luminosity / (4 * M_PI * light_speed * dist_sat_to_sun * dist_sat_to_sun * 1000000.) * r_cg2sun_SC_normalized[0];
	  a_solar_pressure_in_body[1] = a_solar_pressure_in_body[1] - k * INTEGRATOR->surface[sss].solar_radiation_coefficient * INTEGRATOR->surface[sss].area * cos_phi / INTEGRATOR->mass * solar_luminosity / (4 * M_PI * light_speed * dist_sat_to_sun * dist_sat_to_sun * 1000000.) * r_cg2sun_SC_normalized[1];
	  a_solar_pressure_in_body[2] = a_solar_pressure_in_body[2] - k * INTEGRATOR->surface[sss].solar_radiation_coefficient * INTEGRATOR->surface[sss].area * cos_phi / INTEGRATOR->mass * solar_luminosity / (4 * M_PI * light_speed * dist_sat_to_sun * dist_sat_to_sun * 1000000.) * r_cg2sun_SC_normalized[2];
	  // !!!!!!!!!!! END OF THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm).  UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)

	    // !!!!!!!!!! THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS.ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  /* term1 = INTEGRATOR->surface[sss].diffuse_reflectivity / 3.0 + INTEGRATOR->surface[sss].specular_reflectivity * cos_phi; */
	  /* term2 = 1 - INTEGRATOR->surface[sss].specular_reflectivity; */

	  /* a_solar_pressure_in_body[0] = a_solar_pressure_in_body[0] - solar_flux * INTEGRATOR->surface[sss].area * cos_phi * ( 2 * term1 * INTEGRATOR->surface[sss].normal[0] + term2 * r_cg2sun_SC_normalized[0] ) / (light_speed * INTEGRATOR->mass); */
	  /* a_solar_pressure_in_body[1] = a_solar_pressure_in_body[1] - solar_flux * INTEGRATOR->surface[sss].area * cos_phi * ( 2 * term1 * INTEGRATOR->surface[sss].normal[1] + term2 * r_cg2sun_SC_normalized[1] ) / (light_speed * INTEGRATOR->mass); */
	  /* a_solar_pressure_in_body[2] = a_solar_pressure_in_body[2] - solar_flux * INTEGRATOR->surface[sss].area * cos_phi * ( 2 * term1 * INTEGRATOR->surface[sss].normal[2] + term2 * r_cg2sun_SC_normalized[2] ) / (light_speed * INTEGRATOR->mass); */
	    // !!!!!!!!!! END OF THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	}
   
      }
   
    }

    m_x_v(a_solar_pressure_in_LVLH, T_sc_to_lvlh, a_solar_pressure_in_body);
    m_trans(T_lvlh_2_inrtl, T_inrtl_2_lvlh);
    m_x_v(a_solar_pressure_INRTL, T_lvlh_2_inrtl, a_solar_pressure_in_LVLH);  

  } // end oif opengl != 1
  } // end of no collision with VCM as colllision input file

  else{

    double ps = solar_luminosity / (4 * M_PI * light_speed * 1000);// !!!!!!!!!!!!!! vlaue of ps? here i compute it using the other expresssion of the solar radiatino pressure but it might be something different in heujdu and ghrist, 2011 ( see ref below). light_speed coverted from km/s to m/s

    v_scale(a_solar_pressure_INRTL, r_cg2sun_J2000_normalized,   - ps *  INTEGRATOR->srp_vcm * 1000. * 1000.  / ( dist_sat_to_sun * dist_sat_to_sun * 1000 * 1000)); /// used equation (3) of "solar radiation pressure binning for the geosyncrhonous orbit" by M.D. Hejduk and R.W. Ghrist, 2011 (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20110015238.pdf). srp_vcm converted from km2/kg to m2/kg 

    a_solar_pressure_INRTL[0] = a_solar_pressure_INRTL[0] / 1000.; // m/s2 to km/s2
    a_solar_pressure_INRTL[1] = a_solar_pressure_INRTL[1] / 1000.;
    a_solar_pressure_INRTL[2] = a_solar_pressure_INRTL[2] / 1000.;
    //    printf("%d %d %e %e\n", INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb, INTEGRATOR->srp_vcm, dist_sat_to_sun);
  }

  //  v_norm_print(a_solar_pressure_INRTL, "a_solar_pressure_INRTL");
  }
  else{
    a_solar_pressure_INRTL[0] = 0.0;
    a_solar_pressure_INRTL[1] = 0.0;
    a_solar_pressure_INRTL[2] = 0.0;

  }

  //   v_norm_print(a_solar_pressure_INRTL, "a_solar_pressure_INRTL")  ;
  return 0;

}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           compute_earth_pressure
//  Purpose:        Computes earth pressure acceleration
//  Assumptions:    spherical Earth
//  References      /
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|---------------------------------------------------------
//      | C. Bussy-Virat| 01/26/2019    |   ---     | Initial implementation
//
/////////////////////////////////////////////////////////////////////////////////////////

int compute_earth_pressure(double          a_earth_pressure_INRTL[3],
			   double          r_i2cg_INRTL[3],
			   double          v_i2cg_INRTL[3],
			   double          et,
			   PARAMS_T        *PARAMS,
			   INTEGRATOR_T    *INTEGRATOR,
			   CONSTELLATION_T *CONSTELLATION,
			   double          et_initial_epoch,
			   double          et_sc_initial,
			   int             index_in_attitude_interpolated)
		
{


  // Declarations
  int nelev_interpo, nazim_interpo, nzenith_interpo, nradius_interpo;
  	double elev_surf, azim_surf;
	int iradius_arr[2], iazim_surf_arr[2], izenith_arr[2], ielev_surf_arr[2];
	int iradius0, iradius1, iradius2, iradius3;
        int izenith0, izenith1, izenith2, izenith3;
	int ielev_surf0, ielev_surf1, ielev_surf2, ielev_surf3;
        int iazim_surf0, iazim_surf1, iazim_surf2, iazim_surf3;
	double azim_surf_corr;

          double a_earth_pressure_fac_per_surf[3];
          double a_earth_pressure_fac_temp[3];
	  a_earth_pressure_fac_temp[0] = 0; a_earth_pressure_fac_temp[1] = 0; a_earth_pressure_fac_temp[2] = 0;   
          double a_earth_pressure_fac[3];
	  a_earth_pressure_fac[0] = 0; a_earth_pressure_fac[1] = 0; a_earth_pressure_fac[2] = 0;   

	  int ii;
		
		double *xinter_radius, *xinter_zenith, *xinter_elev_surf, *xinter_azim_surf;
    double **yinter;
    int order_interpo_map = 2;
    xinter_radius = malloc(2 * sizeof(double));
    yinter = malloc(3 * sizeof(double *));
    int sss, jj;
    for (jj = 0; jj < 3; jj++){
      yinter[jj] = malloc(2 * sizeof(double));
    }
    xinter_zenith = malloc(2 * sizeof(double));
    xinter_elev_surf = malloc(2 * sizeof(double));
        xinter_azim_surf = malloc(2 * sizeof(double));

	int iel, izen, irad;
	double **y_radius_zenith_elev_surf, **y_radius_zenith, **y_radius;
	y_radius_zenith_elev_surf = malloc(2 * sizeof(double *));
	y_radius_zenith = malloc(2 * sizeof(double *));
	y_radius = malloc(2 * sizeof(double *));
	for (irad = 0; irad < 2; irad++){
	  y_radius[irad] = malloc(3 * sizeof(double));
	}
	for (izen = 0; izen < 2; izen++){
	  y_radius_zenith[izen] = malloc(3 * sizeof(double));
	}
	for (iel = 0; iel < 2; iel++){
	  y_radius_zenith_elev_surf[iel] = malloc(3 * sizeof(double));
	}

  double cr;
  double r_earth_elt_body[3];
  double r_earth_elt_lvlh_norm[3], r_earth_elt_lvlh_norm_minus[3];
  double cos_kappa;

  double T_lvlh_2_inrtl[3][3];
  // Determine the section of Earth visible from the satellite
  // // M: point on the surface of Earth, visible from the satellite
  // // M_max: point on the surface of the Earth at the limit of this section (ie, see the satellite at the horizon)
  // // O: center of the Earth
  // // S: satellite
  // //H: sub-satellite point on the surface of the Earth
  // // P: projection of M on OS line (P_max if M is M_max)
  // // beta: angle (OM, OS) (beta_max = angle  (OM_max, OS)
  // // alpha: angle (SM, SO) (alpha_max = angle  (SM_max, SO)
  // // mp: distance M to P
  // // sm: distance S to M (sm_max = distance S to M_max)
  // // sp: distance S to P (sp_max = distance S to P_max)
  // // sh: distance S to H
  double T_sc_to_lvlh[3][3];
  double earth_center_to_earth_elt_lvlh[3];
  double r_i2cg_lvlh[3];
  double earth_center_to_earth_elt_lvlh_norm[3];
  double cos_zenith;
  double a_earth_pressure_in_body[3];
  double radius_sc, beta, alpha, mp, sm, beta_max, alpha_max, sm_max, sp, oh;
  double delta, sh;
  v_mag( &radius_sc, r_i2cg_INRTL);
  double r_earth_elt_lvlh[3];
  double r_eci_norm[3];
  int row;
  double x[6];
  double lt;

  double earth_elt_to_sun_lvlh[3];
  double T_inrtl_2_lvlh[3][3];
  double r_cg2sun_LVLH[3];
  double r_cg2sun_J2000[3];
  double r_earth2sun_J2000[3];
  double earth_elt_to_sun_lvlh_norm[3];

  double v_angle[3];
  int order_rotation[3];
  double light_speed = 299792.458; // spped of light in km/s
  double solar_luminosity = 3.823e26; // in Watts
  double dist_sat_to_sun;
  double T_lvlh_to_sc[3][3];

  double cos_phi;

  if (INTEGRATOR->coll_vcm != 1){
    if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") == 0){
      v_angle[0] = INTEGRATOR->attitude.pitch_angular_velocity_ensemble * ( et - et_sc_initial ) + INTEGRATOR->attitude.pitch_ini_ensemble;
      v_angle[1] = INTEGRATOR->attitude.roll_angular_velocity_ensemble * ( et - et_sc_initial ) + INTEGRATOR->attitude.roll_ini_ensemble;
      v_angle[2] = INTEGRATOR->attitude.yaw_angular_velocity_ensemble * ( et - et_sc_initial ) + INTEGRATOR->attitude.yaw_ini_ensemble;
      order_rotation[0] = 1; // !!!!!!!! we might want to change that in the future                                                                     
      order_rotation[1] = 2;
      order_rotation[2] = 3;
      INTEGRATOR->attitude.pitch_current = v_angle[0];
      INTEGRATOR->attitude.roll_current = v_angle[1];
      INTEGRATOR->attitude.yaw_current = v_angle[2];

    }

    if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") == 0) {
      v_angle[0] =  INTEGRATOR->attitude.pitch_for_attitude_ensemble  +  INTEGRATOR->attitude.pitch_angular_velocity_constant * ( et - et_sc_initial );
      v_angle[1] =  INTEGRATOR->attitude.roll_for_attitude_ensemble  +  INTEGRATOR->attitude.roll_angular_velocity_constant * ( et - et_sc_initial );
      v_angle[2] =  INTEGRATOR->attitude.yaw_for_attitude_ensemble  +  INTEGRATOR->attitude.yaw_angular_velocity_constant * ( et - et_sc_initial );

      order_rotation[0]  = 1; order_rotation[1]  = 2; order_rotation[2]  = 3;
	   
      INTEGRATOR->attitude.pitch_current = v_angle[0];
      INTEGRATOR->attitude.roll_current = v_angle[1];
      INTEGRATOR->attitude.yaw_current = v_angle[2];
      INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
      INTEGRATOR->attitude.order_roll_current = order_rotation[1];
      INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


	      
    }


    if ( (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "sun_pointed") != 0) ){ // otherwise (atittude is nadir, sun_pointed or manual (from an input file))    
      //    index_in_attitude_interpolated = floor( ( et - et_sc_initial ) / ( INTEGRATOR->dt / 2.0) ) ; 
      if (INTEGRATOR->file_is_quaternion == 0){
	v_angle[0] = INTEGRATOR->attitude.pitch[index_in_attitude_interpolated];
	v_angle[1] = INTEGRATOR->attitude.roll[index_in_attitude_interpolated];
	v_angle[2] = INTEGRATOR->attitude.yaw[index_in_attitude_interpolated];
	order_rotation[0] = INTEGRATOR->attitude.order_pitch[index_in_attitude_interpolated];
	order_rotation[1] = INTEGRATOR->attitude.order_roll[index_in_attitude_interpolated];
	order_rotation[2] = INTEGRATOR->attitude.order_yaw[index_in_attitude_interpolated];
	INTEGRATOR->attitude.pitch_current = v_angle[0];
	INTEGRATOR->attitude.roll_current = v_angle[1];
	INTEGRATOR->attitude.yaw_current = v_angle[2];
        INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
	INTEGRATOR->attitude.order_roll_current = order_rotation[1];
	INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


      }
      else{
	q_copy( INTEGRATOR->attitude.quaternion_current, INTEGRATOR->attitude.quaternion[index_in_attitude_interpolated]);
      }
    }
  } // end of no collision with VCM as colllision input file


	  
  spkez_c(10, et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.

  r_earth2sun_J2000[0] = x[0];
  r_earth2sun_J2000[1] = x[1];
  r_earth2sun_J2000[2] = x[2];
  v_sub(r_cg2sun_J2000, r_earth2sun_J2000, r_i2cg_INRTL);
  v_mag( &dist_sat_to_sun, r_cg2sun_J2000 );
  double prad;
  prad = solar_luminosity / (4 * M_PI * light_speed * 1000 * dist_sat_to_sun * dist_sat_to_sun * 1000000.); // in N/m2. should be about 4.5e-6 N./m2 (the flux prad*c should be about 1350 W/m2)
  double zenith_sc, cos_zenith_sc, r_earth2sun_J2000_mag;
  v_mag(&r_earth2sun_J2000_mag, r_earth2sun_J2000);
  v_dot(&cos_zenith_sc, r_i2cg_INRTL, r_earth2sun_J2000);
  cos_zenith_sc = cos_zenith_sc / ( radius_sc *  r_earth2sun_J2000_mag);
  zenith_sc = acos(cos_zenith_sc); // doesn't matter if get +-zenith angle
  if (INTEGRATOR->coll_vcm != 1){
    compute_T_sc_to_lvlh( T_sc_to_lvlh, v_angle, order_rotation, INTEGRATOR->attitude.attitude_profile, &et,  r_i2cg_INRTL, v_i2cg_INRTL, INTEGRATOR->file_is_quaternion, INTEGRATOR->attitude.quaternion_current, PARAMS);
      double sc_normal_lvlh[3], sc_normal_eci[3], sc_normal_epf_temp[3], sc_normal_epf[3]; //epf: Earth pressure frame (see definition in function compute_T_inrtl_2_earth_pres_frame in prop_math.c)
      double T_lvlh_to_inrtl[3][3];
      double T_inrtl_2_earth_pres_frame[3][3];
      double T_earth_pres_frame_2_inrtl[3][3];
      compute_T_inrtl_2_earth_pres_frame( T_inrtl_2_earth_pres_frame, r_i2cg_INRTL, et);
      m_trans(T_earth_pres_frame_2_inrtl, T_inrtl_2_earth_pres_frame);
      compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, r_i2cg_INRTL, v_i2cg_INRTL);
      m_trans(T_lvlh_to_inrtl, T_inrtl_2_lvlh);
      // determine the bin for the radius
      compute_iradius_earth_pressure_map(xinter_radius, iradius_arr, &nradius_interpo, &PARAMS->EARTH.GRAVITY, radius_sc);
      // determine the bin for the zenith
      compute_izenith_earth_pressure_map(xinter_zenith, izenith_arr, &nzenith_interpo, &PARAMS->EARTH.GRAVITY, zenith_sc);
      for (sss = 0; sss < INTEGRATOR->nb_surfaces_eff; sss++){// !!!!! < INTEGRATOR->nb_surfaces; sss++){
	a_earth_pressure_fac_per_surf[0] = 0; a_earth_pressure_fac_per_surf[1] = 0; a_earth_pressure_fac_per_surf[2] = 0;   
	cr = INTEGRATOR->surface_eff[sss].solar_radiation_coefficient; // shorter notation
    	m_x_v(sc_normal_lvlh, T_sc_to_lvlh, INTEGRATOR->surface_eff[sss].normal);
	m_x_v(sc_normal_eci, T_lvlh_to_inrtl, sc_normal_lvlh);
	m_x_v(sc_normal_epf_temp, T_inrtl_2_earth_pres_frame, sc_normal_eci);
	v_norm(sc_normal_epf, sc_normal_epf_temp);// if sc_normal_epf has a norm not exactly equal to 1 then the acos right below can return nan (for example acos(-1.00000001) = nan)
	elev_surf = acos(sc_normal_epf[2]); // elev_surf varies from 0 to 180
	if (elev_surf*180/M_PI >= PARAMS->EARTH.GRAVITY.min_elev_surf_map){ // otherwise the surface doesn't see any part of the Earthzenith

	azim_surf = atan2(sc_normal_epf[1], sc_normal_epf[0]);
	// determine the bin for the elevation of the normal of the surface
	compute_ielev_surf_earth_pressure_map(xinter_elev_surf, ielev_surf_arr, &nelev_interpo,  &PARAMS->EARTH.GRAVITY, elev_surf);
	// determine the bin for the azimuth of the normal of the surface
	if (azim_surf >= 0){ // azim_surf_corr varies from 0 to 2*M_PI
	  azim_surf_corr = azim_surf;
	}
	else{
	  azim_surf_corr = 2*M_PI + azim_surf;
	}    
	compute_iazim_surf_earth_pressure_map(xinter_azim_surf, iazim_surf_arr, &nazim_interpo, &PARAMS->EARTH.GRAVITY, azim_surf_corr);
	// Interpolat over raidus, zenith of sc, elev and azim of sc surface
	for (irad = 0; irad < nradius_interpo; irad++){
	  for (izen = 0; izen < nzenith_interpo; izen++){
	    order_interpo_map = 2;
	    for (iel = 0; iel < nelev_interpo; iel++){
	      earth_pressure_map_yinter_azim_surf(yinter, &order_interpo_map, azim_surf_corr, &PARAMS->EARTH.GRAVITY, PARAMS->EARTH.GRAVITY.earth_pressure_map[iradius_arr[irad]][izenith_arr[izen]][ielev_surf_arr[iel]], iazim_surf_arr);
	      polynomial_interpo_earth(y_radius_zenith_elev_surf[iel], order_interpo_map, xinter_azim_surf, yinter, azim_surf_corr*180/M_PI);

	    }
	    // // // ->  zenith0 (interpo over elev_surf)

	    earth_pressure_map_yinter_elev_surf(yinter, &order_interpo_map, elev_surf, &PARAMS->EARTH.GRAVITY, y_radius_zenith_elev_surf);
	    polynomial_interpo_earth(y_radius_zenith[izen], order_interpo_map, xinter_elev_surf, yinter, elev_surf*180/M_PI);

	  }
	  // // -> radius0 (interpo over zenith)
	  earth_pressure_map_yinter_zenith(yinter, &order_interpo_map, zenith_sc, &PARAMS->EARTH.GRAVITY, y_radius_zenith);
	  polynomial_interpo_earth(y_radius[irad], order_interpo_map, xinter_zenith, yinter, zenith_sc*180/M_PI);
	}
      // Interpo over radius. function gravity_map_yinter_radius is used because it has the same purpose here
      gravity_map_yinter_radius_earth( yinter, &order_interpo_map,  radius_sc, &PARAMS->EARTH.GRAVITY,  y_radius);
	// Interpolate a_pressure_fac 
      polynomial_interpo_earth(a_earth_pressure_fac_temp, order_interpo_map, xinter_radius, yinter, radius_sc);
	

      
    //    } // end of go over each of the 3 coordinates
	  
	/* printf("\nradius %.0f %.0f %.2f %.0f %.0f\n", xinter_radius[0]-PARAMS->EARTH.GRAVITY.radius, xinter_radius[1]-PARAMS->EARTH.GRAVITY.radius, radius_sc-PARAMS->EARTH.GRAVITY.radius, xinter_radius[2]-PARAMS->EARTH.GRAVITY.radius, xinter_radius[3]-PARAMS->EARTH.GRAVITY.radius); */
	/* printf("zenith %.0f %.0f %.2f %.0f %.0f\n", xinter_zenith[0], xinter_zenith[1], zenith_sc*180/M_PI, xinter_zenith[2], xinter_zenith[3]); */
	/* printf("elev_surf %.0f %.0f %.2f %.0f %.0f\n", xinter_elev_surf[0], xinter_elev_surf[1], elev_surf*180/M_PI, xinter_elev_surf[2], xinter_elev_surf[3]); */
	/* printf("azim_surf %.0f %.0f %.2f %.0f %.0f\n", xinter_azim_surf[0], xinter_azim_surf[1], azim_surf_corr*180/M_PI, xinter_azim_surf[2], xinter_azim_surf[3]); */

    // a_earth_pressure_fac_temp is in unit of km^(-1), INTEGRATOR->surface[sss].area is km^2 so a_earth_pressure_fac_temp * INTEGRATOR->surface[sss].area is in km. prad is in N/m2 so convert a_earth_pressure_fac_temp * INTEGRATOR->surface[sss].area in m -> * 1000 in expresssion below
    v_scale(a_earth_pressure_fac_per_surf,  a_earth_pressure_fac_temp, prad * cr * INTEGRATOR->surface_eff[sss].area * 1000. / INTEGRATOR->mass);
	}
	  
	//    Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][2]

	//    printf("%d: elev %f, azim %f\n",sss,elev_surf*180/M_PI, azim_surf*180/M_PI);
	//  v_print(sc_normal_epf, "sc_normal_epf");
      a_earth_pressure_fac[0] = a_earth_pressure_fac[0] + a_earth_pressure_fac_per_surf[0];
      a_earth_pressure_fac[1] = a_earth_pressure_fac[1] + a_earth_pressure_fac_per_surf[1];
      a_earth_pressure_fac[2] = a_earth_pressure_fac[2] + a_earth_pressure_fac_per_surf[2];
      }

      m_x_v(a_earth_pressure_INRTL, T_earth_pres_frame_2_inrtl, a_earth_pressure_fac);
      //	  v_dot(&cos_phi, r_earth_elt_body_norm, INTEGRATOR->surface[sss].normal);

      //	      a_earth_pressure_in_body[0] = a_earth_pressure_in_body[0] - cr * albedo * cos_zenith * prad  * INTEGRATOR->surface[sss].area*1000000. * cos_phi / INTEGRATOR->mass * area_earth_elt / (M_PI * sm * sm) * r_earth_elt_body_norm[0]/ 1000.; // surface[sss].area in km2 -> m2, a   
    
}  // end of no collision with VCM as colllision input file
   /* v_print(a_earth_pressure_fac, "a_earth_pressure_fac"); */
   /*  v_norm_print(a_earth_pressure_INRTL, "a_earth_pressure_INRTL"); */
   /*  exitf(); */

  return 0;

}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           compute_drag
//  Purpose:        Computes acceleration acceleration using NRLMSIS 2000
//  Assumptions:    None.
//  References      NRL / GSFC
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|---------------------------------------------------------
//      | J. Getchius   | 05/20/2015    |   ---     | Initial Implementation
//      | C. Bussy-Virat| 07/19/2015    |   ---     | Corrections: add Cd as a factor in the drag acceleration
//      | C. Bussy-Virat| 07/28/2015    |   ---     | Add attiude dependence and multi surfaces + the third component of LVLH now points away from the Earth (and not towards the Earth anymore)
//
/////////////////////////////////////////////////////////////////////////////////////////

int compute_drag(       double          adrag_i2cg_INRTL[3],
                        double          r_i2cg_INRTL[3],
                        double          v_i2cg_INRTL[3],
                        double          et,
                        PARAMS_T        *PARAMS,
                        INTEGRATOR_T    *INTEGRATOR,
			double          et_initial_epoch,
			double          et_sc_initial,
			double          *density,
			int             index_in_attitude_interpolated, 
			int             index_in_driver_interpolated,
			CONSTELLATION_T *CONSTELLATION,
			OPTIONS_T       *OPTIONS,
			int iProc,
			int iDebugLevel,
			SPACECRAFT_T *SC)

{


  //  printf("%d %d\n",index_in_attitude_interpolated, index_in_driver_interpolated);

  double total_cross_area = 0;
  int iface, iface_here;
  if (iDebugLevel >= 2){
    printf("--- (compute_drag) Just got in compute_drag ... (iProc %d | iii = %d, eee = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb);
  }

  // if set to 0 then don't write in the given files. if set to 1 the write in given files


  // Declarations
  double cd;
  double cos_angle_v_sc_normal;
  double normal_in_ntw_normalized;
 
  double A_ref = 0;
  double A_ref_tot = 0;
  double cd_tot_norm =  0;

  int aaa_a;
  double nb_index_in_81_days;
  int hhh;
  //  int start_ensemble = 1;
  int eee;

  int ppp;
  double latitude_gitm,  longitude_gitm,  altitude_gitm;
  double r_ecef2cg_ECEF_gitm[3];
  double T_J2000_to_ECEF_gitm[3][3];
  double ballistic_coefficient = 0;
  double v_angle[3];
  int order_rotation[3];
  //      int index_in_attitude_interpolated;
  struct nrlmsise_output output;
  struct nrlmsise_input input;
  char   timestamp[36];
  char   doy_s[256];
  char   year_s[256];
  char   hour_s[256];
  char   min_s[256];
  char   sec_s[256];
  double hour = 0.0;
  double min = 0.0;
  double sec = 0.0;
  double r_ecef2cg_ECEF[3];
  double T_J2000_to_ECEF[3][3];
  double normal_in_lvlh[3];
  double normal_in_inertial[3];
  double normal_in_ntw[3];
  double v_i2cg_NTW_DOT_normal_to_the_surface_NTW;
  double altitude;
  double longitude;
  double latitude;
  double T_sc_to_lvlh[3][3];
  double T_inrtl_2_lvlh[3][3];
  double T_lvlh_2_inrtl[3][3];
  double T_inrtl_2_ntw[3][3];
  double T_ntw_2_inrtl[3][3];
  double v_i2cg_NTW[3];
  double a_i2cg_NTW[3];
  //  int ii;
  int sss;
  char timestamp_isoc[300];
  //   int index_in_driver_interpolated;

  //  double density;
  /////////////////// COMPUTE DENSITY ///////////////////
  if ( ( strcmp(INTEGRATOR->format_density_driver, "density_file") != 0 ) && ( strcmp(INTEGRATOR->format_density_driver, "gitm") != 0 ) ){ // the user chooses f107 and Ap for the density
    // Generate the time inputs
    
    et2utc_c(et, "D", 3, 35, timestamp);
    et2utc_c(et, "ISOC" ,6 ,255 , timestamp_isoc);



    // YEAR
    strcpy(year_s, "");
    strncat(year_s, &timestamp[0], 4);
    input.year = atoi( year_s );

    // DOY
    strcpy(doy_s, "");
    strncat(doy_s, &timestamp[5], 3);
    /* for (ii = 0; ii < 3; ii++ ) { */
    /*   doy_s[ii] = timestamp[ii+5]; */
    /* } */
    input.doy = atoi( doy_s );
    
    //  Hour
    strcpy(hour_s, "");
    strncat(hour_s, &timestamp[12], 2);
    /* for (ii = 0; ii < 2; ii++ ) { */
    /*   hour_s[ii] = timestamp[ii+12]; */
    /* } */
    hour = atof( hour_s );
    
    // Min
    strcpy(min_s, "");
    strncat(min_s, &timestamp[15], 2);
    /* for (ii = 0; ii < 2; ii++ ) { */
    
    /*   min_s[ii] = timestamp[ii+15]; */
    /* } */
    min = atof( min_s );

    // Sec
    strcpy(sec_s, "");
    strncat(sec_s, &timestamp[18], 6);
    /* for (ii = 0; ii < 6; ii++ ) { */
    
    /*   sec_s[ii] = timestamp[ii+18]; */
    
    /* } */
    sec = atof( sec_s );

    input.sec = hour * 3600.0 + min * 60.0 + sec;
    //    printf("%s: %d %d %f\n", timestamp, input.year, input.doy, input.sec);

    /*   double geodetic[3]; */
    /* 	double altitude, latitude, longitude; */
    /* 	eci2lla(r_i2cg_INRTL , et, geodetic ); */
    /* 	altitude = geodetic[2]; */
    /* 	latitude = geodetic[0]; */
    /* 	longitude = geodetic[1];  */

    /* 	// update the planet fixed state		 */
    /* 	  geodetic_to_geocentric(PARAMS->EARTH.flattening,             */
    /* 				 altitude, */
    /* 				 latitude, */
    /* 				 longitude, */
    /* 				 PARAMS->EARTH.radius,        */
    /* 				 r_ecef2cg_ECEF) ; */
    // Geneterate the Geodetic inputs
    pxform_c("J2000", PARAMS->EARTH.earth_fixed_frame, et, T_J2000_to_ECEF);
    m_x_v(r_ecef2cg_ECEF, T_J2000_to_ECEF, r_i2cg_INRTL);
    geocentric_to_geodetic( r_ecef2cg_ECEF,
			    &PARAMS->EARTH.radius,
			    &PARAMS->EARTH.flattening,
			    &altitude,
			    &latitude,
			    &longitude);
    input.g_lat     = latitude * RAD2DEG;
    input.g_long    = longitude * RAD2DEG;
    input.alt       = altitude;
    input.lst       = input.sec/3600. + input.g_long/15.;
    //    printf("%f %f %f %f ||| ", input.g_lat ,input.g_long  ,input.alt , input.lst);
    // Solar / Geomagnetic Activity
    if ( strcmp(INTEGRATOR->format_density_driver, "static") == 0 ){ // if the user chooses a constant f107 and Ap for the density  
      if (iDebugLevel >= 3){
	printf("---- (compute_drag) Computing the static F10.7, F10.7A, and Ap to determine density for drag ... (iProc %d | iii = %d, eee = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb);
      }


      input.f107A = INTEGRATOR->f107A_static;
      input.f107  = INTEGRATOR->f107_static;
      //      printf("%s: %f %f", timestamp_isoc, input.f107,input.f107A);
      if (PARAMS->ATMOSPHERE.flags.switches[9] != -1){ // if daily ap
	
	input.ap    = INTEGRATOR->Ap_static;
	//	printf(" %f", INTEGRATOR->Ap_static);
      }
      else{ //if historical ap
	input.ap_a = malloc(7 * sizeof(double));
	for (ppp = 0; ppp < 7 ; ppp ++){
	  input.ap_a->a[ppp]    = INTEGRATOR->Ap_hist_static[ppp];
	}
	//	printf(" %f", input.ap_a->a[ppp]);
      }
      //      printf("\n");
      if (iDebugLevel >= 3){
	printf("---- (compute_drag) Done computing the static F10.7, F10.7A, and Ap to determine density for drag ... (iProc %d | iii = %d, eee = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb);
      }

    } // end of if the user chooses a constant f107 and Ap for the density  

    else{ // if the user chooses a time varying f107 and Ap for the density  
      //      index_in_driver_interpolated = floor( ( et - et_initial_epoch ) / ( INTEGRATOR->dt / 2.0 ) ) ; // "/ 2.0" because of the Runge Kunta orfer 4 method
      if (iDebugLevel >= 3){
	printf("---- (compute_drag) Computing F10.7, F10.7A, and Ap to determine density for drag ... (iProc %d | iii = %d, eee = %d | index_in_driver_interpolated = %d, CONSTELLATION->aaa_mod[INTEGRATOR->sc_main_nb] = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb, 	  index_in_driver_interpolated, CONSTELLATION->aaa_mod[INTEGRATOR->sc_main_nb]  );
      }

      if ( (OPTIONS->nb_ensembles_density > 0) && (OPTIONS->swpc_need_predictions) && ( et >= OPTIONS->swpc_et_first_prediction ) ) { 
	if (iDebugLevel >= 4){
	  printf("---- (compute_drag) Computing the ensembles on the predictions of F10.7, F10.7A, and Ap ...(iProc %d | iii = %d, eee = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb);
	}

	// if the user chose to run ensembles on the density using data from SWPC 
	// if future predictions of F10.7 and Ap (not only past observations) (past values (value before the current time at running) are perfectly known (result from observations, not predictions))
	// only overwrite predictions of ensembles (not observations). OPTIONS->et_interpo[aaa] corresponds to the first say of predictions
	//      printf("XXXXXXXXXXXXXX\nXXXXXXXXXXXXXXXXXXX\n");
	nb_index_in_81_days =  81. * 24 * 3600 / (OPTIONS->dt/2.) + 1 ;
	if (INTEGRATOR->sc_ensemble_nb == 1 + iProc * OPTIONS->nb_ensemble_min_per_proc){ // initialize array only once per iProc
	  // // Generate nb_ensembles_density normal random values
	  /* if (iProc == 0){ */
	  /* 	etprint(et, ""); */
	  /* 	printf("OPTIONS->sigma_ap[%d] at index %d: %f | %f\n",CONSTELLATION->aaa_sigma[INTEGRATOR->sc_main_nb], index_in_driver_interpolated, OPTIONS->sigma_ap[CONSTELLATION->aaa_sigma[INTEGRATOR->sc_main_nb]], OPTIONS->sigma_f107[CONSTELLATION->aaa_sigma[INTEGRATOR->sc_main_nb]]); */
	  /* 	} */
	  /* if (INTEGRATOR->sc_ensemble_nb == 1){ */
	  /* 	etprint(et, "time"); */
	  /* 	printf("%f %d\n", OPTIONS->sigma_ap[CONSTELLATION->aaa_sigma[INTEGRATOR->sc_main_nb]], INTEGRATOR->sc_main_nb); */
	  /* } */
	    

	  for ( eee = 0; eee < OPTIONS->nb_ensemble_min_per_proc; eee++){
	    CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time[iProc][eee]   = randn( OPTIONS->f107[index_in_driver_interpolated], OPTIONS->sigma_f107[CONSTELLATION->aaa_sigma[INTEGRATOR->sc_main_nb]]);
	    CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time[iProc][eee]   = randn(  OPTIONS->Ap[index_in_driver_interpolated], OPTIONS->sigma_ap[CONSTELLATION->aaa_sigma[INTEGRATOR->sc_main_nb]]);
	    if (CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time[iProc][eee] < 0){ // just in case sigma_ap is too big...
	      CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time[iProc][eee] =  0; // OPTIONS->Ap[index_in_driver_interpolated];
	    }
	    if (CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time[iProc][eee] < 0){// just in case sigma_f107 is too big...   
	      CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time[iProc][eee]  = 0; //OPTIONS->f107[index_in_driver_interpolated];
	    }
	    /* if (eee == 0){ */
	    /* etprint(OPTIONS->et_interpo[aaa], "time"); */
	    /* printf("Ap[%d]: %f | sigma_ap[%d]: %f \n",aaa, OPTIONS->Ap[aaa],CONSTELLATION->aaa_sigma[INTEGRATOR->sc_main_nb], OPTIONS->sigma_ap[CONSTELLATION->aaa_sigma[INTEGRATOR->sc_main_nb]]); */
	    /* } */
	  }

	  // // Order values in ascending order
	  sort_asc_order(CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time_sorted[iProc],  CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time[iProc], OPTIONS->nb_ensemble_min_per_proc);
	  sort_asc_order(CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time_sorted[iProc],  CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time[iProc], OPTIONS->nb_ensemble_min_per_proc);


	

	  // Initialization of swpc_et_first_prediction for the calculation of F10.7A for an ensemble
	  if ( et_sc_initial > OPTIONS->swpc_et_first_prediction){ // if the propagation starts in the future. Otherwise, sum_sigma_for_f107_average = 0 for all ensemble sc

	    if (  CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] == 0 )  {// for the first time that variable et gets bigger than swpc_et_first_prediction 
	      // initialize sum_sigma_for_f107_average as the sum of all sigma on F10.7 from first prediction to inital epoch. There is no mathematical logic in here. The exact solution would be to sum all the deviations between f107_ensemble and f107_refce from first prediciton until intial epoch, and sum them up. This is KIND OF similar. The reason I don't do the correct approach is that I did not calculate f107_ensemble for times before initial epoch
	      for (aaa_a = 0; aaa_a < CONSTELLATION->aaa_sigma[INTEGRATOR->sc_main_nb]; aaa_a++){

		CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] = CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] + OPTIONS->sigma_f107[aaa_a];		    
		//		      ptd( OPTIONS->sigma_f107[aaa_a], "s");
	      }

	      //		      printf("eee: %d | index: %d | iProc: %d | sum: %f\n", INTEGRATOR->sc_ensemble_nb, index_in_driver_interpolated, iProc, CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb]);exit(0);
		   
	    } // end of for the first time that variable et gets bigger than swpc_et_first_prediction
	  }// end of if the propagation starts in the future  
	  // End of initialization of swpc_et_first_prediction for the calculation of F10.7A for an ensemble

	} // end of initialize array only once per iProc
	else if (INTEGRATOR->sc_ensemble_nb == 0){ // also need to calculate f107 and ap for reference sc (no perturbation so directly equal to the values in the prediction files)
	  INTEGRATOR->Ap[index_in_driver_interpolated]            = OPTIONS->Ap[index_in_driver_interpolated];           // magnetic index(daily)
	  if (OPTIONS->use_ap_hist == 1){
	    for (hhh = 0; hhh < 7; hhh++){
	      INTEGRATOR->Ap_hist[hhh][index_in_driver_interpolated]            = OPTIONS->Ap_hist[hhh][index_in_driver_interpolated];           // magnetic index(historical)
	    }
	  }

	  INTEGRATOR->f107[index_in_driver_interpolated]          = OPTIONS->f107[index_in_driver_interpolated];         // Daily average of F10.7 flux
	  /* print_test(); */
	  /* printf("INTEGRATOR->f107[%d]: %f\n", index_in_driver_interpolated, INTEGRATOR->f107[index_in_driver_interpolated]); */

	  INTEGRATOR->f107A[index_in_driver_interpolated]         = OPTIONS->f107A[index_in_driver_interpolated];        // 81 day average of F10.7 flux


	} // end of also need to alculate f107 and ap for reference sc (no perturbation so directly equal to the values in the prediction files)
	if (INTEGRATOR->sc_ensemble_nb != 0) { // don't overwrite previously written values of f107 and Ap for reference sc
	  INTEGRATOR->Ap[index_in_driver_interpolated]            = CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time_sorted[iProc][INTEGRATOR->sc_ensemble_nb-1-iProc*OPTIONS->nb_ensemble_min_per_proc];
	  INTEGRATOR->f107[index_in_driver_interpolated]          = CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time_sorted[iProc][INTEGRATOR->sc_ensemble_nb-1-iProc*OPTIONS->nb_ensemble_min_per_proc]; 
	
	  CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] = CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] +  ( INTEGRATOR->f107[index_in_driver_interpolated]  - OPTIONS->f107[index_in_driver_interpolated] );
	  INTEGRATOR->f107A[index_in_driver_interpolated]         = OPTIONS->f107A[index_in_driver_interpolated] + CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] / nb_index_in_81_days; // derivation of F10.7A considering uncerainties in F10.7
	  //	    printf("eee: %d | index: %d | iProc: %d | sum: %f | 81: %f | opeion %f\n", INTEGRATOR->sc_ensemble_nb, index_in_driver_interpolated, iProc, CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb], nb_index_in_81_days, OPTIONS->dt);
	  //	    print_test();
	  /* if (iProc == 0){ */
	  /*   //	    if ((INTEGRATOR->sc_ensemble_nb == 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + 1)){ */
	  /*   if (index_in_driver_interpolated == 4){ */
	  /*     printf("eee: %d | f107:  %f | index: %d | iProc: %d\n", INTEGRATOR->sc_ensemble_nb, INTEGRATOR->f107[index_in_driver_interpolated], index_in_driver_interpolated, iProc); */
	  /* } */
	  /* } */
	} // end of don't overwrite previously written values of f107 and Ap for reference sc
	if (iDebugLevel >= 4){
	  printf("---- (compute_drag) Done computing the ensembles on the predictions of F10.7, F10.7A, and Ap ...(iProc %d | iii = %d, eee = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb);
	}

      } // end of:
      // if the user chose to run ensembles on the density using data from SWPC AND
      // if future predictions of F10.7 and Ap (not only past observations) (past values (value before the current time at running) are perfectly known (result from observations, not predictions))
      // only overwrite predictions of ensembles (not observations). OPTIONS->et_interpo[aaa] corresponds to the first say of predictions
      else{ // if we don't run ensembles on F10.7/Ap or that we run ensemble on F10.7/Ap but that the current time is before the first prediction (so there is no uncertainty in F10.7/Ap because the time corresponds to an observation, not a prediction)
	//	    	    printf("%d %d %d %d\n", index_in_driver_interpolated, iProc, INTEGRATOR->sc_ensemble_nb, INTEGRATOR->sc_main_nb);
	    //	    print_test();

	if ( (strcmp(OPTIONS->test_omniweb_or_external_file, "swpc_mod") == 0) && (OPTIONS->swpc_need_predictions) && ( et >= OPTIONS->swpc_et_first_prediction ) ){ //if the option "swpc_mod" (so nominal f107 and Ap minus a certain value), that predictions of f107 and ap are used, and that the current time is in the future
	  nb_index_in_81_days =  81. * 24 * 3600 / (OPTIONS->dt/2.) + 1 ;


	  INTEGRATOR->f107[index_in_driver_interpolated]  = OPTIONS->f107[index_in_driver_interpolated] + OPTIONS->mod_f107[CONSTELLATION->aaa_mod[INTEGRATOR->sc_main_nb]];
	  INTEGRATOR->Ap[index_in_driver_interpolated]    = OPTIONS->Ap[index_in_driver_interpolated] + OPTIONS->mod_ap[CONSTELLATION->aaa_mod[INTEGRATOR->sc_main_nb]]; // note: never ap_hist option if spwc prediction
/* 	if (INTEGRATOR->sc_ensemble_nb == 0){ */
/* 	etprint(et, "time"); */
/* 	printf("%d %d %f %d\n", INTEGRATOR->sc_main_nb, index_in_driver_interpolated,  INTEGRATOR->Ap[index_in_driver_interpolated], CONSTELLATION->aaa_mod[INTEGRATOR->sc_main_nb]); */
/* 	} */


	  /* if (iProc == 0){ */
	  /* 	if (INTEGRATOR->sc_main_nb==0){ */
	  /* 	  etprint(et, ""); */
	  /* 	  printf("(%d %d) %f %f (%d)\n", INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb,  INTEGRATOR->f107[index_in_driver_interpolated] , OPTIONS->mod_f107[CONSTELLATION->aaa_mod[INTEGRATOR->sc_main_nb]], CONSTELLATION->aaa_mod[INTEGRATOR->sc_main_nb] ); */

	  /* 	  } */
	  /* } */


	  /* if (iProc == 3){ */
	  /* 	etprint(et, ""); */
	  /* 	printf("%f\n",OPTIONS->mod_ap[CONSTELLATION->aaa_mod[INTEGRATOR->sc_main_nb]]); */
	  /* 	} */
	  // Initialization of swpc_et_first_prediction for the calculation of F10.7A for an ensemble
	  if ( et_sc_initial > OPTIONS->swpc_et_first_prediction){ // if the propagation starts in the future. Otherwise, sum_sigma_for_f107_average = 0 for all ensemble sc

	    if (  CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] == 0 )  {// for the first time that variable et gets bigger than swpc_et_first_prediction 
	      // initialize sum_sigma_for_f107_average as the sum of all mod on F10.7 from first prediction to inital epoch.
	      for (aaa_a = 0; aaa_a < CONSTELLATION->aaa_mod[INTEGRATOR->sc_main_nb]; aaa_a++){

		CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] = CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] + OPTIONS->mod_f107[aaa_a];		    
		//		      ptd( OPTIONS->mod_f107[aaa_a], "s");
	      }

	      //		      printf("eee: %d | index: %d | iProc: %d | sum: %f\n", INTEGRATOR->sc_ensemble_nb, index_in_driver_interpolated, iProc, CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb]);exit(0);
		   
	    } // end of for the first time that variable et gets bigger than swpc_et_first_prediction
	  }// end of if the propagation starts in the future  

	
	  CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] = CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] +  ( INTEGRATOR->f107[index_in_driver_interpolated]  - OPTIONS->f107[index_in_driver_interpolated] );
	  INTEGRATOR->f107A[index_in_driver_interpolated]         = OPTIONS->f107A[index_in_driver_interpolated] + CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] / nb_index_in_81_days; // derivation of F10.7A considerin

	  /* if (iProc == 3){ */
	  /* 	if (INTEGRATOR->sc_main_nb == 1){ */
	  /* 	etprint(et, ""); */
	  /* 	//		printf("%f (%d  %d) \n", INTEGRATOR->f107[index_in_driver_interpolated], index_in_driver_interpolated, INTEGRATOR->sc_ensemble_nb); */
	  /* 		      	printf("%f %f %f | %f  (%d %d)\n",OPTIONS->f107A[index_in_driver_interpolated], INTEGRATOR->f107[index_in_driver_interpolated], OPTIONS->f107[index_in_driver_interpolated] ,CONSTELLATION->sum_sigma_for_f107_average[INTEGRATOR->sc_main_nb][INTEGRATOR->sc_ensemble_nb] ,INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb); */
	  /* } */
	  /* 	} */
	  //	 	  printf("%d %d\n",index_in_driver_interpolated,INTEGRATOR->sc_main_nb);

 
	}


	else{

	  INTEGRATOR->Ap[index_in_driver_interpolated]            = OPTIONS->Ap[index_in_driver_interpolated];           // magnetic index(daily)

	  if (OPTIONS->use_ap_hist == 1){
	    for (hhh = 0; hhh < 7; hhh++){
	      INTEGRATOR->Ap_hist[hhh][index_in_driver_interpolated]            = OPTIONS->Ap_hist[hhh][index_in_driver_interpolated];           // magnetic index(historical)
	    }
	  }

	  INTEGRATOR->f107[index_in_driver_interpolated]          = OPTIONS->f107[index_in_driver_interpolated];         // Daily average of F10.7 flux
	  INTEGRATOR->f107A[index_in_driver_interpolated]         = OPTIONS->f107A[index_in_driver_interpolated];        // 81 day average of F10.7 flux
	  //	  printf("[%d] %f\n", index_in_driver_interpolated, OPTIONS->f107A[index_in_driver_interpolated]);

	}



	//  	    printf("%d %f %d %f\n", INTEGRATOR->sc_ensemble_nb, INTEGRATOR->f107[index_in_driver_interpolated], index_in_driver_interpolated, OPTIONS->f107[index_in_driver_interpolated]);
      } // end of if we don't run ensembles on F10.7/Ap or that we run ensemble on F10.7/Ap but that the current time is before the first prediction (so there is no uncertainty in F10.7/Ap because the time corresponds to an observation, not a prediction)
      //	  if (iProc == 0){
      //if (index_in_driver_interpolated == 4){ 
      // printf("eee: %d | f107:  %f | f107A: %f | Ap: %f | index: %d | iProc: %d\n", INTEGRATOR->sc_ensemble_nb, INTEGRATOR->f107[index_in_driver_interpolated],INTEGRATOR->f107A[index_in_driver_interpolated], INTEGRATOR->Ap[index_in_driver_interpolated] , index_in_driver_interpolated, iProc);
      // }
      //	  	  }
      /* if (iProc == 0){ */
      /*   //	    if ((INTEGRATOR->sc_ensemble_nb == 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + 1)){ */
      /*   if (index_in_driver_interpolated == 4){ */
      /*     printf("eee: %d | f107:  %f | index: %d | iProc: %d\n", INTEGRATOR->sc_ensemble_nb, INTEGRATOR->f107[index_in_driver_interpolated], index_in_driver_interpolated, iProc); */
      /* } */
      /* } */


      input.f107A = INTEGRATOR->f107A[index_in_driver_interpolated];
      input.f107  = INTEGRATOR->f107[index_in_driver_interpolated];

      //      printf(" %f | %f || ", input.f107, input.f107A);
      if (PARAMS->ATMOSPHERE.flags.switches[9] != -1){ // if daily ap
	if (INTEGRATOR->Ap[index_in_driver_interpolated] < 0){ // this can happen when running ensembles on Ap
	  INTEGRATOR->Ap[index_in_driver_interpolated] = 0;
	}

	input.ap    = INTEGRATOR->Ap[index_in_driver_interpolated];

/* 	if (INTEGRATOR->sc_ensemble_nb == 0){ */
/* 	etprint(et, "time"); */
/* 	printf("%d %d %f\n", INTEGRATOR->sc_main_nb, index_in_driver_interpolated,  INTEGRATOR->Ap[index_in_driver_interpolated]); */
/* 	} */

	if ( INTEGRATOR->sc_ensemble_nb == 0 ){
	  if (INTEGRATOR->last_compute_dxdt == 1){
	    if (iDebugLevel >= 5){
	      printf("------ (compute_drag) Writing in file_given_output...(iProc %d | iii = %d, eee = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb);
	    }
	    fprintf(INTEGRATOR->file_given_output, "%s: %f %f %f\n", timestamp_isoc, input.f107,input.f107A, input.ap);
	    if (iDebugLevel >= 5){
	      printf("------ (compute_drag) Done writing in file_given_output...(iProc %d | iii = %d, eee = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb);
	    }
	  }
	}

      }

      else{ //if historical ap

	input.ap_a = malloc(7 * sizeof(double));
	//		printf("%s: ", timestamp);
	if ( INTEGRATOR->sc_ensemble_nb == 0 ){
	  if (SC->INTEGRATOR.write_given_output == 1){
	    if (INTEGRATOR->last_compute_dxdt == 1){

	      fprintf(INTEGRATOR->file_given_output, "%s: %f %f", timestamp_isoc, input.f107,input.f107A);
		    
	      //  printf("%s %d %d %f\n", timestamp_isoc,SC->INTEGRATOR.sc_main_nb,index_in_driver_interpolated, input.f107,input.f107A);
	    }
	  }
	}
	//	printf("%s %d %d %d %f | ", timestamp_isoc,SC->INTEGRATOR.sc_main_nb, INTEGRATOR->sc_ensemble_nb,index_in_driver_interpolated, input.f107,input.f107A);
	//	printf("%d\n",  SC->INTEGRATOR.write_given_output );
	for (ppp = 0; ppp < 7 ; ppp ++){
	  if (INTEGRATOR->Ap_hist[ppp][index_in_driver_interpolated] < 0){ // this can happen when running ensembles on Ap
	    INTEGRATOR->Ap_hist[ppp][index_in_driver_interpolated] = 0;
	  }

	  input.ap_a->a[ppp]    = INTEGRATOR->Ap_hist[ppp][index_in_driver_interpolated];
	  if ( INTEGRATOR->sc_ensemble_nb == 0 ){
	    if (SC->INTEGRATOR.write_given_output == 1){
	      if (INTEGRATOR->last_compute_dxdt == 1){
		fprintf(INTEGRATOR->file_given_output, " %f", input.ap_a->a[ppp]);
	      }
	    }
	  }
	  //	   	  	  printf("[%d] %f | ", ppp, input.ap_a->a[ppp]);
		    
	}
	//	  printf("\n");
	if ( INTEGRATOR->sc_ensemble_nb == 0 ){
       
	  if (SC->INTEGRATOR.write_given_output == 1){
	    if (INTEGRATOR->last_compute_dxdt == 1){
	      fprintf(INTEGRATOR->file_given_output, "\n");
	    }
	  }
	}
      } // end of if historical ap    

      if (iDebugLevel >= 3){
	printf("---- (compute_drag) Done computing F10.7, F10.7A, and Ap to determine density for drag (iProc %d | iii = %d, eee = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb);
      }




    } // end of if the user chooses a time varying f107 and Ap for the density  
    //        printf("\n");
    //  printf("%f %f %f ||| %s\n", input.f107A, input.f107, input.ap, timestamp);
    /* if (input.f107A < 0){ */
    /*   print_test(); */
    /*   exit(0); */
    /* } */
    //  exit(0);


    /*     // !!!!!!!!!!!! remove block below */
    /*     double mean_earth_sun_distance = 149597870.700; */
    /*     double earth_sun_distance; */
    /*   double x[6]; */
    /*   double lt; */
    /*   double r_earth2sun_J2000[3]; */

    /*     spkez_c(10, et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. */

    /*     r_earth2sun_J2000[0] = x[0]; */
    /*     r_earth2sun_J2000[1] = x[1]; */
    /*     r_earth2sun_J2000[2] = x[2]; */
     
    /*     v_mag(&earth_sun_distance, r_earth2sun_J2000); */
    /* /\*     printf("%e %e\n", earth_sun_distance,mean_earth_sun_distance); *\/ */
    /* /\*     printf("b %f %f\n", input.f107, input.f107A ); *\/ */
    /*     input.f107 = input.f107 * ( mean_earth_sun_distance / earth_sun_distance ) * ( mean_earth_sun_distance / earth_sun_distance ); */
    /*     input.f107A = input.f107A * ( mean_earth_sun_distance / earth_sun_distance ) * ( mean_earth_sun_distance / earth_sun_distance ); */
    /* /\*     printf("a %f %f\n", input.f107, input.f107A ); *\/ */
    /*     // !!!!!!!!!!!! end of remove block below */

    // Call the NRL-MSIS-00 Atmosphere Model
    gtd7d(&input, &PARAMS->ATMOSPHERE.flags, &output);
    //    printf("T: %f\n", output.t[1]);
    if (PARAMS->ATMOSPHERE.flags.switches[9] == -1){ // if we used hitorical ap
      free(input.ap_a);
    }
    output.d[5] = (output.d[5] / 1000) * (100 * 100 * 100 ) *1e9; // Convert from Gram/cm^3 to kg/m^3 to kg/km^3 (CBV for the kg/km^3)
    *density = output.d[5];
    //    printf("%f %f %e\n", input.f107, input.f107A, *density);
    SC->INTEGRATOR.Ta = output.t[1];
    //    printf("%f %f\n", output.t[0], output.t[1]);
  } // end of the user chooses f107 and Ap for the density
  else if ( strcmp(INTEGRATOR->format_density_driver, "density_file") == 0 ){ // the user chooses to directly input the density from a file
    //    index_in_attitude_interpolated = floor( ( et - et_initial_epoch ) / ( INTEGRATOR->dt / 2.0) ) ;
    *density = INTEGRATOR->density[index_in_attitude_interpolated];
    SC->INTEGRATOR.Ta = 800.;
  } // end of the user chooses to directly input the density from a file

  else if ( strcmp(INTEGRATOR->format_density_driver, "gitm") == 0 ){ // the user chooses GITM
    // GITM outputs in longitude/latitude/altitude cordinates so the conversion ECI to longitude/latitude/altitude needs to be done before calling GITM
	
    // // ECI to ECEF
    /*     double geodetic[3]; */
    /* 	eci2lla(r_i2cg_INRTL , et, geodetic ); */
    /* 	altitude_gitm = geodetic[2]; */
    /* 	latitude_gitm = geodetic[0]; */
    /* 	longitude_gitm = geodetic[1];  */

    /* 	// update the planet fixed state		 */
    /* 	  geodetic_to_geocentric(PARAMS->EARTH.flattening,             */
    /* 				 altitude_gitm, */
    /* 				 latitude_gitm, */
    /* 				 longitude_gitm, */
    /* 				 PARAMS->EARTH.radius,        */
    /* 				 r_ecef2cg_ECEF_gitm) ; */

    pxform_c("J2000", PARAMS->EARTH.earth_fixed_frame, et, T_J2000_to_ECEF_gitm);
    m_x_v(r_ecef2cg_ECEF_gitm, T_J2000_to_ECEF_gitm, r_i2cg_INRTL);
    
    // // ECEF to longitude/latitude/altitude
    geocentric_to_geodetic( r_ecef2cg_ECEF_gitm,
			    &PARAMS->EARTH.radius,
			    &PARAMS->EARTH.flattening,
			    &altitude_gitm,
			    &latitude_gitm,
			    &longitude_gitm);

    gitm_density(density, et, altitude_gitm, latitude_gitm, longitude_gitm, PARAMS);
    *density = *density * 1e9;
    SC->INTEGRATOR.Ta = 800.;
  } // end of the user chooses GITM
  //				printf("density = %e\n", *density);
  /////////////////// END OF COMPUTE DENSITY ///////////////////

  //  APPLY DENSITY_MOD: FACTOR TO APPLY ON DENSITY AT POSITION OF SATELLITE (CALCULATED BY NRLMSIS, GITM OR FROM DENSITY FILE)

  // 

  // // LATITUDE DEPENDENT: MAX OF SINE AT EQUATOR IF INTEGRATOR.DENSITY_MOD_PHASE = 0
  //  double fac_density =  SC->INTEGRATOR.density_mod + SC->INTEGRATOR.density_mod_amp *  sin(SC->GEODETIC.latitude * M_PI/SC->OE.inclination + (SC->INTEGRATOR.density_mod_phase+0.25)*2*M_PI); // 0.25 so that if phase = 0 then max of sine at the equator and min of sine at "poles" (max/min of lat)
  // // end of LATITUDE DEPENDENT: MAX OF SINE AT EQUATOR IF INTEGRATOR.DENSITY_MOD_PHASE = 0

/*   // // PHASE (ARG PER + TRUE ANO) DEPENDENT: MAX OF SINE AT PERIGEE IF INTEGRATOR.DENSITY_MOD_PHASE = 0 */
/*   double fac_density; */
/*   if ( SC->orbit_number >= 1){ // get in here only if at least one full orbit has been traveled (otherwise orbit average arg per is not defined) */
/*     cart2kep(&SC->OE, r_i2cg_INRTL, v_i2cg_INRTL, et, PARAMS->EARTH.GRAVITY.mu); */
/*     double w_ave_to_sc = fmod(SC->OE.an_to_sc - SC->OE.w_ave, 2*M_PI); */
/*     //    printf("%f %f %f\n", an_to_sc*180/M_PI, SC->OE.w_ave*180/M_PI, w_ave_to_sc*180/M_PI); */
/*     fac_density =  SC->INTEGRATOR.density_mod + SC->INTEGRATOR.density_mod_amp *  sin(w_ave_to_sc + (SC->INTEGRATOR.density_mod_phase+0.25)*2*M_PI); */
/*   } */
/*   else{ */
/*    fac_density = SC->INTEGRATOR.density_mod; */
/*   } */
/*   // // end of PHASE (ARG PER + TRUE ANO) DEPENDENT: MAX OF SINE AT PERIGEE IF INTEGRATOR.DENSITY_MOD_PHASE = 0 */

/* //  // PHASE (ARG PER + TRUE ANO) DEPENDENT: MAX OF SINE AT LOCAL TIME OF 0 DEG ("MIDNIGHT") IF INTEGRATOR.DENSITY_MOD_PHASE = 0 */
/*   int hr, mn, sc; */
/*   char time_local[256], ampm[256]; */
/*   double lon_an_dn = 0 ; */

/*     et2lst_c ( SC->et,  399,  SC->GEODETIC.longitude, "PLANETOCENTRIC", 51, 51, */
/* 	       &hr, &mn,  &sc,  time_local, ampm             ); */

/*     lon_an_dn = fmod( ( hr + mn / 60.0 + sc / 3600.0 ) * 15 , 360) * M_PI / 180.; */
/*     double fac_density; */
/*     fac_density =  SC->INTEGRATOR.density_mod + SC->INTEGRATOR.density_mod_amp *  sin(lon_an_dn + (SC->INTEGRATOR.density_mod_phase+0.25)*2*M_PI); */


/* //  // end of PHASE (ARG PER + TRUE ANO) DEPENDENT: MAX OF SINE AT LOCAL TIME OF 0 DEG ("MIDNIGHT") IF INTEGRATOR.DENSITY_MOD_PHASE = 0 */


//  // PHASE (ARG PER + TRUE ANO) DEPENDENT: MAX OF SINE AT SOLAR ZENITH OF 0 DEG IF INTEGRATOR.DENSITY_MOD_PHASE = 0


  /* double x[6]; */
  /* double lt; */
  /* double r_earth2sun_J2000[3]; */
  /* double r_cg2sun_J2000[3]; */
  /* double r_cg2sun_J2000_normalized[3]; */
  /*   spkez_c(10, et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. */

  /*   r_earth2sun_J2000[0] = x[0]; */
  /*   r_earth2sun_J2000[1] = x[1]; */
  /*   r_earth2sun_J2000[2] = x[2]; */
  /*   double r_i2cg_INRTL_normalized[3]; */
  /*   v_sub(r_cg2sun_J2000, r_earth2sun_J2000, r_i2cg_INRTL); */
  /*   v_norm(r_cg2sun_J2000_normalized, r_cg2sun_J2000); */
  /*   v_norm(r_i2cg_INRTL_normalized, r_i2cg_INRTL); */
  /*   double cos_earthtosc_sctosun; */
  /*   v_dot(&cos_earthtosc_sctosun, r_i2cg_INRTL_normalized, r_cg2sun_J2000_normalized); */

  /*   double cos_v_sctosun; */
  /*   double v_i2cg_INRTL_normalized[3]; */
  /*   v_norm(v_i2cg_INRTL_normalized, v_i2cg_INRTL); */
  /*   v_dot(&cos_v_sctosun, v_i2cg_INRTL_normalized, r_cg2sun_J2000_normalized); */
  /*   if (cos_v_sctosun > 0){ */
  /*   SC->OE.zenith = 2*M_PI - acos(cos_earthtosc_sctosun); */
  /*   } */
  /*   else{ */
  /*   SC->OE.zenith = acos(cos_earthtosc_sctosun); */
  /*   } */
    double fac_density;
    fac_density =  SC->INTEGRATOR.density_mod + SC->INTEGRATOR.density_mod_amp *  sin(SC->OE.zenith + (SC->INTEGRATOR.density_mod_phase+0.25)*2*M_PI);


//  // end of PHASE (ARG PER + TRUE ANO) DEPENDENT: MAX OF SINE AT LOCAL TIME OF 0 DEG ("MIDNIGHT") IF INTEGRATOR.DENSITY_MOD_PHASE = 0
  
  

/*   etprint(et, ""); */
/*   printf("factor %f | %f %f\n", fac_density, SC->OE.an_to_sc*180/M_PI, SC->OE.w_ave*180/M_PI); */

  *density = *density * fac_density;
  //    	              *density = fac_density ;  //!!!!!!!!!!!!! remove line and uncomment right above
  //  printf("%f\n", SC->INTEGRATOR.density_mod + SC->INTEGRATOR.density_mod_amp * sin(2*M_PI*et/SC->OE.period + SC->INTEGRATOR.density_mod_phase));

  compute_T_inrtl_2_ntw(T_inrtl_2_ntw, r_i2cg_INRTL, v_i2cg_INRTL);

  // !!!!!!!!!!!!!!!! the code for the rotating atmo has not been checked yet
  // Take into account the rotating atmosphere: v_rel = v_sat - v_rot_atmo
  /*   //  // Update ECEF state */
  /*   SpiceDouble       xform[6][6]; */
  /*   double estate[6], jstate[6]; */

  /*   estate[0] = r_i2cg_INRTL[0];estate[1] = r_i2cg_INRTL[1];estate[2] = r_i2cg_INRTL[2]; */
  /*   estate[3] = v_i2cg_INRTL[0];estate[4] = v_i2cg_INRTL[1];estate[5] = v_i2cg_INRTL[2]; */
  /*   sxform_c (  "J2000", PARAMS->EARTH.earth_fixed_frame,  et,    xform  );  */
  /*   mxvg_c   (  xform,       estate,   6,  6, jstate );  */
  /*   double r_ecef2cg_ECEF_for_rot_atmo[3], v_ecef2cg_ECEF_for_rot_atmo[3]; */
  /*   r_ecef2cg_ECEF_for_rot_atmo[0] = jstate[0]; r_ecef2cg_ECEF_for_rot_atmo[1] = jstate[1]; r_ecef2cg_ECEF_for_rot_atmo[2] = jstate[2]; */
  /*   v_ecef2cg_ECEF_for_rot_atmo[0] = jstate[3]; v_ecef2cg_ECEF_for_rot_atmo[1] = jstate[4]; v_ecef2cg_ECEF_for_rot_atmo[2] = jstate[5]; */

  /*   double v_rel[3]; */
  /*   double v_rot_atmo[3]; */
  /*   double omega_rot_atmo[3]; */
  /*   omega_rot_atmo[0] = 0; omega_rot_atmo[1] = 0; omega_rot_atmo[2] =  0.0;// !!!!!! should be 0.00007292158553; */
  /*   v_cross(v_rot_atmo, omega_rot_atmo, r_ecef2cg_ECEF_for_rot_atmo); */
  /*   v_sub(v_rel, v_ecef2cg_ECEF_for_rot_atmo, v_rot_atmo); */
  
  /*   // // Back to ECI */
  /*   SpiceDouble       xform_new[6][6]; */
  /*   double estate_new[6], jstate_new[6]; */
  /*   sxform_c (  PARAMS->EARTH.earth_fixed_frame, "J2000",  et,    xform_new  );  */
  /*   estate_new[0] = r_ecef2cg_ECEF_for_rot_atmo[0];estate_new[1] = r_ecef2cg_ECEF_for_rot_atmo[1];estate_new[2] = r_ecef2cg_ECEF_for_rot_atmo[2]; */
  /*   estate_new[3] = v_rel[0];estate_new[4] = v_rel[1];estate_new[5] = v_rel[2]; */
  /*   mxvg_c   (  xform_new,       estate_new,   6,  6, jstate_new );  */
  /*   double v_rel_eci[3]; */
  /*   v_rel_eci[0] = jstate_new[3]; v_rel_eci[1] = jstate_new[4]; v_rel_eci[2] = jstate_new[5]; */
  /*   /\* v_print(v_rot_atmo, "v_rot_atmo"); *\/ */
  /*   /\* v_print(v_rel, "v_rel"); *\/ */
  double v_rel_eci[3];
  double v_rot_atmo[3];
  double omega_rot_atmo_scal = 0.00007292158553;
  v_rel_eci[0] = v_i2cg_INRTL[0] + omega_rot_atmo_scal * r_i2cg_INRTL[1];
  v_rel_eci[1] = v_i2cg_INRTL[1] - omega_rot_atmo_scal * r_i2cg_INRTL[0];
  v_rel_eci[2] = v_i2cg_INRTL[2];
  // End of take into account the rotating atmosphere: v_rel = v_sat - v_rot_atmo


    m_trans(T_ntw_2_inrtl, T_inrtl_2_ntw);
  m_x_v(v_i2cg_NTW, T_inrtl_2_ntw, v_rel_eci); //  !!!!!!!!!!!!!!!! v_rel_eci here is new!Before it was: v_i2cg_INRTL
  // !!!!!!!!!!!!!!!! the code for the rotating atmo has not been checked yet
  //  v_print(v_i2cg_NTW, "v_i2cg_NTW");

  if (INTEGRATOR->coll_vcm != 1){
  if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") == 0){
    v_angle[0] = INTEGRATOR->attitude.pitch_angular_velocity_ensemble * ( et - et_sc_initial ) + INTEGRATOR->attitude.pitch_ini_ensemble;
    v_angle[1] = INTEGRATOR->attitude.roll_angular_velocity_ensemble * ( et - et_sc_initial ) + INTEGRATOR->attitude.roll_ini_ensemble;
    v_angle[2] = INTEGRATOR->attitude.yaw_angular_velocity_ensemble * ( et - et_sc_initial ) + INTEGRATOR->attitude.yaw_ini_ensemble;
    order_rotation[0] = 1; // !!!!!!!! we might want to change that in the future
    order_rotation[1] = 2;
    order_rotation[2] = 3;
    INTEGRATOR->attitude.pitch_current = v_angle[0];
    INTEGRATOR->attitude.roll_current = v_angle[1];
    INTEGRATOR->attitude.yaw_current = v_angle[2];
    INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];

  }

  if (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") == 0) {
    v_angle[0] =  INTEGRATOR->attitude.pitch_for_attitude_ensemble  +  INTEGRATOR->attitude.pitch_angular_velocity_constant * ( et - et_sc_initial );
    v_angle[1] =  INTEGRATOR->attitude.roll_for_attitude_ensemble  +  INTEGRATOR->attitude.roll_angular_velocity_constant * ( et - et_sc_initial );
    v_angle[2] =  INTEGRATOR->attitude.yaw_for_attitude_ensemble  +  INTEGRATOR->attitude.yaw_angular_velocity_constant * ( et - et_sc_initial );

    order_rotation[0]  = 1; order_rotation[1]  = 2; order_rotation[2]  = 3;
	   
    INTEGRATOR->attitude.pitch_current = v_angle[0];
    INTEGRATOR->attitude.roll_current = v_angle[1];
    INTEGRATOR->attitude.yaw_current = v_angle[2];
    INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
    INTEGRATOR->attitude.order_roll_current = order_rotation[1];
    INTEGRATOR->attitude.order_yaw_current = order_rotation[2];


  }
  
  if ( (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_angular_velocity") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "ensemble_initial_attitude") != 0) && (strcmp(INTEGRATOR->attitude.attitude_profile, "sun_pointed") != 0) ){ // otherwise (atittude is nadir, sun_pointed or manual (from an input file)) 
    //  index_in_attitude_interpolated = floor( ( et - et_initial_epoch ) / ( INTEGRATOR->dt / 2.0) ) ;


    //    et2utc_c(et, "D", 3, 35, timestamp);
    //  printf("et = %s ||", timestamp);
    /*  et2utc_c(et_initial_epoch, "D", 3, 35, timestamp); */
    /* printf("et_initial_epoch = %s \n", timestamp); */
    if (INTEGRATOR->file_is_quaternion == 0){
      v_angle[0] = INTEGRATOR->attitude.pitch[index_in_attitude_interpolated];
      v_angle[1] = INTEGRATOR->attitude.roll[index_in_attitude_interpolated];
      v_angle[2] = INTEGRATOR->attitude.yaw[index_in_attitude_interpolated];
      //    printf("%d %d %d\n", INTEGRATOR->attitude.order_pitch[index_in_attitude_interpolated], index_in_attitude_interpolated, INTEGRATOR->sc_ensemble_nb);
      order_rotation[0] = INTEGRATOR->attitude.order_pitch[index_in_attitude_interpolated];
      order_rotation[1] = INTEGRATOR->attitude.order_roll[index_in_attitude_interpolated];
      order_rotation[2] = INTEGRATOR->attitude.order_yaw[index_in_attitude_interpolated];
      // printf("%d %d %d\n", order_rotation[0], order_rotation[1],order_rotation[2]);
      INTEGRATOR->attitude.pitch_current = v_angle[0];
      INTEGRATOR->attitude.roll_current = v_angle[1];
      INTEGRATOR->attitude.yaw_current = v_angle[2];
      INTEGRATOR->attitude.order_pitch_current = order_rotation[0];
      INTEGRATOR->attitude.order_roll_current = order_rotation[1];
      INTEGRATOR->attitude.order_yaw_current = order_rotation[2];

    }
    else{
      q_copy( INTEGRATOR->attitude.quaternion_current, INTEGRATOR->attitude.quaternion[index_in_attitude_interpolated]);
    }
    
    //       printf("%d | %e | %e | %e \n", index_in_attitude_interpolated, v_angle[0], v_angle[1],v_angle[2]);
  }

  compute_T_sc_to_lvlh( T_sc_to_lvlh, v_angle, order_rotation, INTEGRATOR->attitude.attitude_profile, &et,  r_i2cg_INRTL, v_i2cg_INRTL, INTEGRATOR->file_is_quaternion, INTEGRATOR->attitude.quaternion_current, PARAMS);



  //  m_print(T_sc_to_lvlh, "T_sc_to_lvlh");
  //  printf("%d %d %d %d %d\n", index_in_attitude_interpolated, OPTIONS->nb_time_steps*2, order_rotation[0], order_rotation[1], order_rotation[2]);

  /* LVLH to inertial */

  compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, r_i2cg_INRTL, v_i2cg_INRTL);
  m_trans(T_lvlh_2_inrtl, T_inrtl_2_lvlh);

  double v_i2cg_NTW_mag;
  v_mag(&v_i2cg_NTW_mag, v_i2cg_NTW);

  // openGL computation of cross-section area
  if (OPTIONS->opengl == 1){ 
    // // compute relative speed satellite in body frame
    double v_rel_lvlh[3];
    m_x_v(v_rel_lvlh, T_inrtl_2_lvlh, v_rel_eci); 
    double v_rel_body[3];
    double T_lvlh_to_sc[3][3];
    m_trans(T_lvlh_to_sc, T_sc_to_lvlh);
    m_x_v(v_rel_body, T_lvlh_to_sc, v_rel_lvlh); 
    double v_rel_body_normalized[3];
    v_norm(v_rel_body_normalized, v_rel_body);
    // // compute the azimuth and elevation of the relative speed vector in body frame. 
    // // Actually the spherical coordinates are used. 
    // // theta is the angle between z_body and the relative speed vector (from 0 to -180 deg). 
    // // phi is the angle between x_body and the projection of the relative speed vector on the (x_body, y_body) plane (from 0 to 360 deg).
    double x_body[3], y_body[3], z_body[3];
    x_body[0] = 1; x_body[1] = 0; x_body[2] = 0;
    y_body[0] = 0; y_body[1] = 1; y_body[2] = 0;
    z_body[0] = 0; z_body[1] = 0; z_body[2] = 1;
    double theta, phi;
    // // // theta
    theta = acos(v_rel_body_normalized[2]); // v_rel_body_normalized[2] is euqal to the projection of v_rel_body_normalized on z_body
    // // // phi
    phi = acos(v_rel_body_normalized[0]);
    if (v_rel_body_normalized[1] < 0){
      phi = 2*M_PI - phi;
    }


    double theta_deg, phi_deg;
    theta_deg = theta * RAD2DEG;
    phi_deg = phi * RAD2DEG;
    // // get the cross-section area from the file generated by openGL
    // // // compute the indices itheta and iphi
    int itheta, iphi, iface;
    A_ref = 0;
/*     itheta = (int)((theta_deg - CONSTELLATION->area_attitude_opengl_theta0) / CONSTELLATION->area_attitude_opengl_dtheta); */
/*     iphi = (int)((phi_deg - CONSTELLATION->area_attitude_opengl_phi0) / CONSTELLATION->area_attitude_opengl_dphi); */
    itheta = round((theta_deg - CONSTELLATION->area_attitude_opengl_theta0) / CONSTELLATION->area_attitude_opengl_dtheta);
    iphi = round((phi_deg - CONSTELLATION->area_attitude_opengl_phi0) / CONSTELLATION->area_attitude_opengl_dphi);



    // // compute the drag acceleration. 
    // when using a 3d model (opengl option), all surfaces are assigned the same cd (if OPTIONS->new_cd != 1) or the same accomodation coefficient (if OPTIONS->new_cd == 1). This is because as of now, we don’t know how to assign a coefficient (e.g, cd or alpha) to a triangle

      if (OPTIONS->new_cd !=1){
	cd = INTEGRATOR->surface[0].Cd; // take cd of first surface because when using a 3d model (opengl option), all surfaces are assigned the same cd. This is because as of now, we don’t know how to assign a coefficient (e.g, cd or alpha) to a triangle 
	A_ref = CONSTELLATION->area_attitude_opengl_total[itheta][iphi];
	ballistic_coefficient =  cd * A_ref / INTEGRATOR->mass; 
	A_ref_tot =  A_ref * 1000000; //A_ref_tot in m^2
	cd_tot_norm = cd*A_ref * 1000000;
	

	/* if (index_in_attitude_interpolated >= 20){ */
	/* printf("%f %f\n", theta_deg, phi_deg); */
	/* printf("theta %.2f, phi %.2f (%.2f, %.2f)\n", itheta*CONSTELLATION->area_attitude_opengl_dtheta + CONSTELLATION->area_attitude_opengl_theta0, iphi*CONSTELLATION->area_attitude_opengl_dphi + CONSTELLATION->area_attitude_opengl_phi0, theta_deg, phi_deg); */
	/* printf("cd: %f | A_ref %f cm2 | bc %e\n", cd, A_ref * 1e10, ballistic_coefficient * INTEGRATOR->mass) ; */
	/* } */

      }


      else{
	for (iface = 0; iface < CONSTELLATION->nb_faces_theta_phi[itheta][iphi]; iface++){ 
	  iface_here = CONSTELLATION->which_face_theta_phi[itheta][iphi][iface];
	m_x_v(normal_in_lvlh, T_sc_to_lvlh,  CONSTELLATION->normal_face[iface_here]);
	/* LVLH to inertial */
	m_x_v(normal_in_inertial, T_lvlh_2_inrtl, normal_in_lvlh);
	/* inertial to NTW (note: we could simply make the dot product of the speed and the normal vector of the surface in the intertial frame, we don't have to do convert both vectors in the NTW frame. But this was done because it was easy to check that the consistency of the results in the NTW frame (it does not add that much time to the computation)) */
	m_x_v(normal_in_ntw, T_inrtl_2_ntw, normal_in_inertial);
	v_dot(&v_i2cg_NTW_DOT_normal_to_the_surface_NTW, v_i2cg_NTW, normal_in_ntw);

	v_mag(&normal_in_ntw_normalized, normal_in_ntw);// just in case it's not normalized 
	cos_angle_v_sc_normal = v_i2cg_NTW_DOT_normal_to_the_surface_NTW / ( v_i2cg_NTW_mag * normal_in_ntw_normalized );
		if (cos_angle_v_sc_normal > 5e-3){ // 1e-3 for numerical reasons (should be 0 theoritically). Note that in theory opengl touputs only faces that are seen so with cos > 0. But there is a smoall bug in the code that ouptuts a very small area of a surface with cos < 0. So this if here excludes this surface
	  calculate_cd_opengl(&cd,
		       OPTIONS->surface[0].acco_coeff, //when using a 3d model (opengl option), all surfaces are assigned the  same accomodation coefficient (if OPTIONS->new_cd == 1). This is because as of now, we don’t know how to assign a coefficient (here, the accomodation coefficient) to a triangle
		       v_i2cg_NTW_mag, //in km/s
		       SC->INTEGRATOR.Ta, // atmospheric temperature in K, from NRLSMSIS if the user didn't chose "density_file" or "gitm", in which case SC->INTEGRATOR.Ta = 800K
		       CONSTELLATION->area_attitude_opengl[itheta][iphi][iface_here], // in m^2
		       cos_angle_v_sc_normal,
		       16./1000 // in kg/mol!!!!! assuming mainly atom of O
		       );
	  //	  v_print(CONSTELLATION->normal_face[iface_here], "n");
	  
	  A_ref = CONSTELLATION->area_attitude_opengl[itheta][iphi][iface_here];
  /* 	if (index_in_attitude_interpolated >= 20){ */
  /* printf("%f %f\n", theta_deg, phi_deg); */
  /* printf("theta %.2f, phi %.2f | %d faces\n", itheta*CONSTELLATION->area_attitude_opengl_dtheta + CONSTELLATION->area_attitude_opengl_theta0, iphi*CONSTELLATION->area_attitude_opengl_dphi + CONSTELLATION->area_attitude_opengl_phi0,  CONSTELLATION->nb_faces_theta_phi[itheta][iphi]); */
  /* 	  printf("cd[%d]: %f | A_ref %f cm2 | cd * A_ref %e | cos %f  \n", iface_here, cd,  A_ref * 1e10, cd*A_ref* 1000000, cos_angle_v_sc_normal) ; */
  /* 	} */
/* 	  if (iface_here == 3){ */
/* 	    exitf(); */
/*       } */
	  ballistic_coefficient = ballistic_coefficient + cd * A_ref / INTEGRATOR->mass;
	  A_ref_tot = A_ref_tot + A_ref * 1000000; //A_ref_tot in m^2    
	  cd_tot_norm = cd_tot_norm + cd*A_ref * 1000000;
	  	}
	}
	//	printf("cd_tot_norm %f\n", cd_tot_norm / A_ref_tot);
      }
/*     cd = INTEGRATOR->surface[0].Cd; */
/*     ballistic_coefficient = cd * A_ref / INTEGRATOR->mass; */
      
/*       printf("%f %f\n", theta_deg, phi_deg); */
/*       printf("theta %.2f, phi %.2f\n", itheta*CONSTELLATION->area_attitude_opengl_dtheta + CONSTELLATION->area_attitude_opengl_theta0, iphi*CONSTELLATION->area_attitude_opengl_dphi + CONSTELLATION->area_attitude_opengl_phi0);       */
/*       printf("BC %e\n" , ballistic_coefficient * INTEGRATOR->mass); */
    a_i2cg_NTW[0] = -0.5 * ballistic_coefficient * *density * v_i2cg_NTW_mag * v_i2cg_NTW[0];
    a_i2cg_NTW[1] = -0.5 * ballistic_coefficient * *density * v_i2cg_NTW_mag * v_i2cg_NTW[1];
    a_i2cg_NTW[2] = -0.5 * ballistic_coefficient * *density * v_i2cg_NTW_mag * v_i2cg_NTW[2];
    
  } // end of if (OPTIONS->opengl == 1)
  // end of openGL computation of cross-section area
  else{


  //	   	etprint(et, "");

  if (INTEGRATOR->initialize_geo_with_bstar != 1){

    // to not use the new cd calculation, uncomment line below
    if (OPTIONS->new_cd !=1){
      ballistic_coefficient = INTEGRATOR->surface[0].Cd * INTEGRATOR->surface[0].area / INTEGRATOR->mass;
      cd =  INTEGRATOR->surface[0].Cd;
    }
    /* 	if (SC->INTEGRATOR.write_given_output == 1){ */
    /* 	  	  fprintf(INTEGRATOR->file_given_output, "%s ", timestamp_isoc ); */
    /* 	} */
    for (sss = 0; sss < INTEGRATOR->nb_surfaces; sss++){
      //      printf("surface %d\n", sss);
      if (sss == 0){
	/* Conversion of the normal vector from the SC reference system to the NTW reference system */

	/* SC to LVLH */
	m_x_v(normal_in_lvlh, T_sc_to_lvlh, INTEGRATOR->surface[0].normal );
	/* LVLH to inertial */
	m_x_v(normal_in_inertial, T_lvlh_2_inrtl, normal_in_lvlh);



	/* inertial to NTW (note: we could simply make the dot product of the speed and the normal vector of the surface in the intertial frame, we don't have to do convert both vectors in the NTW frame! But this was done because it was easy to check that the consistency of the results in the NTW frame (it does not add that much time at the computation)) */

	m_x_v(normal_in_ntw, T_inrtl_2_ntw, normal_in_inertial);
	v_dot(&v_i2cg_NTW_DOT_normal_to_the_surface_NTW, v_i2cg_NTW, normal_in_ntw);



	// !!!!!!!!!!!!!!!!!!!!!! REMOVE LINES BELOW!!!!!!!!!

	/* v_mag(&v_i2cg_NTW_DOT_normal_to_the_surface_NTW, v_i2cg_NTW); */
	/* ballistic_coefficient = 2.2 * 0.01 / 1000./1000.; */
	// !!!!!!!!!!!!!!!!!!!!!! END OF REMOVE LINES BELOW!!!!!!!!!	

	v_mag(&normal_in_ntw_normalized, normal_in_ntw);// just in case it's not normalized 
	cos_angle_v_sc_normal = v_i2cg_NTW_DOT_normal_to_the_surface_NTW / ( v_i2cg_NTW_mag * normal_in_ntw_normalized );
	//      printf("%e\n", cos_angle_v_sc_normal);

	if (cos_angle_v_sc_normal > 5e-3){ // 1e-3 for numerical reasons (should be 0 theoritically)

	  
	  // new cd
	  if (OPTIONS->new_cd == 1){

	    calculate_cd(&cd,
			 OPTIONS->surface[0].acco_coeff,// !!!!!!!! shoulc be OPTIONS->surface[0].acco_coeff, 
			 v_i2cg_NTW_mag, //in km/s
			 SC->INTEGRATOR.Ta, // atmospheric temperature in K, from NRLSMSIS if the user didn't chose "density_file" or "gitm", in which case SC->INTEGRATOR.Ta = 800K
			 INTEGRATOR->surface[0].area, // in m^2
			 cos_angle_v_sc_normal,
			 16./1000 // in kg/mol!!!!! assuming mainly atom of O
			 );
	    //	    	printf("cd[%d]: %f | %f %f\n", sss, cd,INTEGRATOR->surface[sss].area*cos_angle_v_sc_normal*1e10, OPTIONS->surface[0].acco_coeff );
	    ballistic_coefficient = cd * INTEGRATOR->surface[0].area / INTEGRATOR->mass;
	  }
	  // end of new cd
	  A_ref = INTEGRATOR->surface[0].area * cos_angle_v_sc_normal * 1000000; // convert km^2 to m^2

	  total_cross_area = total_cross_area + INTEGRATOR->surface[0].area* 1000000 * cos_angle_v_sc_normal;
	  /* 	if (SC->INTEGRATOR.write_given_output == 1){ */
	  /* 	  	  fprintf(INTEGRATOR->file_given_output, "%e %d %e (%e) %e | %e || ", cd, sss, cos_angle_v_sc_normal, v_i2cg_NTW_DOT_normal_to_the_surface_NTW , A_ref, cd*A_ref); */
	  /* 	} */

	  A_ref_tot = A_ref_tot + A_ref;
	  cd_tot_norm = cd_tot_norm + cd*A_ref;


	  a_i2cg_NTW[0] =  -0.5 * *density * v_i2cg_NTW_DOT_normal_to_the_surface_NTW * v_i2cg_NTW[0] * ballistic_coefficient; // Correction by CBV: added the factor INTEGRATOR->Cd //  first component of NTW is v normalized-> IN TRACK (CBV)
	  a_i2cg_NTW[1] =  -0.5 * *density * v_i2cg_NTW_DOT_normal_to_the_surface_NTW * v_i2cg_NTW[1] * ballistic_coefficient; // radial cross v (then normalized) -> CROSS-TRACK (CBV)
	  a_i2cg_NTW[2] =  -0.5 * *density * v_i2cg_NTW_DOT_normal_to_the_surface_NTW * v_i2cg_NTW[2] * ballistic_coefficient; // first component cross second component (CBV)

	}
	else{
	  a_i2cg_NTW[0] = 0.0;
	  a_i2cg_NTW[1] = 0.0;
	  a_i2cg_NTW[2] = 0.0;
	}
      } // end of if sss == 0
      else{

	/* Conversion of the normal vector from the SC reference system to the NTW reference system */

	/* SC to LVLH */
	m_x_v(normal_in_lvlh, T_sc_to_lvlh, INTEGRATOR->surface[sss].normal );
      
	/* LVLH to inertial */
	m_x_v(normal_in_inertial, T_lvlh_2_inrtl, normal_in_lvlh);


	/* inertial to NTW (note: we could simply make the dot product of the speed and the normal vector of the surface in the intertial frame, we don't have to do convert both vectors in the NTW frame! But this was done because it was easy to check that the consistency of the results in the NTW frame (it does not add that much time at the computation)) */
	m_x_v(normal_in_ntw, T_inrtl_2_ntw, normal_in_inertial);
	
	v_dot(&v_i2cg_NTW_DOT_normal_to_the_surface_NTW, v_i2cg_NTW, normal_in_ntw);
	v_mag(&normal_in_ntw_normalized, normal_in_ntw);// just in case it's not normalized 
	cos_angle_v_sc_normal = v_i2cg_NTW_DOT_normal_to_the_surface_NTW / ( v_i2cg_NTW_mag * normal_in_ntw_normalized );


	if (cos_angle_v_sc_normal > 5e-3){ // 1e-3 for numerical reasons (should be 0 theoritically)
	  if (INTEGRATOR->initialize_geo_with_bstar != 1){


	    // new cd
	    if (OPTIONS->new_cd == 1){
	      calculate_cd(&cd,
			   OPTIONS->surface[sss].acco_coeff,// !!!!!!!! shoul be OPTIONS->surface[sss].acco_coeff, 
			   v_i2cg_NTW_mag, //in km/s
			   SC->INTEGRATOR.Ta, // atmospheric temperature in K, from NRLSMSIS if the user didn't chose "density_file" or "gitm", in which case SC->INTEGRATOR.Ta = 800K
			   INTEGRATOR->surface[sss].area, // in m^2
			   cos_angle_v_sc_normal,
			   16./1000 // in kg/mol!!!!! assuming mainly atom of O
			   );
	      //	      	printf("cd[%d]: %f | %f %f\n", sss, cd,INTEGRATOR->surface[sss].area*cos_angle_v_sc_normal*1e10, OPTIONS->surface[sss].acco_coeff );
	      ballistic_coefficient = cd * INTEGRATOR->surface[sss].area / INTEGRATOR->mass;
	    }
	    // end of new cd
	    A_ref = INTEGRATOR->surface[sss].area * cos_angle_v_sc_normal * 1000000; // convert km^2 to m^2

	    /* 	if (SC->INTEGRATOR.write_given_output == 1){ */
	    /* 	  	  total_cross_area = total_cross_area + INTEGRATOR->surface[sss].area* 1000000 * cos_angle_v_sc_normal; */
	    /* 		  //		  fprintf(INTEGRATOR->file_given_output, "%e %d %e (%e) %e | %e || ", cd, sss, cos_angle_v_sc_normal,v_i2cg_NTW_DOT_normal_to_the_surface_NTW, A_ref,  cd*A_ref ); */
	    /* 	} */



	    if (OPTIONS->new_cd !=1){
	      ballistic_coefficient = INTEGRATOR->surface[sss].Cd * INTEGRATOR->surface[sss].area / INTEGRATOR->mass;
	      cd = INTEGRATOR->surface[sss].Cd;
	/* if (index_in_attitude_interpolated >= 20){ */
	/*       	printf("cd[%d]: %f | %f %f\n", sss, cd,INTEGRATOR->surface[sss].area*cos_angle_v_sc_normal*1e10, OPTIONS->surface[sss].acco_coeff ); */
	/* } */
	    }

	    A_ref_tot = A_ref_tot + A_ref;
	    cd_tot_norm = cd_tot_norm + cd*A_ref;

	  }


	  a_i2cg_NTW[0] = a_i2cg_NTW[0]  -0.5 * *density * v_i2cg_NTW_DOT_normal_to_the_surface_NTW * v_i2cg_NTW[0] * ballistic_coefficient; // Correction by CBV: added the factor INTEGRATOR->Cd //  first component of NTW is v normalized-> IN TRACK (CBV)
	  a_i2cg_NTW[1] = a_i2cg_NTW[1] -0.5 * *density * v_i2cg_NTW_DOT_normal_to_the_surface_NTW * v_i2cg_NTW[1] * ballistic_coefficient; // radial cross v (then normalized) -> CROSS-TRACK (CBV)
	  a_i2cg_NTW[2] = a_i2cg_NTW[2] -0.5 * *density * v_i2cg_NTW_DOT_normal_to_the_surface_NTW * v_i2cg_NTW[2] * ballistic_coefficient; // first component cross second component (CBV)

	}
	
      }
   
    }
  }

  else if (INTEGRATOR->initialize_geo_with_bstar == 1){

    ballistic_coefficient = 6378. * 2.461 * 0.00001 / (2 * INTEGRATOR->bstar ) ;
    ballistic_coefficient = 1/ballistic_coefficient / 1000000; // "/1000000" to convert from m^2/kg to km^2/kg

    a_i2cg_NTW[0] =  -0.5 * *density *  v_i2cg_NTW[0] * v_i2cg_NTW[0] * ballistic_coefficient; // Correction by CBV: added the factor INTEGRATOR->Cd //  first component of NTW is v normalized-> IN TRACK (CBV)
    a_i2cg_NTW[1] =  -0.5 * *density * v_i2cg_NTW[1] * v_i2cg_NTW[1] * ballistic_coefficient; // radial cross v (then normalized) -> CROSS-TRACK (CBV)
    a_i2cg_NTW[2] =  -0.5 * *density * v_i2cg_NTW[2] * v_i2cg_NTW[2] * ballistic_coefficient; // first component cross second component (CBV)

  }
  } // end of openGL != 1


  // !!!!!!!!!!!!!!!!!! ERASE LINES BELOW
  /*  v_norm_print(v_i2cg_NTW, "v_i2cg_NTW"); */
  /*  v_norm_print(v_i2cg_INRTL, "v_i2cg_INRTL"); */
  /*   double v_mag_inrtl;  */
  /* v_mag(&v_mag_inrtl, v_rel_eci); */
  /* v_scale(adrag_i2cg_INRTL, v_rel_eci, -1/2* *density * 2.2 * 0.01 * M_PI / 1000. / 1000.* v_mag_inrtl); */
  /* v_norm_print(adrag_i2cg_INRTL, "adrag_i2cg_INRTL"); */
  //  MPI_Finalize();exit(0);
  // !!!!!!!!!!!!!!!!!! END OF ERASE LINES BELOW

  if ( strcmp(INTEGRATOR->format_density_driver, "density_file") != 0 ){ // the user does not choose to directly input the density from a file
    INTEGRATOR->density[index_in_driver_interpolated] = *density;
  }

  m_x_v(adrag_i2cg_INRTL, T_ntw_2_inrtl, a_i2cg_NTW);
  } // end of no collision with VCM as colllision input file
  else{
    double v_rel_eci_normalized[3];
    v_norm(v_rel_eci_normalized, v_rel_eci);
    double v_rel_eci_mag;
  v_mag(&v_rel_eci_mag, v_rel_eci);
  v_scale(adrag_i2cg_INRTL, v_rel_eci_normalized,  -0.5 * INTEGRATOR->bc_vcm * (*density) *v_rel_eci_mag*v_rel_eci_mag );
  //  printf("%d %d %e %e %e\n", INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb, INTEGRATOR->bc_vcm , *density, v_rel_eci_mag);
  /* v_print(adrag_i2cg_INRTL, "adrag_i2cg_INRTL"); */
  /* v_print(v_rel_eci_normalized, "v_rel_eci_normalized"); */
  

  }

  v_copy(SC->a_i2cg_INRTL_drag, adrag_i2cg_INRTL);
  if (iDebugLevel >= 2){
    printf("--- (compute_drag) Just got out of compute_drag ... (iProc %d | iii = %d, eee = %d)\n", iProc, INTEGRATOR->sc_main_nb, INTEGRATOR->sc_ensemble_nb   );
  }
					       //v_norm_print(adrag_i2cg_INRTL, "adrag_i2cg_INRTL");
  

  /*   // # TOTAL NORMALIZED DRAG COEFFICIENT */
  /*   double A_ref_tot = 0; */
  /*   double A_ref; */
  /*   double cd_tot_norm = 0; */
  /*   for (sss = 0; sss < INTEGRATOR->nb_surfaces; sss++){ */
  /*     if (cos_angle_v_sc_normal[sss] >= 0){ */
  /*       A_ref = INTEGRATOR->surface[sss].area * cos_angle_v_sc_normal[sss] * 1000000; // convert km^2 to km */

  /*         A_ref_tot = A_ref_tot + A_ref; */
  /* 	cd_tot_norm = cd_tot_norm + cd[sss]*A_ref; */
  /* 	//	printf("cd[%d] * A_ref: %e\n",sss, cd[sss]*A_ref); */
  /*       } */
  /*   } */

  /* 	etprint(et, "time"); */
  //	   printf("%e\n",cd_tot_norm);
  //  printf("A_ref_tot %f\n", A_ref_tot);
  SC->density_here = *density;
  if (INTEGRATOR->coll_vcm != 1){
  SC->INTEGRATOR.sum_cd_a_cos = cd_tot_norm / 1000000;// used in the kalman filter. convert back to use km^2 on the area
  SC->INTEGRATOR.cd_tot_norm = cd_tot_norm / A_ref_tot;// # see sutton09a equation 9
  SC->INTEGRATOR.A_ref_tot =  A_ref_tot; // A_ref_tot in m^2 but in generate_ephemerides i convert it to km^2 for the output
  cd_tot_norm = cd_tot_norm / A_ref_tot;
  } // end of no collision with VCM as colllision input file
    else{
    SC->INTEGRATOR.sum_cd_a_cos = 0;// non applicable
    SC->INTEGRATOR.cd_tot_norm = 0;// non applicable
    SC->INTEGRATOR.A_ref_tot =  0;// non applicable

  }

  if (SC->INTEGRATOR.write_given_output == 1){
    /* 	  fprintf(INTEGRATOR->file_given_output, "%e ", cd_tot_norm ); */
    /* fprintf(INTEGRATOR->file_given_output, "%e \n", total_cross_area ); */
    //	    printf( "%s %e\n", timestamp_isoc, SC->INTEGRATOR.cd_tot_norm );
  }
  /* 	if (cd_tot_norm < 2.9){ */
  /* 	  printf("XXXXXXXXXXXXXX\nXXXXXXXXXXXXX\n"); */
  /* 	} */
  //  printf("Cd total normalized: %f %e\n",cd_tot_norm, A_ref_tot);
  //	printf("cd[%d]: %f | %e %f\n", sss, cd,cd*INTEGRATOR->surface[sss].area*cos_angle_v_sc_normal, OPTIONS->surface[0].acco_coeff );






  return 0;
}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           shadow_light
//  Purpose:        Returns if the SC is in the umbra/penumbra of Earth
//  Assumptions:    None.
//  References      Vallado version 3 section 5.3.2
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/04/2015    |   ---     | Initial Implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int shadow_light( char      shadow[256],
		  double    r_i2cg_INRTL[3],
		  double    et,
		  PARAMS_T  *PARAMS) 
{

  /* Declarations */
  double x[6];
  double lt;
  double r_earth2sun_J2000[3];
  double mag_r_earth2sun_J2000;
  double r_earth2sun_J2000_normalized[3];
  double r_sun2earth_J2000_normalized[3];
  double r_i2cg_INRTL_normalized[3];
  double mag_r_i2cg_INRTL;
  double angle_minus_sun_to_sc;
  double cos_angle_minus_sun_to_sc;
  double sc_horiz, sc_vert; // satellite horizontal and vertical distances from the Sun-Earth line
  double pen_vert, umb_vert; // vertical lengths of the penumbra and the umbra
  double alpha_pen;
  double alpha_umb;
  
  /* Algorithm */
  
  /* Earth to Sun vector */
  spkez_c(10, et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. 
  r_earth2sun_J2000[0] = x[0];
  r_earth2sun_J2000[1] = x[1];
  r_earth2sun_J2000[2] = x[2];

  /* Umbra and penumbra angles */
  v_mag(&mag_r_earth2sun_J2000, r_earth2sun_J2000);
  alpha_umb = asin( ( PARAMS->SUN.radius - PARAMS->EARTH.radius ) / mag_r_earth2sun_J2000 ); // see errata of Vallado3
  alpha_pen = asin( ( PARAMS->SUN.radius + PARAMS->EARTH.radius ) / mag_r_earth2sun_J2000 ); // see errata of Vallado3

  /* Angle between the Sun-Earth direction and the Earth-satellite direction */
  v_norm(r_earth2sun_J2000_normalized, r_earth2sun_J2000);
  v_scale(r_sun2earth_J2000_normalized, r_earth2sun_J2000_normalized, -1.0);
  v_norm(r_i2cg_INRTL_normalized, r_i2cg_INRTL);
  v_dot(&cos_angle_minus_sun_to_sc, r_sun2earth_J2000_normalized, r_i2cg_INRTL_normalized);

  angle_minus_sun_to_sc = acos(cos_angle_minus_sun_to_sc);

  /* Calculate if the SC is in the umbra or penumbra */
  if (cos_angle_minus_sun_to_sc > 0){

    v_mag(&mag_r_i2cg_INRTL, r_i2cg_INRTL);
    sc_horiz = mag_r_i2cg_INRTL * cos( angle_minus_sun_to_sc );
    sc_vert =  mag_r_i2cg_INRTL * sin( angle_minus_sun_to_sc );
    pen_vert = PARAMS->EARTH.radius + tan( alpha_pen ) * sc_horiz;
  
    strcpy(shadow, "none");
    if ( sc_vert < pen_vert ){
      strcpy(shadow, "penumbra");
      umb_vert = PARAMS->EARTH.radius - tan( alpha_umb ) * sc_horiz;
      if ( sc_vert < umb_vert ){
	strcpy(shadow, "umbra");
      }
    }
    else{
      strcpy(shadow, "light");
    }
 
  }

  else{
    strcpy(shadow, "light");
  }

  return 0;
}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           shadow_light
//  Purpose:        Returns if the SC is in the umbra/penumbra of Moon
//  Assumptions:    None.
//  References      Vallado version 3 section 5.3.2
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/04/2015    |   ---     | Initial Implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int shadow_light_moon( char      shadow[256],
		  double    r_i2cg_INRTL[3],
		  double    et,
		  PARAMS_T  *PARAMS) 
{

  /* Declarations */
  double x[6];
  double lt;
  double r_moon2sun_J2000[3];
  double mag_r_moon2sun_J2000;
  double r_moon2sun_J2000_normalized[3];
  double r_sun2moon_J2000_normalized[3];
  double mag_r_moon2sc_J2000;
  double angle_minus_sun_to_sc;
  double cos_angle_minus_sun_to_sc;
  double sc_horiz, sc_vert; // satellite horizontal and vertical distances from the Sun-Moon line
  double pen_vert, umb_vert; // vertical lengths of the penumbra and the umbra
  double alpha_pen;
  double alpha_umb;
  
  /* Algorithm */
  
  /* Moon to Sun vector */
  spkez_c(10, et, "J2000", "NONE", 301, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. 

  r_moon2sun_J2000[0] = x[0];
  r_moon2sun_J2000[1] = x[1];
  r_moon2sun_J2000[2] = x[2];

  /* Umbra and penumbra angles */
  v_mag(&mag_r_moon2sun_J2000, r_moon2sun_J2000);
  alpha_umb = asin( ( PARAMS->SUN.radius - PARAMS->MOON.radius ) / mag_r_moon2sun_J2000 ); // see errata of Vallado3
  alpha_pen = asin( ( PARAMS->SUN.radius + PARAMS->MOON.radius ) / mag_r_moon2sun_J2000 ); // see errata of Vallado3

  /* Angle between the Sun-Moon direction and the Moon-satellite direction */
  v_norm(r_moon2sun_J2000_normalized, r_moon2sun_J2000);
  v_scale(r_sun2moon_J2000_normalized, r_moon2sun_J2000_normalized, -1.0);

  double x_earth[6];
  double lt_earth;
spkez_c(10, et, "J2000", "NONE", 399, x_earth, &lt_earth); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. 

 double r_earth2sun_J2000[3];
  r_earth2sun_J2000[0] = x_earth[0];
  r_earth2sun_J2000[1] = x_earth[1];
  r_earth2sun_J2000[2] = x_earth[2];
  
  double r_moon2earth_J2000[3];
  v_sub(r_moon2earth_J2000, r_moon2sun_J2000, r_earth2sun_J2000);

  double r_moon2sc_J2000[3];
  v_add(r_moon2sc_J2000, r_moon2earth_J2000, r_i2cg_INRTL);
  
  double r_moon2sc_J2000_normalized[3];
  v_norm(r_moon2sc_J2000_normalized, r_moon2sc_J2000);
  v_dot(&cos_angle_minus_sun_to_sc, r_sun2moon_J2000_normalized, r_moon2sc_J2000_normalized);

  angle_minus_sun_to_sc = acos(cos_angle_minus_sun_to_sc);

  /* Calculate if the SC is in the umbra or penumbra */
  if (cos_angle_minus_sun_to_sc > 0){

    v_mag(&mag_r_moon2sc_J2000, r_moon2sc_J2000);
    sc_horiz = mag_r_moon2sc_J2000 * cos( angle_minus_sun_to_sc );
    sc_vert =  mag_r_moon2sc_J2000 * sin( angle_minus_sun_to_sc );
    pen_vert = PARAMS->MOON.radius + tan( alpha_pen ) * sc_horiz;
  
    strcpy(shadow, "none");
    if ( sc_vert < pen_vert ){
      strcpy(shadow, "penumbra");
      umb_vert = PARAMS->MOON.radius - tan( alpha_umb ) * sc_horiz;
      if ( sc_vert < umb_vert ){
	strcpy(shadow, "umbra");
      }
    }
    else{
      strcpy(shadow, "light");
    }
 
  }

  else{
    strcpy(shadow, "light");
  }

  return 0;
}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           load_params
//  Purpose:        Loads params. Ellipsoid based on WGS84
//  Assumptions:    None.
//  References      Various
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/20/2015    |   ---     | Initial Implementation
//      | C. Bussy-Virat| 10/14/2015    |   ---     | Change the inputs
//
/////////////////////////////////////////////////////////////////////////////////////////
//newstructure
//int load_params( PARAMS_T *PARAMS, char main_directory_location[256], int iDebugLevel, char earth_fixed_frame[100],   double use_ap_hist, int iProc) {
int load_params( PARAMS_T *PARAMS,  int iDebugLevel, char earth_fixed_frame[100],   double use_ap_hist, int iProc, char path_to_spice[256], int degree, int gravity_map_use, int earth_pressure) {
//newstructure


      /* Set up the geophysical quantities.  At last check these were the values used by Space Command and SGP4 */
      PARAMS->geophs[ 0 ] =    1.082616e-3;   // J2
      PARAMS->geophs[ 1 ] =   -2.53881e-6;    // J3
      PARAMS->geophs[ 2 ] =   -1.65597e-6;    // J4
      PARAMS->geophs[ 3 ] =    7.43669161e-2; // KE
      PARAMS->geophs[ 4 ] =    120.0;         // QO
      PARAMS->geophs[ 5 ] =    78.0;          // SO
      PARAMS->geophs[ 6 ] =    6378.135;      // ER
      PARAMS->geophs[ 7 ] =    1.0;           // AE
  
  int ii;

  PARAMS->EARTH.flattening    = 1/298.257223560;
  PARAMS->EARTH.radius        = 6378.137;
  PARAMS->EARTH.GRAVITY.mu    = 398600.4418; // km^3/s^2
  PARAMS->EARTH.GRAVITY.j2    = 1.081874e-3;
  strcpy(PARAMS->EARTH.earth_fixed_frame, earth_fixed_frame);

  PARAMS->EARTH.GRAVITY.dlat_map = 0.5;
  PARAMS->EARTH.GRAVITY.dlon_map = 1.;
  PARAMS->EARTH.GRAVITY.dradius_map = 10.;//;22. ;//PARAMS->EARTH.GRAVITY.dlon_map * 100.;
  PARAMS->EARTH.GRAVITY. min_lat_map = -90;//-0.11;
  PARAMS->EARTH.GRAVITY.max_lat_map = 90;//0.11;
  PARAMS->EARTH.GRAVITY.min_radius_map = PARAMS->EARTH.radius + 400;//200.;
  PARAMS->EARTH.GRAVITY.max_radius_map = PARAMS->EARTH.radius + 600;//501;//40000.;
  PARAMS->EARTH.GRAVITY.earth_min_radius_map = PARAMS->EARTH.radius + 450;//200.;
  PARAMS->EARTH.GRAVITY.earth_max_radius_map = PARAMS->EARTH.radius + 600;//501;//40000.;
  
  PARAMS->EARTH.GRAVITY.nlon_map = (int)(ceil( 360./PARAMS->EARTH.GRAVITY.dlon_map)) + 1;// nb lon bins
  PARAMS->EARTH.GRAVITY.nlat_map = (int)(ceil( (PARAMS->EARTH.GRAVITY.max_lat_map- PARAMS->EARTH.GRAVITY.min_lat_map)/PARAMS->EARTH.GRAVITY.dlat_map)) + 1;// nb lat bins
  PARAMS->EARTH.GRAVITY.nradius_map = (int)(ceil( (PARAMS->EARTH.GRAVITY.max_radius_map - PARAMS->EARTH.GRAVITY.min_radius_map)/PARAMS->EARTH.GRAVITY.dradius_map)) + 1;


  // for the earth pressure map (the radius maping is the same as for the gravity map)
  // // zenith
  PARAMS->EARTH.GRAVITY.dzenith_map = 3.;//;22. ;//PARAMS->EARTH.GRAVITY.dlon_map * 100.;
  PARAMS->EARTH.GRAVITY.min_zenith_map = 0;
  PARAMS->EARTH.GRAVITY.max_zenith_map = 180.;

  // // azim and elevation of the Earth element (these azim and elev are not technically azimuth or elevation angles but these nominations are used)
  PARAMS->EARTH.GRAVITY.dazim_elt_map = 30.;
  PARAMS->EARTH.GRAVITY.min_azim_elt_map = 0;
  PARAMS->EARTH.GRAVITY.max_azim_elt_map = 360.;
  
  PARAMS->EARTH.GRAVITY.delev_elt_map = 1.;
  PARAMS->EARTH.GRAVITY.min_elev_elt_map = 0;
  PARAMS->EARTH.GRAVITY.max_elev_elt_map = acos(PARAMS->EARTH.radius/PARAMS->EARTH.GRAVITY.max_radius_map) * 180./M_PI;

  // // azim and elevation of the surface normal vector (these azim and elev are not technically azimuth or elevation angles but these nominations are used)
  PARAMS->EARTH.GRAVITY.dazim_surf_map = 30.;
  PARAMS->EARTH.GRAVITY.min_azim_surf_map = 0;
  PARAMS->EARTH.GRAVITY.max_azim_surf_map = 360.;

  PARAMS->EARTH.GRAVITY.delev_surf_map = 3.;
  PARAMS->EARTH.GRAVITY.min_elev_surf_map = PARAMS->EARTH.GRAVITY.max_elev_elt_map; // any surface which elev_surf_map lower than this value won't see any Earth element
  PARAMS->EARTH.GRAVITY.max_elev_surf_map = 180;

    PARAMS->EARTH.GRAVITY.nzenith_map = (int)(ceil( (PARAMS->EARTH.GRAVITY.max_zenith_map - PARAMS->EARTH.GRAVITY.min_zenith_map)/PARAMS->EARTH.GRAVITY.dzenith_map)) + 1;
    PARAMS->EARTH.GRAVITY.nazim_elt_map = (int)(ceil( (PARAMS->EARTH.GRAVITY.max_azim_elt_map - PARAMS->EARTH.GRAVITY.min_azim_elt_map)/PARAMS->EARTH.GRAVITY.dazim_elt_map)); //not +1 here because we wawnt to exlude 360 otherwise the pressure from the elemetn at 0 deg and the element at 360 deg are both counted, while they correspond to a single and same element
        PARAMS->EARTH.GRAVITY.nelev_elt_map = (int)(ceil( (PARAMS->EARTH.GRAVITY.max_elev_elt_map - PARAMS->EARTH.GRAVITY.min_elev_elt_map)/PARAMS->EARTH.GRAVITY.delev_elt_map)) + 1;
        PARAMS->EARTH.GRAVITY.nazim_surf_map = (int)(ceil( (PARAMS->EARTH.GRAVITY.max_azim_surf_map - PARAMS->EARTH.GRAVITY.min_azim_surf_map)/PARAMS->EARTH.GRAVITY.dazim_surf_map)) + 1;
        PARAMS->EARTH.GRAVITY.nelev_surf_map = (int)(ceil( (PARAMS->EARTH.GRAVITY.max_elev_surf_map - PARAMS->EARTH.GRAVITY.min_elev_surf_map)/PARAMS->EARTH.GRAVITY.delev_surf_map)) + 1;

  // end of for the earth pressure map 



  
  PARAMS->MOON.flattening     = 0.0012;//https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
  PARAMS->MOON.radius         = 1737.0;//https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
  PARAMS->MOON.GRAVITY.mu     = 4.902801076000e+003; // (value from STK) 
  PARAMS->MOON.GRAVITY.j2     = 202.7e-6;//https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
  strcpy(PARAMS->MOON.earth_fixed_frame, "");

  PARAMS->SUN.flattening      = 0.0;
  PARAMS->SUN.radius          = 696300.0;
  PARAMS->SUN.GRAVITY.mu      = 1.327122000000E+011; // in km^3 / s^2 (value from STK)
  PARAMS->SUN.GRAVITY.j2      = 0.0;
  strcpy(PARAMS->SUN.earth_fixed_frame, "");
    
  if (iDebugLevel >= 2){
        if (iProc == 0) printf("--- (load_params) Loading gravity...\n");
  }
  //newstructure
  load_gravity(  &PARAMS->EARTH.GRAVITY,  path_to_spice );
  //  load_gravity(  &PARAMS->EARTH.GRAVITY, main_directory_location );
  //newstructure
  PARAMS->EARTH.GRAVITY.radius = 6378.137;
   

  // Flags for NRLMSIS-00e
  PARAMS->ATMOSPHERE.flags.switches[0]=0;
  for (ii=1;ii<24;ii++){
        PARAMS->ATMOSPHERE.flags.switches[ii]=1;
  }
    // this allow to use msis with historical ap rather than daily ap
  if (use_ap_hist == 1){
      PARAMS->ATMOSPHERE.flags.switches[9]=-1;
  }
  // !!!!!!!!!!!!!!!!!! TO ERASE AND UNCOMMENT BLOCK ABOVE
  /* FILE *fp_switch = NULL; */
  /* char *line = NULL; */
  /* size_t len = 0; */
  /* fp_switch = fopen("./input/switch_for_msis.txt","r"); */
  /* for (ii=1;ii<24;ii++){ */
  /*   getline(&line,&len,fp_switch); */
  /*   sscanf(line, "%d", &PARAMS->ATMOSPHERE.flags.switches[ii]); */
  /*      printf("[%d]: %d\n",ii, PARAMS->ATMOSPHERE.flags.switches[ii]); */
  /* } */
  /* fclose(fp_switch); */
  // !!!!!!!!!!!END OF TO ERASE AND UNCOMMENT BLOCK ABOVE





  if (earth_pressure == 1){
  char text[200];
  strcpy(text, "");
  strcpy(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, "earthPres");
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, "_alt");
  sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.earth_min_radius_map-PARAMS->EARTH.GRAVITY.radius);
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, text);
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, "-");
  sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.earth_max_radius_map-PARAMS->EARTH.GRAVITY.radius);
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, text);
    strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, "-");
    sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.dradius_map );
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, text);
    strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, "_zenith");
    sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.dzenith_map );
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, text);
    strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, "_elev_surf");
    sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.delev_surf_map );
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, text);
    strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, "_azim_surf");
    sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.dazim_surf_map );
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, text);
    strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, "_elev_elt");
    sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.delev_elt_map );
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, text);
    strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, "_azim_elt");
    sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.dazim_elt_map );
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, text);
  
  //strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, ".txt");  // txt
  strcat(PARAMS->EARTH.GRAVITY.filename_earth_pressure_map, ".bin");
  printf("Map: <%s>\n", PARAMS->EARTH.GRAVITY.filename_earth_pressure_map);

  //      build_earth_pressure_map(&(PARAMS->EARTH.GRAVITY), iProc);
	read_earth_pressure_map(&(PARAMS->EARTH.GRAVITY), iProc);

  }

  
  if (gravity_map_use == 1){

    
  char text[200];
  strcpy(text, "");
  strcpy(PARAMS->EARTH.GRAVITY.filename_gravity_map, "grav");
  sprintf(text, "%d", degree);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, text);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, "_alt");
  sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.min_radius_map-PARAMS->EARTH.GRAVITY.radius);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, text);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, "-");
  sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.max_radius_map-PARAMS->EARTH.GRAVITY.radius);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, text);
    strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, "-");
  sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.dradius_map );
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, text);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, "_lat");
  sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.min_lat_map);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, text);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, "-");
  sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.max_lat_map);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, text);
    strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, "-");
  sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.dlat_map);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, text);
    strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, "_lon");
  sprintf(text, "%.1f", PARAMS->EARTH.GRAVITY.dlon_map);
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, text);
  //strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, ".txt");  // txt
  strcat(PARAMS->EARTH.GRAVITY.filename_gravity_map, ".bin");
  //printf("<%s>\n", PARAMS->EARTH.GRAVITY.filename_gravity_map);


  


  //  build_gravity_map( &(PARAMS->EARTH.GRAVITY), degree,  iProc);
  	  read_gravity_map( &(PARAMS->EARTH.GRAVITY), degree,  iProc);
  }


  //  printf("AAAAAAA = %e\n", PARAMS->EARTH.GRAVITY.gravity_map[iradius][ilat][ilon][0]);
  
  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           test_print
//  Purpose:        Very useless function thatpr ints something to see if the code runs fine up to this line
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 06/08/2015    |   ---     | Initial implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int test_print( char to_print[1]  )

{

  char to_print_more[10];
  int i;
  strcpy(to_print_more, "");
  //  for (i = 0; i < 10; i++ ){
      
    strcat(to_print_more, to_print);
    //  }
  printf("**********\n");
  for (i = 0; i < 3; i++){
    printf("%s\n", to_print_more);
  }
  printf("**********\n");
  return 0;
}



int test_print_iproc( int iProc, char to_print[1] )

{
  char to_print_more[20];
  int i;
  strcpy(to_print_more, "");
  for (i = 0; i < 10; i++ ){
    strcat(to_print_more, to_print);
  }    
  char iproc_str[5];
  sprintf(iproc_str, "%d", iProc);
  strcat(to_print_more, " ");
  strcat(to_print_more, iproc_str);
  printf("**********\n");
  for (i = 0; i < 3; i++){
    printf("%s\n", to_print_more);
  }
  printf("**********\n");
  return 0;
}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           print_test
//  Purpose:        Prints something to see if the code runs fine up to this line
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 12/09/2015    |   ---     | Initial implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int print_test(  )

{
  printf("\n****************** TEST ******************\n****************** TEST ******************\n\n");
  return 0;
}



int print_teste(  )

{
  printf("\n****************** TEST ******************\n****************** TEST ******************\n\n");
  MPI_Finalize();exit(0);
  return 0;
}

int set_attitude( ATTITUDE_T attitude, int index_in_attitude_interpolated, OPTIONS_T *OPTIONS, int file_is_quaternion ){
  // this function sets the attitude in case no ensemble at all is run on the attitude


  /* Declarations */
  //  printf("%d\n", index_in_attitude_interpolated);
	  //	 printf("%d %f %f\n", index_in_attitude_interpolated, OPTIONS->pitch[index_in_attitude_interpolated], attitude.pitch[index_in_attitude_interpolated]);
  if (file_is_quaternion == 0){
  attitude.pitch[index_in_attitude_interpolated] = OPTIONS->pitch[index_in_attitude_interpolated];
  attitude.roll[index_in_attitude_interpolated] = OPTIONS->roll[index_in_attitude_interpolated];
  attitude.yaw[index_in_attitude_interpolated] = OPTIONS->yaw[index_in_attitude_interpolated];
  attitude.order_pitch[index_in_attitude_interpolated] = OPTIONS->order_pitch[index_in_attitude_interpolated];
  attitude.order_roll[index_in_attitude_interpolated] = OPTIONS->order_roll[index_in_attitude_interpolated];
  attitude.order_yaw[index_in_attitude_interpolated] = OPTIONS->order_yaw[index_in_attitude_interpolated];
  }
  else{
    attitude.quaternion[index_in_attitude_interpolated][0] = OPTIONS->quaternion[index_in_attitude_interpolated][0];
    attitude.quaternion[index_in_attitude_interpolated][1] = OPTIONS->quaternion[index_in_attitude_interpolated][1];
    attitude.quaternion[index_in_attitude_interpolated][2] = OPTIONS->quaternion[index_in_attitude_interpolated][2];
    attitude.quaternion[index_in_attitude_interpolated][3] = OPTIONS->quaternion[index_in_attitude_interpolated][3];
  }

  /* Algorithm */

  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           gitm_density
//  Purpose:        Calculate the density at the position of the satellite using the GITM thermospheric model 
//  Assumptions:    The epoch time has to be newer than 2000. Other assumptions are commented in the code
//  References      Various
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 01/24/2016    |   ---     | Initial implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int gitm_density( double *density, double et, double altitude, double latitude, double longitude, PARAMS_T *PARAMS) {

  // Declarations
  int pos_for_lon_gitm_right_before, iDataLength_gitm_right_before, iHeaderLength_gitm_right_before;
  int time_gitm_right_before_temp[7];
  int j,k;
  FILE *gitm_file_gitm_right_before;

  int i;
  int pos_for_lon_gitm_right_after, iDataLength_gitm_right_after, iHeaderLength_gitm_right_after;
  int time_gitm_right_after_temp[7];


  FILE *gitm_file_gitm_right_after;

  char name_gitm_file_right_before_epoch_of_sc[256], name_gitm_file_right_after_epoch_of_sc[256];
  char time_sc[256];
  double et_gitm_file;
  double delta_et, delta_et_right_before, delta_et_right_after;
  int found_2_closest_GITM_file_date;
  int need_to_read_GITM_before, need_to_read_GITM_after;
 
  et2utc_c(et, "ISOC" ,0 ,255 , time_sc);

  // The goal of this block is to find the 2 GITM files that have dates around the epoch of the sc
  if (PARAMS->ATMOSPHERE.is_first_step_of_run == 1){ //if it is the first step of the run then find the two cloest GITM file dates to the epoch start tie of the sc
    i = 0;
    found_2_closest_GITM_file_date = 0;
    while ( (i < PARAMS->ATMOSPHERE.nb_gitm_file) && (found_2_closest_GITM_file_date == 0) ){
      et_gitm_file = PARAMS->ATMOSPHERE.array_gitm_date[i][0];
      // Find the 2 GITM files that have dates around the epoch of the sc
      delta_et = et_gitm_file - et;
      if (delta_et <= 0){ // the date of the GITM file is older than the epoch of the sc
	found_2_closest_GITM_file_date = 0;
      }
      else {//if (delta_et > 0){
	if (i > 0){
	  PARAMS->ATMOSPHERE.gitm_index_lower = i-1;
	  PARAMS->ATMOSPHERE.gitm_index_higher = i;
	  found_2_closest_GITM_file_date = 1;
	}
	else{
	  printf("The oldest GITM file has to be older than epoch start time of the satellite. The program will stop.\n");
	  exit(0);
	}
      }
      i = i+1;
    }

    need_to_read_GITM_before = 1;
    need_to_read_GITM_after = 1;
  } // end of if first step of run
  else{ //if it is not the first step of run then just look if the epoch time of the sc got newer than the subsequent GITM file date
    if(et >= PARAMS->ATMOSPHERE.array_gitm_date[PARAMS->ATMOSPHERE.gitm_index_higher][0]){ // if the epoch got newer than the time of the GITM subsequent file then the new subsequent GITM file is the next file in the array PARAMS->ATMOSPHERE.array_gitm_file since this latter is arranged in ascending order of times. Then, read the new file
      PARAMS->ATMOSPHERE.gitm_index_lower = PARAMS->ATMOSPHERE.gitm_index_higher;
      PARAMS->ATMOSPHERE.gitm_index_higher = PARAMS->ATMOSPHERE.gitm_index_higher + 1; // the GITM files are arranged in the ascending order

      need_to_read_GITM_before = 0;
      need_to_read_GITM_after = 1;
    }
    else{ // otherwise no need to do anything: we already read the 2 GITM files
      need_to_read_GITM_before = 0;
      need_to_read_GITM_after = 0;
    }
  } // end of if it is not the first step of run 
  // end of the goal of this block is to find the 2 GITM files that have dates around the epoch of the sc


  strcpy(name_gitm_file_right_before_epoch_of_sc, PARAMS->ATMOSPHERE.array_gitm_file[PARAMS->ATMOSPHERE.gitm_index_lower]);
  strcpy(name_gitm_file_right_after_epoch_of_sc, PARAMS->ATMOSPHERE.array_gitm_file[PARAMS->ATMOSPHERE.gitm_index_higher]);
  delta_et_right_before = PARAMS->ATMOSPHERE.array_gitm_date[PARAMS->ATMOSPHERE.gitm_index_lower][0] - et; 
  delta_et_right_after = PARAMS->ATMOSPHERE.array_gitm_date[PARAMS->ATMOSPHERE.gitm_index_higher][0] - et;
  
  if (need_to_read_GITM_before == 1){ // read the GITM file that has a date right before the epoch date of the sc. This happens only at the first step of the run

    /*   /\******************************************************************\/ */
    /*   /\*************************** READ GITM RIGHT BEFORE ***************\/ */
    /*   /\******************************************************************\/ */

  
    gitm_file_gitm_right_before = fopen(name_gitm_file_right_before_epoch_of_sc,"rb");


    // TIME
    fseek(gitm_file_gitm_right_before,48+ 4+PARAMS->ATMOSPHERE.nVars_gitm*(4+4+40), SEEK_CUR	 );
    fread(&time_gitm_right_before_temp,sizeof(time_gitm_right_before_temp),1,gitm_file_gitm_right_before);
    /* printf("\n\n"); */
    /*  for (i = 0; i < 7; i++){ */
    /*     printf("time: %d\n", time_gitm_right_before_temp[i]); */
    /*    } */

    // SKIP VARIABLES BEFORE THE DENSITY
    iHeaderLength_gitm_right_before = 8L + 4+4 +	3*4 + 4+4 + 4 + 4+4 + PARAMS->ATMOSPHERE.nVars_gitm*40 + PARAMS->ATMOSPHERE.nVars_gitm*(4+4) +  7*4 + 4+4;
    iDataLength_gitm_right_before = PARAMS->ATMOSPHERE.nLons_gitm*PARAMS->ATMOSPHERE.nLats_gitm*PARAMS->ATMOSPHERE.nAlts_gitm*8 + 4+4;
    pos_for_lon_gitm_right_before = iHeaderLength_gitm_right_before + 3*iDataLength_gitm_right_before ;// density is the 4th variable in the GITM file
    fseek(gitm_file_gitm_right_before, pos_for_lon_gitm_right_before, SEEK_SET);


    // READ THE DENSITY AT THE GIVEN LON/LAT/ALT
    fseek(gitm_file_gitm_right_before, 4, SEEK_CUR );
    for (k = 0; k < PARAMS->ATMOSPHERE.nAlts_gitm; k++){
      for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
	for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
	  fread(&PARAMS->ATMOSPHERE.density_gitm_right_before[i][j][k], sizeof(PARAMS->ATMOSPHERE.density_gitm_right_before[i][j][k]), 1, gitm_file_gitm_right_before);
	}
      }
    }

    /******************************************************************/
    /******************** END OF READ GITM RIGHT BEFORE ***************/
    /******************************************************************/

    fclose(gitm_file_gitm_right_before);
  }  // end of read the GITM file that has a date right before the epoch date of the sc. This happens only at the first step of the run 



  if (need_to_read_GITM_after == 1){ // this happens at the first step of the run or when a new subsequent GITM file needs to be found (see previous comments above)

    if (PARAMS->ATMOSPHERE.is_first_step_of_run == 0){

      for (i = 0; i < 7; i++){
	time_gitm_right_before_temp[i] = time_gitm_right_after_temp[i];
      }

      for (k = 0; k < PARAMS->ATMOSPHERE.nAlts_gitm; k++){
	for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
	  for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
	    PARAMS->ATMOSPHERE.density_gitm_right_before[i][j][k] = PARAMS->ATMOSPHERE.density_gitm_right_after[i][j][k];
	  }
	}
      }

    }


    /******************************************************************/
    /**************** NOW READ THE NEW GITM RIGHT AFTER ***************/
    /******************************************************************/
    gitm_file_gitm_right_after = fopen(name_gitm_file_right_after_epoch_of_sc,"rb");

    // TIME
    fseek(gitm_file_gitm_right_after,48+ 4+PARAMS->ATMOSPHERE.nVars_gitm*(4+4+40), SEEK_CUR	 );
    fread(&time_gitm_right_after_temp,sizeof(time_gitm_right_after_temp),1,gitm_file_gitm_right_after);
    /* printf("\n\n"); */
    /*  for (i = 0; i < 7; i++){ */
    /*     printf("time: %d\n", time_gitm_right_after_temp[i]); */
    /*    } */

    // SKIP VARIABLES AFTER THE DENSITY
    iHeaderLength_gitm_right_after = 8L + 4+4 +	3*4 + 4+4 + 4 + 4+4 + PARAMS->ATMOSPHERE.nVars_gitm*40 + PARAMS->ATMOSPHERE.nVars_gitm*(4+4) +  7*4 + 4+4;
    iDataLength_gitm_right_after = PARAMS->ATMOSPHERE.nLons_gitm*PARAMS->ATMOSPHERE.nLats_gitm*PARAMS->ATMOSPHERE.nAlts_gitm*8 + 4+4;
    pos_for_lon_gitm_right_after = iHeaderLength_gitm_right_after + 3*iDataLength_gitm_right_after ;// density is the 4th variable in the GITM file
    fseek(gitm_file_gitm_right_after, pos_for_lon_gitm_right_after, SEEK_SET);


    // READ THE DENSITY AT THE GIVEN LON/LAT/ALT
    fseek(gitm_file_gitm_right_after, 4, SEEK_CUR );
    for (k = 0; k < PARAMS->ATMOSPHERE.nAlts_gitm; k++){
      for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
	for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
	  fread(&PARAMS->ATMOSPHERE.density_gitm_right_after[i][j][k], sizeof(PARAMS->ATMOSPHERE.density_gitm_right_after[i][j][k]), 1, gitm_file_gitm_right_after);

	}
      }
    }

    /******************************************************************/
    /******************** END OF READ GITM RIGHT AFTER ***************/
    /******************************************************************/
    fclose(gitm_file_gitm_right_after);
  }

  /* printf("\n%s | %s | %s\n", name_gitm_file_right_before_epoch_of_sc, time_sc, name_gitm_file_right_after_epoch_of_sc);   */
  /* printf("longitude = %f \n", PARAMS->ATMOSPHERE.longitude_gitm[11][40][21]); */
  /* printf("latitude = %f \n", PARAMS->ATMOSPHERE.latitude_gitm[11][40][21]); */
  /* printf("altitude = %f \n", PARAMS->ATMOSPHERE.altitude_gitm[11][40][21]); */
  /* printf("density = %e | %e\n", PARAMS->ATMOSPHERE.density_gitm_right_before[11][40][21], PARAMS->ATMOSPHERE.density_gitm_right_after[11][40][21]); */



  /*   /\******************************************************************\/ */
  /*   /\********** INTERPOLATION LON/LAT/ALT WITH GITM RIGHT BEFORE ******\/ */
  /*   /\******************************************************************\/ */

  /*      // find the 2 closest longitudes in GITM file to the sc longitude */
  /*      double delta_longitude_GITM_file_right_before, delta_longitude_right_lower_sc_lon_in_GITM_file_right_before = -10*M_PI, delta_longitude_right_higher_sc_lon_in_GITM_file_right_before = 10*M_PI; */
  /*      double longitude_right_lower_sc_lon_in_GITM_file_right_before = 0, longitude_right_higher_sc_lon_in_GITM_file_right_before = 0; */
  /*      int index_longitude_right_lower_sc_lon_in_GITM_file_right_before = -1, index_longitude_right_higher_sc_lon_in_GITM_file_right_before = -1; */
  /*      int index_latitude_right_lower_sc_lat_in_GITM_file_right_before = -1, index_latitude_right_higher_sc_lat_in_GITM_file_right_before = -1; */
  /*      int index_altitude_right_lower_sc_alt_in_GITM_file_right_before = -1, index_altitude_right_higher_sc_alt_in_GITM_file_right_before = -1; */
  /*      for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){ */
  /*        delta_longitude_GITM_file_right_before = PARAMS->ATMOSPHERE.longitude_gitm[i][0][0] - longitude; */
  /*        if (delta_longitude_GITM_file_right_before < 0){ // GITM longitude is smaller than sc longitude */
  /* 	 if (delta_longitude_GITM_file_right_before > delta_longitude_right_lower_sc_lon_in_GITM_file_right_before){ */
  /* 	   delta_longitude_right_lower_sc_lon_in_GITM_file_right_before = delta_longitude_GITM_file_right_before; */
  /* 	   longitude_right_lower_sc_lon_in_GITM_file_right_before = PARAMS->ATMOSPHERE.longitude_gitm[i][0][0]; */
  /* 	   index_longitude_right_lower_sc_lon_in_GITM_file_right_before = i; */
  /* 	 } */
  /*        } */
  /*        if (delta_longitude_GITM_file_right_before > 0){ // GITM longitude is greater than sc longitude */
  /* 	 if (delta_longitude_GITM_file_right_before < delta_longitude_right_higher_sc_lon_in_GITM_file_right_before){ */
  /* 	   delta_longitude_right_higher_sc_lon_in_GITM_file_right_before = delta_longitude_GITM_file_right_before; */
  /* 	   longitude_right_higher_sc_lon_in_GITM_file_right_before = PARAMS->ATMOSPHERE.longitude_gitm[i][0][0]; */
  /* 	   index_longitude_right_higher_sc_lon_in_GITM_file_right_before = i; */
  /* 	 } */
  /*        } */
  /*      } */


  /* // find the 2 closest latitudes in GITM file to the sc latitude */
  /*      double delta_latitude_GITM_file_right_before, delta_latitude_right_lower_sc_lat_in_GITM_file_right_before = -10*M_PI, delta_latitude_right_higher_sc_lat_in_GITM_file_right_before = 10*M_PI; */
  /*      double latitude_right_lower_sc_lat_in_GITM_file_right_before = 0, latitude_right_higher_sc_lat_in_GITM_file_right_before = 0; */
  /*      for (i = 0; i < PARAMS->ATMOSPHERE.nLats_gitm; i++){ */
  /*        delta_latitude_GITM_file_right_before = PARAMS->ATMOSPHERE.latitude_gitm[0][i][0] - latitude; */
  /*        if (delta_latitude_GITM_file_right_before < 0){ // GITM latitude is smaller than sc latitude */
  /* 	 if (delta_latitude_GITM_file_right_before > delta_latitude_right_lower_sc_lat_in_GITM_file_right_before){ */
  /* 	   delta_latitude_right_lower_sc_lat_in_GITM_file_right_before = delta_latitude_GITM_file_right_before; */
  /* 	   latitude_right_lower_sc_lat_in_GITM_file_right_before = PARAMS->ATMOSPHERE.latitude_gitm[0][i][0]; */
  /* 	   index_latitude_right_lower_sc_lat_in_GITM_file_right_before = i; */
  /* 	 } */
  /*        } */
  /*        if (delta_latitude_GITM_file_right_before > 0){ // GITM latitude is greater than sc latitude */
  /* 	 if (delta_latitude_GITM_file_right_before < delta_latitude_right_higher_sc_lat_in_GITM_file_right_before){ */
  /* 	   delta_latitude_right_higher_sc_lat_in_GITM_file_right_before = delta_latitude_GITM_file_right_before; */
  /* 	   latitude_right_higher_sc_lat_in_GITM_file_right_before = PARAMS->ATMOSPHERE.latitude_gitm[0][i][0]; */
  /* 	   index_latitude_right_higher_sc_lat_in_GITM_file_right_before = i; */
  /* 	 } */
  /*        } */
  /*      } */
  /*      /\* print_test(); *\/ */
  /*      /\* printf("%f | %f | %f\n", latitude_right_lower_sc_lat_in_GITM_file_right_before, latitude, latitude_right_higher_sc_lat_in_GITM_file_right_before); *\/ */
       


  // find the 2 closest altitudes in GITM file to the sc altitude
  double delta_altitude_GITM;
  int index_altitude_lower;
  int have_found_index_altitude_lower;
  i = PARAMS->ATMOSPHERE.index_altitude_right_below_perigee;  // index in altitude_gitm of the altitude that is right below the perigee altitude (calculated at the initialization). Note: if the duration of the run is long (so that the sc looses a good amount of altitude, so like 6 months) then this index is set to 0 so that we go over all the altitudes and not only the ones that start below the perigee calcualted at the initialization.
  have_found_index_altitude_lower = 0;
  while ( (i < PARAMS->ATMOSPHERE.nAlts_gitm) && (have_found_index_altitude_lower == 0) ) {

    delta_altitude_GITM = PARAMS->ATMOSPHERE.altitude_gitm[0][0][i] - altitude;
    if (delta_altitude_GITM >= 0){ // GITM altitude is smaller than sc altitude	 
      index_altitude_lower = i-1;
      have_found_index_altitude_lower = 1;
    }
    i = i + 1;
  }
  if (have_found_index_altitude_lower == 0){
    printf("The altitude of the satellite is not in the range of altitudes covered by the GITM files. The programs will stop.\n");
    exit(0);
  }


  int index_longitude_lower;
  index_longitude_lower = rint( longitude / (PARAMS->ATMOSPHERE.longitude_gitm[1][0][0] - PARAMS->ATMOSPHERE.longitude_gitm[0][0][0]) ) + 1;
  int index_latitude_lower;
  index_latitude_lower = rint( ( latitude + M_PI/2) / (PARAMS->ATMOSPHERE.latitude_gitm[0][1][0] - PARAMS->ATMOSPHERE.latitude_gitm[0][0][0]) ) + 1;
  double dist_lon, dist_lat, dist_alt;
  dist_lon = ( longitude - PARAMS->ATMOSPHERE.longitude_gitm[index_longitude_lower][0][0] ) / ( PARAMS->ATMOSPHERE.longitude_gitm[index_longitude_lower+1][0][0] - PARAMS->ATMOSPHERE.longitude_gitm[index_longitude_lower][0][0] );
  dist_lat = ( latitude - PARAMS->ATMOSPHERE.latitude_gitm[0][index_latitude_lower][0] ) / ( PARAMS->ATMOSPHERE.latitude_gitm[0][index_latitude_lower+1][0] - PARAMS->ATMOSPHERE.latitude_gitm[0][index_latitude_lower][0] );
  dist_alt = ( altitude - PARAMS->ATMOSPHERE.altitude_gitm[0][0][index_altitude_lower] ) / ( PARAMS->ATMOSPHERE.altitude_gitm[0][0][index_altitude_lower+1] - PARAMS->ATMOSPHERE.altitude_gitm[0][0][index_altitude_lower] );

     
  double density_lat_lon_lower_plane_gitm_right_before;
  density_lat_lon_lower_plane_gitm_right_before = (1-dist_lon) * (1-dist_lat) * PARAMS->ATMOSPHERE.density_gitm_right_before[index_longitude_lower][index_latitude_lower][index_altitude_lower] +
    dist_lon * (1-dist_lat) * PARAMS->ATMOSPHERE.density_gitm_right_before[index_longitude_lower+1][index_latitude_lower][index_altitude_lower] +
    dist_lon * dist_lat * PARAMS->ATMOSPHERE.density_gitm_right_before[index_longitude_lower+1][index_latitude_lower+1][index_altitude_lower] +
    (1-dist_lon) * dist_lat * PARAMS->ATMOSPHERE.density_gitm_right_before[index_longitude_lower][index_latitude_lower+1][index_altitude_lower];

  double density_lat_lon_higher_plane_gitm_right_before;
  density_lat_lon_higher_plane_gitm_right_before = (1-dist_lon) * (1-dist_lat) * PARAMS->ATMOSPHERE.density_gitm_right_before[index_longitude_lower][index_latitude_lower][index_altitude_lower+1] +
    dist_lon * (1-dist_lat) * PARAMS->ATMOSPHERE.density_gitm_right_before[index_longitude_lower+1][index_latitude_lower][index_altitude_lower+1] +
    dist_lon * dist_lat * PARAMS->ATMOSPHERE.density_gitm_right_before[index_longitude_lower+1][index_latitude_lower+1][index_altitude_lower+1] +
    (1-dist_lon) * dist_lat * PARAMS->ATMOSPHERE.density_gitm_right_before[index_longitude_lower][index_latitude_lower+1][index_altitude_lower+1];

  double density_interpolated_right_before;
  double a_right_before, b_right_before;
  a_right_before = ( log(density_lat_lon_higher_plane_gitm_right_before) - log(density_lat_lon_lower_plane_gitm_right_before) ) / ( PARAMS->ATMOSPHERE.altitude_gitm[0][0][index_altitude_lower+1] - PARAMS->ATMOSPHERE.altitude_gitm[0][0][index_altitude_lower] );
  b_right_before = log(density_lat_lon_lower_plane_gitm_right_before) - a_right_before*PARAMS->ATMOSPHERE.altitude_gitm[0][0][index_altitude_lower];
  density_interpolated_right_before = exp(a_right_before * altitude + b_right_before);
  //     printf("\n%e | %e | %e\n", density_lat_lon_lower_plane_gitm_right_before, density_interpolated_right_before, density_lat_lon_higher_plane_gitm_right_before);


  double density_lat_lon_lower_plane_gitm_right_after;
  density_lat_lon_lower_plane_gitm_right_after = (1-dist_lon) * (1-dist_lat) * PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower][index_latitude_lower][index_altitude_lower] +
    dist_lon * (1-dist_lat) * PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower+1][index_latitude_lower][index_altitude_lower] +
    dist_lon * dist_lat * PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower+1][index_latitude_lower+1][index_altitude_lower] +
    (1-dist_lon) * dist_lat * PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower][index_latitude_lower+1][index_altitude_lower];

  double density_lat_lon_higher_plane_gitm_right_after;
  density_lat_lon_higher_plane_gitm_right_after = (1-dist_lon) * (1-dist_lat) * PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower][index_latitude_lower][index_altitude_lower+1] +
    dist_lon * (1-dist_lat) * PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower+1][index_latitude_lower][index_altitude_lower+1] +
    dist_lon * dist_lat * PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower+1][index_latitude_lower+1][index_altitude_lower+1] +
    (1-dist_lon) * dist_lat * PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower][index_latitude_lower+1][index_altitude_lower+1];

  /* printf("\n%e | %e\n", density_lat_lon_lower_plane_gitm_right_after,density_lat_lon_higher_plane_gitm_right_after); */
  /* printf("%e | %e\n", PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower][index_latitude_lower][index_altitude_lower],PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower][index_latitude_lower][index_altitude_lower+1] ); */
  /* printf("%e | %e\n", PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower+1][index_latitude_lower][index_altitude_lower],PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower+1][index_latitude_lower][index_altitude_lower+1] ); */
  /* printf("%e | %e\n", PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower][index_latitude_lower+1][index_altitude_lower],PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower][index_latitude_lower+1][index_altitude_lower+1] ); */
  /* printf("%e | %e\n", PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower+1][index_latitude_lower+1][index_altitude_lower],PARAMS->ATMOSPHERE.density_gitm_right_after[index_longitude_lower+1][index_latitude_lower+1][index_altitude_lower+1] ); */


  double density_interpolated_right_after;
  double a_right_after, b_right_after;
  a_right_after = ( log(density_lat_lon_higher_plane_gitm_right_after) - log(density_lat_lon_lower_plane_gitm_right_after) ) / ( PARAMS->ATMOSPHERE.altitude_gitm[0][0][index_altitude_lower+1] - PARAMS->ATMOSPHERE.altitude_gitm[0][0][index_altitude_lower] );
  b_right_after = log(density_lat_lon_lower_plane_gitm_right_after) - a_right_after*PARAMS->ATMOSPHERE.altitude_gitm[0][0][index_altitude_lower];
  density_interpolated_right_after = exp(a_right_after * altitude + b_right_after);
  //     printf("%e | %e | %e\n", density_lat_lon_lower_plane_gitm_right_after, density_interpolated_right_after, density_lat_lon_higher_plane_gitm_right_after);


  /*   /\******************************************************************\/ */
  /*   /\***************** INTERPOLATION IN TIME **************************\/ */
  /*   /\******************************************************************\/ */
  double a_time;
  double b_time;
  a_time = ( density_interpolated_right_after - density_interpolated_right_before ) / (PARAMS->ATMOSPHERE.array_gitm_date[PARAMS->ATMOSPHERE.gitm_index_higher][0] - PARAMS->ATMOSPHERE.array_gitm_date[PARAMS->ATMOSPHERE.gitm_index_lower][0]) ;
  b_time = density_interpolated_right_before - a_time * PARAMS->ATMOSPHERE.array_gitm_date[PARAMS->ATMOSPHERE.gitm_index_lower][0];
  *density = a_time * et + b_time;

  /* printf("\n%s | %s | %s\n",name_gitm_file_right_before_epoch_of_sc, time_sc, name_gitm_file_right_after_epoch_of_sc ); */
  /*     printf("%e | %e | %e\n", density_interpolated_right_before, *density, density_interpolated_right_after); */


  PARAMS->ATMOSPHERE.is_first_step_of_run = 0;
  return 0;

  /*   // Free memory */
  /*   for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){ */
  /*     for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){ */
  /*       free(PARAMS->ATMOSPHERE.longitude_gitm[i][j]); */
  /*       } */
  /*       free(PARAMS->ATMOSPHERE.longitude_gitm[i]); */
  /*   } */
  /*   free(PARAMS->ATMOSPHERE.longitude_gitm); */
  /*   for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){ */
  /*     for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){ */
  /*       free(PARAMS->ATMOSPHERE.latitude_gitm[i][j]); */
  /*       } */
  /*   free(PARAMS->ATMOSPHERE.latitude_gitm[i]); */
  /*   } */
  /*   free(PARAMS->ATMOSPHERE.latitude_gitm ); */
  /*   for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){ */
  /*     for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){ */
  /*       free(PARAMS->ATMOSPHERE.altitude_gitm[i][j]); */
  /*       } */
  /*   free(PARAMS->ATMOSPHERE.altitude_gitm[i]); */
  /*   } */
  /*   free(PARAMS->ATMOSPHERE.altitude_gitm           ); */
  /*   for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){ */
  /*     for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm  ; j++){ */
  /*       free(PARAMS->ATMOSPHERE.density_gitm_right_before[i][j]); */
  /*       } */
  /*   free(PARAMS->ATMOSPHERE.density_gitm_right_before[i]); */
  /*   } */
  /*   free(PARAMS->ATMOSPHERE.density_gitm_right_before); */

  /*   for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){ */
  /*     for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){ */
  /*       free(PARAMS->ATMOSPHERE.longitude_gitm[i][j]); */
  /*       } */
  /*       free(PARAMS->ATMOSPHERE.longitude_gitm[i]); */
  /*   } */
  /*   free(PARAMS->ATMOSPHERE.longitude_gitm              ); */
  /*   for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){ */
  /*     for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){ */
  /*       free(PARAMS->ATMOSPHERE.latitude_gitm[i][j]); */
  /*       } */
  /*   free(PARAMS->ATMOSPHERE.latitude_gitm[i]); */
  /*   } */
  /*   free(PARAMS->ATMOSPHERE.latitude_gitm       ); */
  /*   for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){ */
  /*     for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){ */
  /*       free(PARAMS->ATMOSPHERE.altitude_gitm[i][j]); */
  /*       } */
  /*   free(PARAMS->ATMOSPHERE.altitude_gitm[i]); */
  /*   } */
  /*   free(PARAMS->ATMOSPHERE.altitude_gitm       ); */
  /*   for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm ; i++){ */
  /*     for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm ; j++){ */
  /*       free(PARAMS->ATMOSPHERE.density_gitm_right_after[i][j]); */
  /*       } */
  /*   free(PARAMS->ATMOSPHERE.density_gitm_right_after[i]); */
  /*   } */
  /*   free(PARAMS->ATMOSPHERE.density_gitm_right_after ); */










  /* if ((dir = opendir ("/raid3/Gitm/Runs/20150317/data")) != NULL) { */
  /*   while ((ent = readdir (dir)) != NULL) { */
  /*     // find the extension of each file. We only care about the file if it is a .bin file */
  /*     strcpy(extension_gitm_file, ""); */
  /* next = &(ent->d_name)[0]; */
  /*   find_extension = (int)(strchr(next, '.') - next); */
  /*   strncat(extension_gitm_file, next+find_extension,4); */

  /*   if (strcmp(extension_gitm_file, ".bin") == 0) { // We only care about the file if it is a .bin file */
      
  /*     // We only care about the file if its name starts with "3DALL_t" (see IMPORTANT assumption above) */

  /*     strcpy(start_filename, ""); */
  /*     strncat(start_filename, &(ent->d_name)[0], 7); */
  /*     if (strcmp(start_filename, "3DALL_t") == 0) { // We only care about the file if its name starts with "3 *\/DALL_t" (see IMPORTANT assumption above) */

  /* 	// Convert the date format of the GITM file name to the Julian Epoch format */
  /* 	strcpy(date_gitm_file, "20"); */
  /* 	strncat(date_gitm_file, &(ent->d_name)[7],2); */
  /* 	strcat(date_gitm_file, "-"); */
  /* 	strncat(date_gitm_file, &(ent->d_name)[9],2); */
  /* 	strcat(date_gitm_file, "-"); */
  /* 	strncat(date_gitm_file, &(ent->d_name)[11],2); */
  /* 	strcat(date_gitm_file, "T"); */
  /* 	strncat(date_gitm_file, &(ent->d_name)[14],2); */
  /* 	strcat(date_gitm_file, ":"); */
  /* 	strncat(date_gitm_file, &(ent->d_name)[16],2); */
  /* 	strcat(date_gitm_file, ":"); */
  /* 	strncat(date_gitm_file, &(ent->d_name)[18],2); */
  /* 	str2et_c(date_gitm_file, &et_gitm_file); */

  /* 	// Find the 2 GITM files that have dates around the epoch of the sc */
  /* 	delta_et = et_gitm_file - et; */
  /* 	if (delta_et < 0){ // the date of the GITM file is older than the epoch of the sc */
  /* 	  if (delta_et > delta_et_right_before) { // the date of the GITM file is newer than the previous record */
  /* 	    delta_et_right_before = delta_et; */
  /* 	    et_right_before = et_gitm_file; */
  /* 	    strcpy(name_gitm_file_right_before_epoch_of_sc, "/raid3/Gitm/Runs/20150317/data/"); */
  /* 	    strcat(name_gitm_file_right_before_epoch_of_sc, ent->d_name); */
  /* 	  } */
  /* 	} */
  /* 	if (delta_et > 0){ // the date of the GITM file is newer than the epoch of the sc */
  /* 	  if (delta_et < delta_et_right_after) { // the date of the GITM file is older than the previous record */
  /* 	    delta_et_right_after = delta_et; */
  /* 	    et_right_after = et_gitm_file; */
  /* 	    strcpy(name_gitm_file_right_after_epoch_of_sc, "/raid3/Gitm/Runs/20150317/data/"); */
  /* 	    strcat(name_gitm_file_right_after_epoch_of_sc, ent->d_name); */
  /* 	  } */
  /* 	} */
  /* 	if (delta_et == 0){ // the date of the GITM file is exactly the epoch of the sc */
  /* 	  do_not_need_to_interpolate = 1; */
  /* 	  strcpy(name_gitm_file_right_before_epoch_of_sc, "/raid3/Gitm/Runs/20150317/data/"); */
  /* 	  strcat(name_gitm_file_right_before_epoch_of_sc, ent->d_name); */
  /* 	} */

  /*     } // end of we only care about the file if its name starts with "3DALL_t" (see IMPORTANT assumption above)   */
  /*   } // end of we only care about the file if it is a .bin file */
   
  /*   } */
  /*   closedir (dir); */
  /* } else { */
  /*   /\* could not open directory *\/ */
  /*   printf("Could not open the GITM data directory. The program will stop.\n"); */
  /*   exit(0); */
  /* } */





  //printf("\n\n%s | %s | %s \n\n", name_gitm_file_right_before_epoch_of_sc, time_sc, name_gitm_file_right_after_epoch_of_sc);
  //  char time_gitm_right_before_sc[256], time_gitm_right_after_sc[256];
  /* et2utc_c(et_right_before, "ISOC" ,0 ,255 , time_gitm_right_before_sc);     */
  /* et2utc_c(et_right_after, "ISOC" ,0 ,255 , time_gitm_right_after_sc);     */
  /* printf("\n\n%s | %s | %s || %d\n\n", time_gitm_right_before_sc, time_sc, time_gitm_right_after_sc, do_not_need_to_interpolate  ); */

}




/* int eclipse_sun_moon_sc(){ */






/* /\*   //  spkez_c(10, et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. *\/ */
/* /\*   double x_earth_sun[6], x_moon_sun[6]; *\/ */
/* /\*   double lt_earth_sun, lt_moon_sun; *\/ */
/* /\*     double r_earth2sun_J2000[3], r_moon2sun_J2000[3]; *\/ */
/* /\*     double moon_radius = 1737. km; *\/ */
/* /\*   spkez_c(10, et, "J2000", "NONE", 399, x_earth_sun, &lt_earth_sun); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. *\/ */
/* /\*   spkez_c(10, et, "J2000", "NONE", 301, x_moon_sun, &lt_moon_sun); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. *\/ */
/* /\*   r_earth2sun_J2000[0] = x_earth_sun[0]; *\/ */
/* /\*   r_earth2sun_J2000[1] = x_earth_sun[1]; *\/ */
/* /\*   r_earth2sun_J2000[2] = x_earth_sun[2]; *\/ */
  
/* /\*   r_moon2sun_J2000[0] = x_moon_sun[0]; *\/ */
/* /\*   r_moon2sun_J2000[1] = x_moon_sun[1]; *\/ */
/* /\*   r_moon2sun_J2000[2] = x_moon_sun[2]; *\/ */

/* /\*   r_edge_moon3sun_J2000[0] = jj; *\/ */
  

/* /\*   // sc in shadow of Moon if: *\/ */
/* /\*   // - angle (sun to moon, sun to sc) < angle (sun to moon, sun to edge of moon); and *\/ */
/* /\*   // - distance sun to moon < distance sun to sc *\/ */


/*   // sc in eclipse by moon if: */
/*   // - sc in shadow of Moon; and */
/*   // - sc not in shadow of Earth */
  
/* /\*   // Returns if the satellite is or not in the shadow of Earth *\/ */
/* /\*   shadow_light( SC->INTEGRATOR.shadow, SC->r_i2cg_INRTL, SC->et, PARAMS); *\/ */

/*   // Returns if the satellite is or not in the shadow of Moon */
/*   shadow_light_moon( SC->INTEGRATOR.shadow, SC->r_i2cg_INRTL, SC->et, PARAMS); */
  
/*   return 0; */

/* } */





/* ///////////////////////////////////////////////////////////////////////////////////////// */
/* // */
/* //  Name:           read_gitm_bin */
/* //  Purpose:        Read a GITM bin file   */
/* //  Assumptions:    None. */
/* //  References      Various. Note: the full code is in PropSim_v3/idl_propagator/read_gitm_bin.c */
/* // */
/* //  Change Log: */
/* //      |   Developer   |       Date    |   SCR     |   Notes */
/* //      | --------------|---------------|-----------|------------------------------- */
/* //      | C. Bussy-Virat| 01/24/2016    |   ---     | Initial implementation */
/* // */
/* ///////////////////////////////////////////////////////////////////////////////////////// */
/* int read_gitm_bin( double small_array[5], double ***longitude_gitm, double ***latitude_gitm, double ***altitude_gitm, double ***density_gitm, char time_gitm[256], int nLons, int nLats, int nAlts, char name_gitm_file_to_read[256]) { */

/*   // Declarations */
/*      int pos_for_lon, iDataLength, iHeaderLength; */
/*    int time_gitm_temp[7]; */
/*      int j,k; */
/*   int i; */
/*   int nVars = 0; */
/*   FILE *gitm_file; */

/*   //  small_array = malloc(sizeof(double) * 3); */
/*   small_array[0] = 0; */
/*   small_array[1] = 1; */
/*   small_array[2] = 2; */

/*   gitm_file = fopen(name_gitm_file_to_read,"rb"); */

/* // NLONS, NLATS, NALTS */
/* fseek(gitm_file, 20, SEEK_SET	 ); */
/*   fread(&nLons,sizeof(nLons),1,gitm_file); */
/*   fread(&nLats,sizeof(nLats),1,gitm_file); */
/*   fread(&nAlts,sizeof(nAlts),1,gitm_file); */
/* fseek(gitm_file, 4, SEEK_CUR	 ); */
/*  /\*   printf("nLons = %d \n", nLons); *\/ */
/*  /\* printf("nLats = %d \n", nLats); *\/ */
/*  /\* printf("nAlts = %d \n", nAlts); *\/ */

/* // NVARS */
/* fseek(gitm_file, 4, SEEK_CUR	 ); */
/*   fread(&nVars,sizeof(nVars),1,gitm_file); */
/* fseek(gitm_file, 4, SEEK_CUR	 ); */
/* // printf("nVars = %d \n", nVars); */
  
/*    // TIME */
/*    fseek(gitm_file, 4+nVars*(4+4+40), SEEK_CUR	 );     */
/*   fread(&time_gitm_temp,sizeof(time_gitm_temp),1,gitm_file); */
/*    /\* for (i = 0; i < 7; i++){ *\/ */
/*    /\*    printf("time: %d\n", time_gitm_temp[i]); *\/ */
/*    /\*   } *\/ */

/*      // SKIP VARIABLES BEFORE THE DENSITY */
/*      iHeaderLength = 8L + 4+4 +	3*4 + 4+4 + 4 + 4+4 + nVars*40 + nVars*(4+4) +  7*4 + 4+4; */
/*      iDataLength = nLons*nLats*nAlts*8 + 4+4; */
/*      pos_for_lon = iHeaderLength  ;// density is the 4th variable in the GITM file */
/*      fseek(gitm_file, pos_for_lon, SEEK_SET); */



/*      longitude_gitm = malloc(  nLons * sizeof(double **) ); */
/*       for (i = 0; i < nLons; i++){ */
/* 	longitude_gitm[i] = malloc(nLats * sizeof(double *)); */
/* 	for (j = 0; j < nLats; j++){ */
/* 	  longitude_gitm[i][j] = malloc(nAlts * sizeof(double)); */
/* 	} */
/*       }  */
/*      latitude_gitm = malloc(  nLons * sizeof(double **) ); */
/*       for (i = 0; i < nLons; i++){ */
/* 	latitude_gitm[i] = malloc(nLats * sizeof(double *)); */
/* 	for (j = 0; j < nLats; j++){ */
/* 	  latitude_gitm[i][j] = malloc(nAlts * sizeof(double)); */
/* 	} */
/*       } */
/*       altitude_gitm = malloc(  nLons * sizeof(double **) ); */
/*       for (i = 0; i < nLons; i++){ */
/* 	altitude_gitm[i] = malloc(nLats * sizeof(double *)); */
/* 	for (j = 0; j < nLats; j++){ */
/* 	  altitude_gitm[i][j] = malloc(nAlts * sizeof(double)); */
/* 	} */
/*       } */


/*       density_gitm = malloc(  nLons * sizeof(double **) ); */
/*       for (i = 0; i < nLons; i++){ */
/* 	density_gitm[i] = malloc(nLats * sizeof(double *)); */
/* 	for (j = 0; j < nLats; j++){ */
/* 	  density_gitm[i][j] = malloc(nAlts * sizeof(double)); */
/* 	} */
/*       } */
/*       if (  longitude_gitm == NULL ){ */
/*   	printf("Could not allow memory for longitude_gitm. The program will stop.\n"); */
/*   	exit(0); */
/*       } */
/*       if (  latitude_gitm == NULL ){ */
/*   	printf("Could not allow memory for latitude_gitm. The program will stop.\n"); */
/*   	exit(0); */
/*       } */
/*       if (  altitude_gitm == NULL ){ */
/*   	printf("Could not allow memory for altitude_gitm. The program will stop.\n"); */
/*   	exit(0); */
/*       } */

/*        if (  density_gitm == NULL ){ */
/*   	printf("Could not allow memory for density_gitm. The program will stop.\n"); */
/*   	exit(0); */
/*       } */

/*   // READ THE LONGITUDE */
/*      fseek(gitm_file, 4, SEEK_CUR ); */
/*      for (k = 0; k < nAlts; k++){ */
/*        for (j = 0; j < nLats; j++){ */
/*   	 for (i = 0; i < nLons; i++){ */
/*      	   fread(&longitude_gitm[i][j][k], sizeof(longitude_gitm[i][j][k]), 1, gitm_file); */
/*      	 } */
/*        } */
/*      } */
/*      fseek(gitm_file, 4, SEEK_CUR ); */
/*      //    printf("longitude = %f\n", longitude_gitm[11][40][21]); */


/*   // READ THE LATITUDE */
/*      fseek(gitm_file, 4, SEEK_CUR ); */
/*      for (k = 0; k < nAlts; k++){ */
/*        for (j = 0; j < nLats; j++){ */
/*   	 for (i = 0; i < nLons; i++){ */
/*      	   fread(&latitude_gitm[i][j][k], sizeof(latitude_gitm[i][j][k]), 1, gitm_file); */
/*      	 } */
/*        } */
/*      } */
/*      fseek(gitm_file, 4, SEEK_CUR ); */
/*      //     printf("latitude = %f\n", latitude_gitm[11][40][21]); */

/*   // READ THE ALTITUDE */

/*      fseek(gitm_file, 4, SEEK_CUR ); */
/*      for (k = 0; k < nAlts; k++){ */
/*        for (j = 0; j < nLats; j++){ */
/*   	 for (i = 0; i < nLons; i++){ */
/*      	   fread(&altitude_gitm[i][j][k], sizeof(altitude_gitm[i][j][k]), 1, gitm_file); */
/*      	 } */
/*        } */
/*      } */
/*      fseek(gitm_file, 4, SEEK_CUR ); */
/*      //   printf("altitude = %f\n", altitude_gitm[11][40][21]); */

/*      // READ THE DENSITY AT THE GIVEN LON/LAT/ALT */

/*      fseek(gitm_file, 4, SEEK_CUR ); */
/*      for (k = 0; k < nAlts; k++){ */
/*        for (j = 0; j < nLats; j++){ */
/*   	 for (i = 0; i < nLons; i++){ */
/*      	   fread(&density_gitm[i][j][k], sizeof(density_gitm[i][j][k]), 1, gitm_file); */
/*      	 } */
/*        } */
/*      } */
/*      //               printf("density = %e\n", density_gitm[11][40][21]); */

/* /\* int read_gitm_bin( double *longitude_gitm, double *latitude_gitm, double *altitude_gitm, double *density_gitm, char time_gitm[256], char name_gitm_file_to_read[256]) { *\/ */
/*   return 0; */
/* } */


int calculate_cd(double *cd,
		 double acco_coeff, 
		 double v_sc_mag, //in km/s
		 double Ta, // atmospheric temperature in K, from NRLSMSIS
		 double surface_area, // in km^2
		 double gamma,
		 double Ma // in kg/mol!!!!!from NRLMSIS?
		 ){
  //Ta = 800;  // !!!!!!!!!!!!!!!!! REMOVE
  //v_sc_mag  = 7.3; // !!!!!!!!!!!!!!!!! REMOVE

//printf("%e\n",surface_area*gamma);
  double Tw = 300.; // temperature of the satellite surface // assumed at 300 K in Moe04, 273 in Sutton07 and Bruisma03
  v_sc_mag = v_sc_mag * 1000.; // km/s to m/s
  surface_area = surface_area * 1000000.; // km^2 to m^2
  //    printf("%e\n", v_sc_mag);

      double A_ref =  surface_area * gamma;
      double R = 8.31; // universal gas constant in J/mol.K
	double s = v_sc_mag / sqrt( 2 * R * Ta / Ma) ;
	double Q = 1 + 1 / ( 2 * s*s );
	double T = Ma * v_sc_mag * v_sc_mag / ( 3 * R );  // kinetic temperature
				      
	double Vr = v_sc_mag * sqrt( 2./3 * ( 1 + acco_coeff * ( Tw / T - 1 )  ) );
	
	double P = exp( -( gamma*gamma ) * ( s*s ) ) / s;
	double Z = 1 + erf( gamma * s );

	double    term1 = P / sqrt(M_PI);
	double term2 = gamma * Q * Z;
	double fac1 = gamma * Vr / (2 * v_sc_mag);
	  double fac2 = gamma * sqrt(M_PI) * Z + P;
	  double term3 = fac1 * fac2;
	  double term = term1 + term2 + term3;
	  cd[0] = surface_area / A_ref * term ;


	  //	  printf("%e %e %e\n",term1 , term2 , term3);
/*     // TOTAL NORMALIZED DRAG COEFFICIENT */
/*     A_ref_tot = 0 */
/*     for isurf in range(nb_surf): */
/*   if gamma[isurf] >= 0: */
/*   A_ref_tot = A_ref_tot + A_ref[isurf] */
/*     cd_tot_norm = np.sum(cd * A_ref) / A_ref_tot //     see sutton09a equation 9 */

//	  printf("%f\n", cd[0]);
    return 0;

}



// same as caclulate_cd but inputs are a bit different because used a 3d model (opengl)
int calculate_cd_opengl(double *cd,
		 double acco_coeff, 
		 double v_sc_mag, //in km/s
		 double Ta, // atmospheric temperature in K, from NRLSMSIS
		 double surface_area_projected, // in km^2
		 double gamma,
		 double Ma // in kg/mol!!!!!from NRLMSIS?
		 ){
  //Ta = 800;  // !!!!!!!!!!!!!!!!! REMOVE
  //v_sc_mag  = 7.3; // !!!!!!!!!!!!!!!!! REMOVE

//printf("%e\n",surface_area*gamma);
  double Tw = 300.; // temperature of the satellite surface // assumed at 300 K in Moe04, 273 in Sutton07 and Bruisma03
  v_sc_mag = v_sc_mag * 1000.; // km/s to m/s
  surface_area_projected = surface_area_projected * 1000000.; // km^2 to m^2
  //    printf("%e\n", v_sc_mag);


      double R = 8.31; // universal gas constant in J/mol.K
	double s = v_sc_mag / sqrt( 2 * R * Ta / Ma) ;
	double Q = 1 + 1 / ( 2 * s*s );
	double T = Ma * v_sc_mag * v_sc_mag / ( 3 * R );  // kinetic temperature
				      
	double Vr = v_sc_mag * sqrt( 2./3 * ( 1 + acco_coeff * ( Tw / T - 1 )  ) );
	
	double P = exp( -( gamma*gamma ) * ( s*s ) ) / s;
	double Z = 1 + erf( gamma * s );

	double    term1 = P / sqrt(M_PI);
	double term2 = gamma * Q * Z;
	double fac1 = gamma * Vr / (2 * v_sc_mag);
	  double fac2 = gamma * sqrt(M_PI) * Z + P;
	  double term3 = fac1 * fac2;
	  double term = term1 + term2 + term3;
	  cd[0] = term / gamma;


	  //	  printf("%e %e %e\n",term1 , term2 , term3);
/*     // TOTAL NORMALIZED DRAG COEFFICIENT */
/*     A_ref_tot = 0 */
/*     for isurf in range(nb_surf): */
/*   if gamma[isurf] >= 0: */
/*   A_ref_tot = A_ref_tot + A_ref[isurf] */
/*     cd_tot_norm = np.sum(cd * A_ref) / A_ref_tot //     see sutton09a equation 9 */

//	  printf("%f\n", cd[0]);
    return 0;

}


  double factorial(unsigned long f)
  {
    if ( f == 0 ) 
      return 1;
      return(f * factorial(f - 1));
	}



/* Assumptions: */
/* - the zenith angle of the satellite is approximated to be the same as the angle between the direction "center of the Earth to satelltie" and the direction "center of the Earth to Sun" (instead of the angle between the direction "center of the Earth to satelltie" and the direction "satelltie to Sun"). This approximation is valid since the distance Earth to Sun is much larger than the distance Earth to satellite */
/* - Similarly, the zenith angle of the Earth element is approximated to be the same as the angle between the direction "center of the Earth to Earth element" and the direction "center of the Earth to Sun" (instead of the angle between the direction "center of the Earth to Earth element" and the direction "Earth element to Sun"). */
/* - albedo and infra-red emmissivities are uniform (don't depend on latitude) */
int build_earth_pressure_map(GRAVITY_T  *Gravity, int iProc){

  /*  Notations: */
  /* - O center of the Earth */
  /* - C position of the satellite */
  /* - H projection of the satellite lcoation on the surface of the Earth (ie sub-satellite point) */
  /* - M position of the Earth surface elements */
  /* - S position of the Sun */
  /* - P: projection of M on the direction OC */
  /* - nvec: normal to the satellite surface in the reference frame defined as: */
  // - the z vector is the direction OH
  // - the y vector is the vector in the plane OHS and in the direction of the Sun
  // - the x vector completes the orthogonal basis
  // Declarations
  double mp, elt_area;
  int ccc= 0;
  double a_earth_pressure_fac[3];
    double a_earth_pressure_fac_ir[3];
  double nvec[3];
  	  double cm[3], oc[3], om[3];
	  double cos_nvec_cm, cm_mag;
	  double albedo = 0.34; 
	  double emiss = 0.68;
    double max_radius = Gravity->max_radius_map; // max radius of all sc (km)
    double min_radius = Gravity->min_radius_map; // min radius of all sc (km)
    double dradius, radius;
    int nradius = Gravity->nradius_map;
    int iradius;

    double max_zenith = Gravity->max_zenith_map; 
    double min_zenith = Gravity->min_zenith_map; 
    double dzenith, zenith;
    int nzenith = Gravity->nzenith_map;
    int izenith;

    double max_elev_elt = Gravity->max_elev_elt_map; 
    double min_elev_elt = Gravity->min_elev_elt_map; 
    double delev_elt, elev_elt;
    int nelev_elt = Gravity->nelev_elt_map;
    int ielev_elt;

    double max_azim_elt = Gravity->max_azim_elt_map; 
    double min_azim_elt = Gravity->min_azim_elt_map; 
    double dazim_elt, azim_elt;
    int nazim_elt = Gravity->nazim_elt_map;
    int iazim_elt;

    double max_elev_surf = Gravity->max_elev_surf_map; 
    double min_elev_surf = Gravity->min_elev_surf_map; 
    double delev_surf, elev_surf;
    int nelev_surf = Gravity->nelev_surf_map;
    int ielev_surf;

    double max_azim_surf = Gravity->max_azim_surf_map; 
    double min_azim_surf = Gravity->min_azim_surf_map; 
    double dazim_surf, azim_surf;
    int nazim_surf = Gravity->nazim_surf_map;
    int iazim_surf;

    double om_norm[3], os[3];
    double osmag = 149597871. ;// average distance Earth-Sun
    double zenith_elt, cos_zenith_elt;

    //		    printf("nradius %d, nzenith %d, nelev_surf %d, nazim_surf %d, nelev_elt %d, nazim_elt %d | max elt %f %f | total size %.1fM\n", nradius, nzenith, nelev_surf, nazim_surf, nelev_elt, nazim_elt, Gravity->max_elev_elt_map, Gravity->max_elev_surf_map, nradius*nzenith*nelev_surf*nazim_surf*3*4./1024./1024.);
    		    printf("nradius %d, nzenith %d, nelev_surf %d, nazim_surf %d, nelev_elt %d, nazim_elt %d | max elt %f %f | total size %.1fM\n", nradius, nzenith, nelev_surf, nazim_surf, nelev_elt, nazim_elt, Gravity->max_elev_elt_map, Gravity->max_elev_surf_map, nradius*nzenith*nelev_surf*nazim_surf*3*8./1024./1024.);
    		    printf("Map: max elevation elt %.1f, min elevation surface normal %.1f | total size %.1fM\n",  Gravity->max_elev_elt_map, Gravity->min_elev_surf_map, nradius*nzenith*nelev_surf*nazim_surf*3*4./1024./1024.);
  Gravity->file_earth_pressure_map = fopen(Gravity->filename_earth_pressure_map, "wb");		    
    for (iradius = 0; iradius < nradius; iradius++){ // go over radii (distance center of the Earth to satellite)


      	if (iradius == nradius - 1){
	  dradius = max_radius - (min_radius + (nradius - 2)*Gravity->dradius_map);//not used 
	  radius = max_radius;
	}
	else{
	  dradius =  Gravity->dradius_map;
	  radius = min_radius + iradius * dradius;
	}
	for (izenith = 0; izenith < nzenith; izenith++){ // go over zeniths	  
	  if (izenith == nzenith - 1){
	    dzenith = max_zenith - (min_zenith + (nzenith - 2)*Gravity->dzenith_map);//not used 
	    zenith = max_zenith;
	  }
	  else{
	    dzenith =  Gravity->dzenith_map;
	    zenith = min_zenith + izenith * dzenith;
	  }
	  zenith = zenith * M_PI / 180.;
	for (ielev_surf = 0; ielev_surf < nelev_surf; ielev_surf++){ // go over all "elevations" of satellite normal vector. See the definition of this "elevation" in comments a bit below
	  
	  if (ielev_surf == nelev_surf - 1){
	    delev_surf = max_elev_surf - (min_elev_surf + (nelev_surf - 2)*Gravity->delev_surf_map);//not used 
	    elev_surf = max_elev_surf;
	  }
	  else{
	    delev_surf =  Gravity->delev_surf_map;
	    elev_surf = min_elev_surf + ielev_surf * delev_surf             ;
	  }

	  elev_surf = -180;//!!!!!!!!!!!!!!!!!!! remove line
	  elev_surf = elev_surf * M_PI / 180.;
	for (iazim_surf = 0; iazim_surf < nazim_surf; iazim_surf++){ // go over all "azimuth" of satellite normal vector. See the definition of this "azimuth" in comments a bit below

	  
	  if (iazim_surf == nazim_surf - 1){
	    dazim_surf = max_azim_surf - (min_azim_surf + (nazim_surf - 2)*Gravity->dazim_surf_map);//not used 
	    azim_surf = max_azim_surf;
	  }
	  else{
	    dazim_surf =  Gravity->dazim_surf_map;
	    azim_surf = min_azim_surf + iazim_surf * dazim_surf             ;
	  }
	  azim_surf = azim_surf * M_PI / 180.;
	  // the "elevation" of the satellite normal vector is the angle between nvec and CO (direction satellite to center of the Earth, ie -z direction where z has been previously defined)

	  a_earth_pressure_fac[0] = 0; a_earth_pressure_fac[1] = 0; a_earth_pressure_fac[2] = 0;
	  a_earth_pressure_fac_ir[0] = 0; a_earth_pressure_fac_ir[1] = 0; a_earth_pressure_fac_ir[2] = 0;
	  double om_sum[3], cm_sum[3];
	  om_sum[0] = 0; om_sum[1] = 0; om_sum[2] = 0;
	  cm_sum[0] = 0; cm_sum[1] = 0; cm_sum[2] = 0;
	for (ielev_elt = 0; ielev_elt < nelev_elt; ielev_elt++){ // go over all "elevations" of Earth element. See the definition of this "elevation" in comments a bit below
	  if (ielev_elt == nelev_elt - 1){
	    delev_elt = max_elev_elt - (min_elev_elt + (nelev_elt - 2)*Gravity->delev_elt_map);//not used 
	    elev_elt = max_elev_elt;
	  }
	  else{
	    delev_elt =  Gravity->delev_elt_map;
	    elev_elt = min_elev_elt + ielev_elt * delev_elt             ;
	  }
	  elev_elt = elev_elt * M_PI / 180.;
	for (iazim_elt = 0; iazim_elt < nazim_elt; iazim_elt++){ // go over all "azimuth" of Earth element. See the definition of this "azimuth" in comments a bit below
	  // COMMENT block below because we wawnt to exlude 360 otherwise the pressure from the elemetn at 0 deg and the element at 360 deg are both counted, while they correspond to a single and same element
	  /* if (iazim_elt == nazim_elt - 1){ */
	  /*   dazim_elt = max_azim_elt - (min_azim_elt + (nazim_elt - 2)*Gravity->dazim_elt_map);//not used  */
	  /*   azim_elt = max_azim_elt; */
	  /* } */
	  /* else{ */
	    dazim_elt =  Gravity->dazim_elt_map;
	    azim_elt = min_azim_elt + iazim_elt * dazim_elt             ;
	    //	  }
	  azim_elt = azim_elt * M_PI / 180.;
	  nvec[0] = sin(elev_surf) * cos (azim_surf);
	  nvec[1] = sin(elev_surf) * sin (azim_surf);
	  nvec[2] = cos(elev_surf);

	  // Compute the solar zenith angle of the Earth element as a function of the solar zenith angle of the satellite
	  // The reference frame is defined in the comments at the beginning of this function
	  // The "elevation" of the surface element (elev_elt) is defined as the angle (OH, OM).
	  // The "elevation" of the surface element (elev_elt) is defined as the angle (OH, OM), ie the angle (z, OM) (where z is previously defined)
	  // The "azimuth" of the surface element (azim_elt) is defined as the angle (x, OM_proj), where M_proj is the projection of M on the xy plane previously defined.

	  //	  if ((ielev_surf == 0) && (iazim_surf == 0)){ // only need to calculate the zenith angle of the Earth element for the first surface since the zenith angle of the Earth element is not dependent on the satellite surface
	  om_norm[0] = sin(elev_elt) * cos (azim_elt);
	  om_norm[1] = sin(elev_elt) * sin (azim_elt);
	  om_norm[2] = cos (elev_elt);
	  os[0] = 0;
	  os[1] = osmag * sin(zenith);
	  os[2] = osmag * cos(zenith);
	  v_dot(&cos_zenith_elt, om_norm, os);
	  cos_zenith_elt = cos_zenith_elt / osmag;
	  //	  zenith_elt = acos(cos_zenith_elt); // doesn't matter if you get it right by + or minus - (ie if it's 40 deg instead of -40 deg, it's the same in terms of results for the earth radiation pressure).
	  //	  }
	  oc[0] = 0; oc[1] = 0; oc[2] = radius; // by definiton of the reference frame (see comments at the beginning of the script), the satellite (C) is along the z direction
	  v_scale(om, om_norm, Gravity->radius); // assume spherical Earth
	  v_sub(cm, om, oc);
	  v_add(om_sum, om_sum, om);
	  v_add(cm_sum, cm_sum, cm);
	  v_mag(&cm_mag, cm);
	  v_dot(&cos_nvec_cm, cm, nvec);
	  cos_nvec_cm = cos_nvec_cm  / cm_mag;
	  if (cos_nvec_cm > 0){ // if the satellite surface sees the Earth element
	    mp = Gravity->radius * sin(elev_elt);
	    elt_area = Gravity->radius * delev_elt * M_PI/180. * mp * dazim_elt * M_PI/180.;
	    if (cos_zenith_elt > 0){ // if the Earth element sees the Sun
	    a_earth_pressure_fac[0] = a_earth_pressure_fac[0] + albedo * cos_zenith_elt * cos_nvec_cm * elt_area / (M_PI * cm_mag * cm_mag) * cm[0] / cm_mag;
	    a_earth_pressure_fac[1] = a_earth_pressure_fac[1] + albedo * cos_zenith_elt * cos_nvec_cm * elt_area / (M_PI * cm_mag * cm_mag) * cm[1] / cm_mag;
	    a_earth_pressure_fac[2] = a_earth_pressure_fac[2] + albedo * cos_zenith_elt * cos_nvec_cm * elt_area / (M_PI * cm_mag * cm_mag) * cm[2] / cm_mag;
	    } // end of if the Earth element sees the Sun
	    a_earth_pressure_fac[0] = a_earth_pressure_fac[0] + 0.25 * emiss * cos_nvec_cm * elt_area / (M_PI * cm_mag * cm_mag) * cm[0] / cm_mag;
	    a_earth_pressure_fac[1] = a_earth_pressure_fac[1] + 0.25 * emiss * cos_nvec_cm * elt_area / (M_PI * cm_mag * cm_mag) * cm[1] / cm_mag;
	    a_earth_pressure_fac[2] = a_earth_pressure_fac[2] + 0.25 * emiss * cos_nvec_cm * elt_area / (M_PI * cm_mag * cm_mag) * cm[2] / cm_mag;

	    	    a_earth_pressure_fac_ir[0] = a_earth_pressure_fac_ir[0] + 0.25 * emiss * cos_nvec_cm * elt_area / (M_PI * cm_mag * cm_mag) * cm[0] / cm_mag;
	    a_earth_pressure_fac_ir[1] = a_earth_pressure_fac_ir[1] + 0.25 * emiss * cos_nvec_cm * elt_area / (M_PI * cm_mag * cm_mag) * cm[1] / cm_mag;
	    a_earth_pressure_fac_ir[2] = a_earth_pressure_fac_ir[2] + 0.25 * emiss * cos_nvec_cm * elt_area / (M_PI * cm_mag * cm_mag) * cm[2] / cm_mag;

	  } // end of if the satellite surface sees the Earth element  
	} // go over all azimuth of Earth element	
	
	} // go over all elevation of Earth element
		  	fwrite(&(a_earth_pressure_fac[0]), sizeof(a_earth_pressure_fac[0]), 1, Gravity->file_earth_pressure_map);
		fwrite(&(a_earth_pressure_fac[1]), sizeof(a_earth_pressure_fac[1]), 1, Gravity->file_earth_pressure_map);
		fwrite(&(a_earth_pressure_fac[2]), sizeof(a_earth_pressure_fac[2]), 1, Gravity->file_earth_pressure_map); 
		//		v_print(a_earth_pressure_fac_ir, "a_earth_pressure_fac_ir");// should be only along the 3rd coordinate since the IR pressure is purely radial
	} // go over all azimuth of satellite normal vector

	} // go over all elevation of satellite normal vector


	
	}// go over zeniths
	
  } // go over radii
    //    printf("%d\n",ccc,  ccc / 1024./1024);
    fclose(Gravity->file_earth_pressure_map);
    exitf();
  return 0;
}
/* Assumptions: */
/* - the zenith angle of the satellite is approximated to be the same as the angle between the direction "center of the Earth to satelltie" and the direction "center of the Earth to Sun" (instead of the angle between the direction "center of the Earth to satelltie" and the direction "satelltie to Sun"). This approximation is valid since the distance Earth to Sun is much larger than the distance Earth to satellite */
/* - Similarly, the zenith angle of the Earth element is approximated to be the same as the angle between the direction "center of the Earth to Earth element" and the direction "center of the Earth to Sun" (instead of the angle between the direction "center of the Earth to Earth element" and the direction "Earth element to Sun"). */
/* - albedo and infra-red emmissivities are uniform (don't depend on latitude) */
int read_earth_pressure_map(GRAVITY_T  *Gravity, int iProc){

  /*  Notations: */
  /* - O center of the Earth */
  /* - C position of the satellite */
  /* - H projection of the satellite lcoation on the surface of the Earth (ie sub-satellite point) */
  /* - M position of the Earth surface elements */
  /* - S position of the Sun */
  /* - P: projection of M on the direction OC */
  /* - nvec: normal to the satellite surface in the reference frame defined as: */
  // - the z vector is the direction OH
  // - the y vector is the vector in the plane OHS and in the direction of the Sun
  // - the x vector completes the orthogonal basis
  // Declarations

  double nvec[3];
  	  double cm[3], oc[3], om[3];
	  double cos_nvec_cm, cm_mag;
	  double albedo = 0.34; 
	  double emiss = 0.68;
    double max_radius = Gravity->max_radius_map; // max radius of all sc (km)
    double min_radius = Gravity->min_radius_map; // min radius of all sc (km)
    double dradius, radius;
    int nradius = Gravity->nradius_map;
    int iradius;

    double max_zenith = Gravity->max_zenith_map; 
    double min_zenith = Gravity->min_zenith_map; 
    double dzenith, zenith;
    int nzenith = Gravity->nzenith_map;
    int izenith;

    double max_elev_elt = Gravity->max_elev_elt_map; 
    double min_elev_elt = Gravity->min_elev_elt_map; 
    double delev_elt, elev_elt;
    int nelev_elt = Gravity->nelev_elt_map;
    int ielev_elt;

    double max_azim_elt = Gravity->max_azim_elt_map; 
    double min_azim_elt = Gravity->min_azim_elt_map; 
    double dazim_elt, azim_elt;
    int nazim_elt = Gravity->nazim_elt_map;
    int iazim_elt;

    double max_elev_surf = Gravity->max_elev_surf_map; 
    double min_elev_surf = Gravity->min_elev_surf_map; 
    double delev_surf, elev_surf;
    int nelev_surf = Gravity->nelev_surf_map;
    int ielev_surf;

    double max_azim_surf = Gravity->max_azim_surf_map; 
    double min_azim_surf = Gravity->min_azim_surf_map; 
    double dazim_surf, azim_surf;
    int nazim_surf = Gravity->nazim_surf_map;
    int iazim_surf;

    double om_norm[3], os[3];
    double osmag = 149597871. ;// average distance Earth-Sun
    double zenith_elt, cos_zenith_elt;
		    FILE *file_earth_pressure_map = NULL;
    
		    printf("nradius %d, nzenith %d, nelev_surf %d, nazim_surf %d, nelev_elt %d, nazim_elt %d | max elt %f %f | total size %.1fM\n", nradius, nzenith, nelev_surf, nazim_surf, nelev_elt, nazim_elt, Gravity->max_elev_elt_map, Gravity->max_elev_surf_map, nradius*nzenith*nelev_surf*nazim_surf*3*8./1024./1024.);

      file_earth_pressure_map = fopen(Gravity->filename_earth_pressure_map, "r");

  if (iProc == 0){
  if (file_earth_pressure_map == NULL){
    printf("!***!\nThe 4D Earth pressure map file:\n%s\ncould not be found. The program will stop.\n!***!\n", Gravity->filename_earth_pressure_map); MPI_Finalize(); exit(0);
  }
  }

      
  Gravity->earth_pressure_map = malloc(nradius  * sizeof(double ****) );
  Gravity->radius_map = malloc(nradius  * sizeof(double) );
    Gravity->zenith_map = malloc(nzenith  * sizeof(double) );
      Gravity->elev_surf_map = malloc(nelev_surf  * sizeof(double) );
        Gravity->azim_surf_map = malloc(nazim_surf  * sizeof(double) );
     
  if ( Gravity->earth_pressure_map == NULL){
    printf("***! Could not allow memory to Gravity->earth_pressure_map. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  if ( Gravity->radius_map == NULL){
    printf("***! Could not allow memory to Gravity->radius_map. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  if ( Gravity->zenith_map == NULL){
    printf("***! Could not allow memory to Gravity->zenith_map. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  if ( Gravity->elev_surf_map == NULL){
    printf("***! Could not allow memory to Gravity->elev_surf_map. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  if ( Gravity->azim_surf_map == NULL){
    printf("***! Could not allow memory to Gravity->azim_surf_map. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }

    for (iradius = 0; iradius < nradius; iradius++){ // go over radii (distance center of the Earth to satellite)
      Gravity->earth_pressure_map[iradius] = malloc(nzenith * sizeof(double ***));
          if ( Gravity->earth_pressure_map[iradius] == NULL){
      printf("***! Could not allow memory to Gravity->earth_pressure_map[iradius]. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }

  
      	if (iradius == nradius - 1){
	  dradius = max_radius - (min_radius + (nradius - 2)*Gravity->dradius_map);//not used 
	  radius = max_radius;
	}
	else{
	  dradius =  Gravity->dradius_map;
	  radius = min_radius + iradius * dradius;
	}

		Gravity->radius_map[iradius] = radius;
	for (izenith = 0; izenith < nzenith; izenith++){ // go over zeniths

	  	  Gravity->earth_pressure_map[iradius][izenith] = malloc(nelev_surf  * sizeof(double**) );
	  if ( Gravity->earth_pressure_map[iradius][izenith] == NULL){ 
	    printf("***! Could not allow memory to Gravity->earth_pressure_map[iradius][izenith]. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
	  }

	  if (izenith == nzenith - 1){
	    dzenith = max_zenith - (min_zenith + (nzenith - 2)*Gravity->dzenith_map);//not used 
	    zenith = max_zenith;
	  }
	  else{
	    dzenith =  Gravity->dzenith_map;
	    zenith = min_zenith + izenith * dzenith;
	  }

	Gravity->zenith_map[izenith] = zenith;	  
	  zenith = zenith * M_PI / 180.;
	for (ielev_surf = 0; ielev_surf < nelev_surf; ielev_surf++){ // go over all "elevations" of satellite normal vector. See the definition of this "elevation" in comments a bit below
	  	  Gravity->earth_pressure_map[iradius][izenith][ielev_surf] = malloc(nazim_surf  * sizeof(double*) );
	  if ( Gravity->earth_pressure_map[iradius][izenith][ielev_surf] == NULL){ 
	    printf("***! Could not allow memory to Gravity->earth_pressure_map[iradius][izenith][ielev_surf]. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
	  }


	  if (ielev_surf == nelev_surf - 1){
	    delev_surf = max_elev_surf - (min_elev_surf + (nelev_surf - 2)*Gravity->delev_surf_map);//not used 
	    elev_surf = max_elev_surf;
	  }
	  else{
	    delev_surf =  Gravity->delev_surf_map;
	    elev_surf = min_elev_surf + ielev_surf * delev_surf             ;
	  }
	Gravity->elev_surf_map[ielev_surf] = elev_surf;	  
	  
	  elev_surf = elev_surf * M_PI / 180.;
	for (iazim_surf = 0; iazim_surf < nazim_surf; iazim_surf++){ // go over all "azimuth" of satellite normal vector. See the definition of this "azimuth" in comments a bit below
	  	  Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf] = malloc(3  * sizeof(double) );
	  if ( Gravity->earth_pressure_map[iradius][izenith][ielev_surf] == NULL){ 
	    printf("***! Could not allow memory to Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf]. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
	  }

	  if (iazim_surf == nazim_surf - 1){
	    dazim_surf = max_azim_surf - (min_azim_surf + (nazim_surf - 2)*Gravity->dazim_surf_map);//not used 
	    azim_surf = max_azim_surf;
	  }
	  else{
	    dazim_surf =  Gravity->dazim_surf_map;
	    azim_surf = min_azim_surf + iazim_surf * dazim_surf             ;
	  }
	Gravity->azim_surf_map[iazim_surf] = azim_surf;	  
	  
	  azim_surf = azim_surf * M_PI / 180.;
	  // the "elevation" of the satellite normal vector is the angle between nvec and CO (direction satellite to center of the Earth, ie -z direction where z has been previously defined)
 
	  	fread(&Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][0],
	      sizeof(Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][0]),1,file_earth_pressure_map);
	fread(&Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][1],
	      sizeof(Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][1]),1,file_earth_pressure_map);
	fread(&Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][2],
	      sizeof(Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][2]),1,file_earth_pressure_map);

	//	printf("earth[%f][%f][%f][%f]: (%.2e, %.2e, %.2e)\n", Gravity->radius_map[iradius] - Gravity->radius, Gravity->zenith_map[izenith], Gravity->elev_surf_map[ielev_surf], Gravity->azim_surf_map[iazim_surf], Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][0], Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][1], Gravity->earth_pressure_map[iradius][izenith][ielev_surf][iazim_surf][2]);
	
	} // go over all azimuth of satellite normal vector

	} // go over all elevation of satellite normal vector


	
	}// go over zeniths
	
  } // go over radii

    fclose(file_earth_pressure_map);

  return 0;
}


int build_gravity_map(GRAVITY_T  *Gravity, int degree,  int iProc){

  double dlat = Gravity->dlat_map;// size of latitude bins for the map (in degrees)
  double dlon = Gravity->dlon_map;// size of longitude bins for the map (in degrees)
  double dradius;// size of radius bins for the map (in km);

  int  iradius, ilon, ilat;
  double radius, lon,lat;
  
  double min_lat = Gravity->min_lat_map;// min latitude of all sc
  double max_lat = Gravity->max_lat_map;// max latitude of all sc
  double max_radius = Gravity->max_radius_map; // max radius of all sc (km)
  double min_radius = Gravity->min_radius_map; // min radius of all sc (km)
  // ex: if dlon = 1 deg then draf = 100 km.

  double dlonrad;
  double dlatrad;
  double min_latrad = min_lat * M_PI/ 180;
  double max_latrad = max_lat * M_PI/ 180;
 
  double slat;
  double radius2;
  double radius3;
  double rad_over_r;
  double Plm;
  double Plm_plus1;
  double term;
  double term2;
  double coeff;
  int l;
  int m;


 
  double dUdr     = 0.0;
  double dUdlat   = 0.0;
  double dUdlong  = 0.0;

  
  int nlon = Gravity->nlon_map;// nb lon bins
  int nlat = Gravity->nlat_map;// nb lat bins
  int nradius = Gravity->nradius_map;

  printf("3D Gravity map: dradius: %.1f km (%.1f to %.1f km, %d), dlat: %.1f deg (%d), dlon: %.1f deg (%d)\n\n",  Gravity->dradius_map, min_radius-Gravity->radius, max_radius-Gravity->radius, nradius, dlat, nlat, dlon, nlon);
  
  Gravity->file_gravity_map = NULL;
  Gravity->file_gravity_map = fopen(Gravity->filename_gravity_map, "wb");

  for (iradius = 0; iradius < nradius; iradius++){ // go over all radii
      	if (iradius == nradius - 1){
	  dradius = max_radius - (min_radius + (nradius - 2)*Gravity->dradius_map);//not used 
	  radius = max_radius;
	}
	else{
	  dradius =  Gravity->dradius_map;
	  radius = min_radius + iradius * dradius;
	}

    
      if (iProc == 0){ // print progress
	printf("\033[A\33[2K\rBuilding the 3D gravity map... %.0f%%\n",  iradius*100.  / ( nradius-1 ) );
      }



			printf("\n");
    for (ilat = 0; ilat < nlat; ilat++){ // go over all latitudes
      	if (ilat == nlat - 1){
	  dlatrad = max_latrad - (min_latrad + (nlat - 2)*dlatrad);//not used 
	  lat = max_latrad;
	}
	else{
	  dlatrad = Gravity->dlat_map * M_PI/ 180;
	  lat = min_latrad + ilat * dlatrad; // GEOCENTRIC latitude
	}



      
      if (iProc == 0){ // print progress
	//	      	printf("\033[A\33[2K\rBuilding the 3D gravity map... %.0f%%\n",  ilat*100.  / ( nlat-1 ) );
		      	printf("\033[A\33[2K\r... %.0f%%\n",  ilat*100.  / ( nlat-1 ) );
      }
      

      for (ilon = 0; ilon < nlon; ilon++){ // go over all longitude
	if (ilon == nlon - 1){
	  dlonrad = 2*M_PI - (nlon - 2)*dlonrad;//not used 
	  lon = 2*M_PI;
	}
	else{
	  dlonrad = Gravity->dlon_map * M_PI/ 180;
	lon = 0 + ilon * dlonrad;
	}

	//	printf("%f (%d), %f (%d), %f (%d)\n", radius, iradius, lat*180/M_PI, ilat, lon*180/M_PI, ilon);

	
	// Declarations

	dUdr     = 0.0;
	dUdlat   = 0.0;
	dUdlong  = 0.0;
    
	// Some intermediate
	rad_over_r  = Gravity->radius / radius;
	slat     = sin( lat );

	radius2       = pow( radius, 2);
	radius3       = pow( radius, 3);

	//   Compute Partial of potential function wrt range, lat, long
	for (l = 2; l <= degree; l++) {

	  for ( m = 0; m <= l; m++) { // !!!!!!!!!!!!! replace 0 with l

	    Plm         = pow(-1,m) * gsl_sf_legendre_Plm( l, m, slat ); // in ~/gsl-1.16/specfunc/legendre_poly.c // pow(-1,m) has been added to agree with Vallado's theory (the C library uses a pow(-1,m) that the theory does not))

	    if ((m+1) > l) {
            
	      Plm_plus1 = 0.0;
            
	    } else {
                
	      Plm_plus1   = pow(-1,m+1) * gsl_sf_legendre_Plm( l, (m+1), slat ); // !!!!!!!!!!  pow(-1,m+1) has been added to agree with Vallado's theory (the C library uses a pow(-1,m+1) that the theory does not))
            
	    }


	    coeff       = pow(rad_over_r, l);

	    term        = (Gravity->Clm[l][m] * cos( m*lon ) + Gravity->Slm[l][m] * sin( m*lon ));
	    term2       = (-Gravity->Clm[l][m] * sin( m*lon ) + Gravity->Slm[l][m] * cos( m*lon )); // Modif by CBV 07-19-2015 from Vallado3 p. 548
            
	    dUdr    += coeff * (l + 1.0) * Plm * term;
	    dUdlat  += coeff * ( Plm_plus1  - m * tan( lat) * Plm ) * term;
	    dUdlong += coeff * m * Plm * term2;
            
	  }
	}


	dUdr    = -dUdr     * Gravity->mu / radius2;
	dUdlat  = dUdlat    * Gravity->mu / radius;
	dUdlong = dUdlong   * Gravity->mu / radius;


	// on one line: all longitudes for a given latitude 
	fwrite(&(dUdr), sizeof(dUdr), 1, Gravity->file_gravity_map);
	fwrite(&(dUdlat), sizeof(dUdlat), 1, Gravity->file_gravity_map);
	fwrite(&(dUdlong), sizeof(dUdlong), 1, Gravity->file_gravity_map); 
      } // go over all longitudes
      //fprintf(Gravity->file_gravity_map, "\n"); // move to next line for next latitude // txt
    } // go over all latitudes
    //    fprintf(Gravity->file_gravity_map, "\n"); // move to next line for next radius // txt
  }
  fclose(Gravity->file_gravity_map);
	
  exitf();

  return 0;
}
 


int read_gravity_map(GRAVITY_T  *Gravity, int degree,  int iProc){


  double dlat = Gravity->dlat_map;// size of latitude bins for the map (in degrees)
  double dlon = Gravity->dlon_map;// size of longitude bins for the map (in degrees)
  double dradius = Gravity->dradius_map;// size of radius bins for the map (in km);

  int  iradius, ilon, ilat;
  double radius, lon,lat;
  
  double min_lat = Gravity->min_lat_map;// min latitude of all sc
  double max_lat = Gravity->max_lat_map;// max latitude of all sc
  double max_radius = Gravity->max_radius_map; // max radius of all sc (km)
  double min_radius = Gravity->min_radius_map; // min radius of all sc (km)
  // ex: if dlon = 1 deg then draf = 100 km.

  double dlonrad;
  double dlatrad;
  double min_latrad = min_lat * M_PI/ 180;
  double max_latrad = max_lat * M_PI/ 180;
 
  double slat;
  double radius2;
  double radius3;
  double rad_over_r;
  double Plm;
  double Plm_plus1;
  double term;
  double term2;
  double coeff;
  int l;
  int m;


 
  double dUdr     = 0.0;
  double dUdlat   = 0.0;
  double dUdlong  = 0.0;


  FILE *file_gravity_map = NULL;

  int nlon = Gravity->nlon_map;// nb lon bins
  int nlat = Gravity->nlat_map;// nb lat bins
  int nradius = Gravity->nradius_map;

  
  file_gravity_map = fopen(Gravity->filename_gravity_map, "r");

  if (iProc == 0){
  if (file_gravity_map == NULL){
    printf("!***!\nThe 3D gravity map file:\n%s\ncould not be found. The program will stop.\n!***!\n", Gravity->filename_gravity_map); MPI_Finalize(); exit(0);
  }
  }
  printf("3D map name: <%s>\n", Gravity->filename_gravity_map);  
  printf("3D Gravity map: dradius: %.1f km (%.1f to %.1f km, %d), dlat: %.1f deg (%d), dlon: %.1f deg (%d)\n\n",  Gravity->dradius_map, min_radius-Gravity->radius, max_radius-Gravity->radius, nradius, dlat, nlat, dlon, nlon);  

  Gravity->gravity_map = malloc(nradius  * sizeof(double ***) );
  Gravity->radius_map = malloc(nradius  * sizeof(double) );
    Gravity->lat_map = malloc(nlat  * sizeof(double) );
      Gravity->lon_map = malloc(nlon  * sizeof(double) );
     
  if ( Gravity->gravity_map == NULL){
    printf("***! Could not allow memory to Gravity->gravity_map. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  if ( Gravity->radius_map == NULL){
    printf("***! Could not allow memory to Gravity->radius_map. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  if ( Gravity->lat_map == NULL){
    printf("***! Could not allow memory to Gravity->lat_map. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  if ( Gravity->lon_map == NULL){
    printf("***! Could not allow memory to Gravity->lon_map. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }

  
  for (iradius = 0; iradius < nradius; iradius++){ // go over all radii
      if (iProc == 0){ // print progress
					              printf("\033[A\33[2K\rReading the 3D gravity map... %.0f%%\n",  iradius*100.  / ( nradius-1 ) );
      }
      	if (iradius == nradius - 1){
	  dradius = max_radius - (min_radius + (nradius - 2)*Gravity->dradius_map);//not used 
	  radius = max_radius;
	}
	else{
	  dradius =  Gravity->dradius_map;
	  radius = min_radius + iradius * dradius;
	}

	Gravity->radius_map[iradius] = radius;
    Gravity->gravity_map[iradius] = malloc(nlat  * sizeof(double **) );
    if ( Gravity->gravity_map[iradius] == NULL){
      printf("***! Could not allow memory to Gravity->gravity_map[iradius]. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }
    printf("\n");
    for (ilat = 0; ilat < nlat; ilat++){ // go over all latitudes
      	if (ilat == nlat - 1){
	  dlatrad = max_latrad - (min_latrad + (nlat - 2)*dlatrad);//not used 
	  lat = max_latrad;
	}
	else{
	  dlatrad = Gravity->dlat_map * M_PI/ 180;
	  lat = min_latrad + ilat * dlatrad; // GEOCENTRIC latitude
	}
	Gravity->lat_map[ilat] = lat * 180/M_PI;
	//	printf(, )
      if (iProc == 0){ // print progress
		      	printf("\033[A\33[2K\r... %.0f%%\n",  ilat*100.  / ( nlat-1 ) );
      }
      
      Gravity->gravity_map[iradius][ilat] = malloc(nlon  * sizeof(double*) );
      
      if ( Gravity->gravity_map[iradius][ilat] == NULL){ 
	printf("***! Could not allow memory to Gravity->gravity_map[iradius][ilat]. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
      }
      for (ilon = 0; ilon < nlon; ilon++){ // go over all longitude
	if (ilon == nlon - 1){
	  dlonrad = 2*M_PI - (nlon - 2)*dlonrad;//not used 
	  lon = 2*M_PI;
	}
	else{
	  dlonrad = Gravity->dlon_map * M_PI/ 180;
	lon = 0 + ilon * dlonrad;
	}
	Gravity->lon_map[ilon] = lon*180/M_PI;
	
	//	printf("%f (%d), %f (%d), %f (%d)\n", radius, iradius, lat*180/M_PI, ilat,  Gravity->lon_map[ilon]*180/M_PI, ilon);
	Gravity->gravity_map[iradius][ilat][ilon] = malloc(3  * sizeof(double) );
	if ( Gravity->gravity_map[iradius][ilat][ilon] == NULL){ 
	  printf("***! Could not allow memory to Gravity->gravity_map[iradius][ilat][ilon]. \
The program will stop. !***\n"); MPI_Finalize(); exit(0);
	}
	fread(&Gravity->gravity_map[iradius][ilat][ilon][0],
	      sizeof(Gravity->gravity_map[iradius][ilat][ilon][0]),1,file_gravity_map);
	fread(&Gravity->gravity_map[iradius][ilat][ilon][1],
	      sizeof(Gravity->gravity_map[iradius][ilat][ilon][1]),1,file_gravity_map);
	fread(&Gravity->gravity_map[iradius][ilat][ilon][2],
	      sizeof(Gravity->gravity_map[iradius][ilat][ilon][2]),1,file_gravity_map);

      } // go over all longitudes

    } // go over all latitudes

  }
	


  return 0;
}




double compute_earth_albedo(){ // see http://ocw.upm.es/ingenieria-aeroespacial/modeling-the-space-environment/contenidos/material-declase/mse07_radiationpressure.pdf or Knocke, P.C. et al, “Earth Radiation Pressure Effects on Satellites.” Proceedings of the AIAA/AAS Astrodynamics Specialist Conference, Washington DC, 1988, pp. 577-586.
  double albedo;
  albedo = 0.34; // simplification (complete equation at the reference above)
  return albedo;
}

double compute_earth_emissivity(){ // see http://ocw.upm.es/ingenieria-aeroespacial/modeling-the-space-environment/contenidos/material-declase/mse07_radiationpressure.pdf or Knocke, P.C. et al, “Earth Radiation Pressure Effects on Satellites.” Proceedings of the AIAA/AAS Astrodynamics Specialist Conference, Washington DC, 1988, pp. 577-586.
  double emiss;
  emiss = 0.68; // simplification (complete equation at the reference above)
  return emiss;
}

int polynomial_interpo(double *val, int n, double *x, double *y, double to_inter){
  //  printf("%d\n", n);
  float s=1,t=1;
  int i,j;
  *val = 0;
  /* if (printVar == 1){ */
  /*   for(i=0; i<n; i++){ */
  /*     for(j=i+1; j<n; j++){ */
  /* 	  if (x[j] == x[i]){ */
  /* 	    printf("same %d %d\n", i, j); */
  /* 	  } */
  /*     } */
  /*   } */
  /* } */
  for(i=0; i<n; i++){
    s=1;
    t=1;
    for(j=0; j<n; j++)
      {
	if(j!=i)
	  {
	    s=s*(to_inter-x[j]);
	    t=t*(x[i]-x[j]);
	    /* if (printVar == 1){ */
	    /*   printf("%d %d %f\n", i,j, x[i]-x[j]); */
	    /* } */
	  }
      }
    
    /* if (printVar == 1){ */
    /*   printf("t %f %d\n", t, i);   */
    /* } */
    
    *val=*val+((s/t)*y[i]);

  }
  return 0;

}

int polynomial_interpo_earth(double *val, int n, double *x, double **y, double to_inter){
  //  printf("%d\n", n);
  float s=1,t=1;
  int i,j, c;
  for (c = 0; c < 3; c++){
    val[c] = 0;
  }
  for(i=0; i<n; i++){
    s=1;
    t=1;
    for(j=0; j<n; j++)
      {
	if(j!=i)
	  {
	    s=s*(to_inter-x[j]);
	    t=t*(x[i]-x[j]);
	  }
      }

    for (c = 0; c < 3; c++){ 
    val[c]=val[c]+((s/t)*y[c][i]);
    }

  }
  return 0;

}



int compute_iradius_gravity_map(int iradius_arr[4], GRAVITY_T *Gravity, double rmag){

  int iradius0, iradius1, iradius2, iradius3;
  
  
    if (rmag < Gravity->min_radius_map){ // this happens if Gravity->min_radius_map was not set correctly, ie too large compared to the min radius of the satellite
      iradius0 = 0;
      iradius1 = 0;
      iradius2 = 0;
      iradius3 = 0;
      //      xradius = 0;
    }
    else if ((rmag < Gravity->min_radius_map+Gravity->dradius_map) && (rmag >= Gravity->min_radius_map)){
      iradius0 = 0;
      iradius1 = (int)(( rmag - Gravity->min_radius_map ) / Gravity->dradius_map); // should be 0
      iradius2 = iradius1 + 1;  // 1
      iradius3 = iradius2 + 1; // 2

    }
    else if (rmag >= Gravity->max_radius_map + Gravity->dradius_map){
      iradius1 = Gravity->nradius_map - 1;
      iradius2 = iradius1;
      iradius3 = iradius2;
      iradius0 = iradius1;
      //      xradius = 0;
    }
    else if ((rmag >= Gravity->max_radius_map) && (rmag < Gravity->max_radius_map + Gravity->dradius_map)){
      iradius0 = Gravity->nradius_map - 2;
      iradius1 = Gravity->nradius_map - 1;
      iradius2 = iradius1;
      iradius3 = iradius2;
      //      xradius = 0;
    }
    
    else if ((rmag >= Gravity->max_radius_map-Gravity->dradius_map) && (rmag < Gravity->max_radius_map)){
      iradius1 = (int)(( rmag - Gravity->min_radius_map ) / Gravity->dradius_map);
      iradius2 = iradius1 + 1; 
      iradius3 = iradius2;
      iradius0 = iradius1-1;
	}
    else{
      iradius1 = (int)(( rmag - Gravity->min_radius_map ) / Gravity->dradius_map);
      iradius2 = iradius1 + 1; 
      iradius3 = iradius2 + 1; 
      iradius0 = iradius1 - 1;
      //      xradius = (rmag - Gravity->radius_map[iradius1]) / Gravity->dradius_map;
    }

    iradius_arr[0] = iradius0; iradius_arr[1] = iradius1; iradius_arr[2] = iradius2; iradius_arr[3] = iradius3; 

    return 0;
    
    }

int compute_iradius_earth_pressure_map(double *xinter_radius, int iradius_arr[2], int *nradius_interpo, GRAVITY_T *Gravity, double rmag){
  int iradius0, iradius1;
  if (rmag < Gravity->min_radius_map){ // iradius1 won't be used
    iradius0 = 0;
    xinter_radius[0] = Gravity->radius_map[iradius0];
    *nradius_interpo = 1;
  }
  else if (rmag >= Gravity->max_radius_map){ // iradius1 won't be used
    iradius0 = Gravity->nradius_map - 1;
    xinter_radius[0] = Gravity->radius_map[iradius0];
    *nradius_interpo = 1;
  }
  else{
    iradius0 = (int)(( rmag - Gravity->min_radius_map ) / Gravity->dradius_map);
    iradius1 = iradius0 + 1;
    xinter_radius[0] = Gravity->radius_map[iradius0];   xinter_radius[1] = Gravity->radius_map[iradius1];
    *nradius_interpo = 2;
  }

  iradius_arr[0] = iradius0; iradius_arr[1] = iradius1; 

  return 0;    
}



int compute_ilat_gravity_map(int ilat_arr[4], GRAVITY_T *Gravity, double lat_gc){

  int ilat0, ilat1, ilat2, ilat3;
  
    // determine the bin for the lat
    if (lat_gc*180/M_PI < Gravity->min_lat_map){ // this happens if Gravity->min_lat_map was not set correctly, ie too large compared to the min latitude of the satellite
      ilat0 = 0;
      ilat1 = 0;
      ilat2 = 0;
      ilat3 = 0;
      //           xlat = 0;
    }
    else if ((lat_gc*180/M_PI < Gravity->min_lat_map + Gravity->dlat_map) && (lat_gc*180/M_PI >= Gravity->min_lat_map)){
      ilat0 = 0;
      ilat1 = (int)(( lat_gc*180/M_PI - Gravity->min_lat_map ) / Gravity->dlat_map); // should be 0
      ilat2 = ilat1 + 1; // 1
      ilat3 = ilat2 + 1; // 2
    }
    else if (lat_gc*180/M_PI >= Gravity->max_lat_map + Gravity->dlat_map){
      ilat1 = Gravity->nlat_map - 1;
      ilat2 = ilat1;
      ilat3 = ilat1;
      ilat0 = ilat1;
      //            xlat = 0;
    }
    else if ((lat_gc*180/M_PI >= Gravity->max_lat_map) && (Gravity->max_lat_map < Gravity->max_lat_map + Gravity->dlat_map)){// this happens only if the Gravity->max_lat_map was not set correctly, ie too small compared to the max latitude of the satellite
      ilat1 = Gravity->nlat_map - 1;
      ilat2 = ilat1;
      ilat3 = ilat1;
      ilat0 = ilat1-1;
      //            xlat = 0;
    }

    else if ( (lat_gc*180/M_PI >= Gravity->max_lat_map-Gravity->dlat_map) && (lat_gc*180/M_PI < Gravity->max_lat_map)){
      ilat1 = (int)(( lat_gc*180/M_PI - Gravity->min_lat_map ) / Gravity->dlat_map);
      ilat2 = ilat1 + 1;
      ilat3 = ilat2;
      ilat0 = ilat1 - 1;
      //      xlat = (lat_gc*180/M_PI - Gravity->lat_map[ilat1]) / Gravity->dlat_map;
    }

    else{
      ilat1 = (int)(( lat_gc*180/M_PI - Gravity->min_lat_map ) / Gravity->dlat_map);
      ilat2 = ilat1 + 1; 
      ilat3 = ilat2 + 1; 
      ilat0 = ilat1 - 1;
      // so we're safe doing: ilat2 = ilat1 + 1
      //           xlat = (lat_gc*180/M_PI - Gravity->lat_map[ilat1]) / Gravity->dlat_map;
    }


    ilat_arr[0] = ilat0; ilat_arr[1] = ilat1; ilat_arr[2] = ilat2; ilat_arr[3] = ilat3; 
    return 0;
    
    }




int compute_izenith_earth_pressure_map(double *xinter_zenith, int izenith_arr[2], int *nzenith_interpo, GRAVITY_T *Gravity, double zenith){

  int izenith0, izenith1;
  
    // determine the bin for the zenith

  if (zenith*180/M_PI < Gravity->min_zenith_map){ // izenith1 won't be used
      izenith0 = 0; 
      xinter_zenith[0] = Gravity->zenith_map[izenith0];
      *nzenith_interpo = 1;
    }

    else if (zenith*180/M_PI >= Gravity->max_zenith_map){ // izenith1 won't be used
      izenith0 = Gravity->nzenith_map - 1;;
      xinter_zenith[0] = Gravity->zenith_map[izenith0];
      *nzenith_interpo = 1;
    }

    else{
      izenith0 = (int)(( zenith*180/M_PI - Gravity->min_zenith_map ) / Gravity->dzenith_map);
      izenith1 = izenith0 + 1;
      xinter_zenith[0] = Gravity->zenith_map[izenith0]; xinter_zenith[1] = Gravity->zenith_map[izenith1];
      *nzenith_interpo = 2;
    }
    izenith_arr[0] = izenith0; izenith_arr[1] = izenith1; 
    return 0;
    
    }


int compute_ielev_surf_earth_pressure_map(double *xinter_elev_surf, int ielev_surf_arr[2], int *nelev_interpo, GRAVITY_T *Gravity, double elev_surf){
    int ielev_surf0, ielev_surf1;
    // determine the bin for the elev_surf
    if (elev_surf*180/M_PI < Gravity->min_elev_surf_map ){ // ielev_surf1 wont' eb sued outside this function anyway
      ielev_surf0 = 0;
      xinter_elev_surf[0] = Gravity->elev_surf_map[ielev_surf0];
      *nelev_interpo = 1;
    }
   else if (elev_surf*180/M_PI >= Gravity->max_elev_surf_map){ // ielev_surf1 wont' eb sued outside this function anyway
      ielev_surf0 = Gravity->nelev_surf_map-1;
      xinter_elev_surf[0] = Gravity->elev_surf_map[ielev_surf0];
      *nelev_interpo = 1;
   }
    else{
      ielev_surf0 = (int)(( elev_surf*180/M_PI - Gravity->min_elev_surf_map ) / Gravity->delev_surf_map);
      ielev_surf1 = ielev_surf0 + 1;
      xinter_elev_surf[0] = Gravity->elev_surf_map[ielev_surf0]; xinter_elev_surf[1] = Gravity->elev_surf_map[ielev_surf1];
      *nelev_interpo = 2;
    }

    ielev_surf_arr[0] = ielev_surf0; ielev_surf_arr[1] = ielev_surf1; 
    return 0;
    
    }






int compute_ilon_gravity_map(int ilon_arr[4], GRAVITY_T *Gravity, double long_gc_corr){

  int ilon0, ilon1, ilon2, ilon3;


    ilon1 = (int)(( long_gc_corr*180/M_PI ) / Gravity->dlon_map);
    // the longitude array of the 3d map varies frm 0 (1st bin) to 360 (last bin)
    // so it's impossible that ilon1 == nlon - 1 (since long_gc_corr can't be
    // exactly equal to 360). So we're safe to do: ilon2 = ilon1+1
        ilon2 = ilon1+1;
	if (long_gc_corr*180/M_PI >= 360 - Gravity->dlon_map){
	  ilon3 = 1;//  ilon2;
    }
    else{
    ilon3 = ilon2+1;
    }
    if (long_gc_corr*180/M_PI <= Gravity->dlon_map){ 
      ilon0 = Gravity->nlon_map-2; // example: ilon1 = 0 and ilon0 = 359 (Gravity->nlon_map-1 = 360 so neeed to have ilon0 = Gravity->nlon_map-2, not ilon0 = Gravity->nlon_map-1)
    }
    else{
    ilon0 = ilon1-1;
    }
  
    ilon_arr[0] = ilon0; ilon_arr[1] = ilon1; ilon_arr[2] = ilon2; ilon_arr[3] = ilon3; 
    return 0;
    
    }



int compute_iazim_surf_earth_pressure_map(double *xinter_azim_surf, int iazim_surf_arr[2], int *nazim_interpo, GRAVITY_T *Gravity, double azim_surf_corr){
  int iazim_surf0, iazim_surf1, iazim_surf2, iazim_surf3;
  iazim_surf0 = (int)(( azim_surf_corr*180/M_PI ) / Gravity->dazim_surf_map);
  iazim_surf1 = iazim_surf0+1; // the azim_surf_map includes 360 so if azim_surf_corr is for example 359.999 then azim_surf_map[iazim_surf0] is 359 and azim_surf_map[iazim_surf1] is 360.
  *nazim_interpo = 2;
  xinter_azim_surf[0] = Gravity->azim_surf_map[iazim_surf0]; xinter_azim_surf[1] = Gravity->azim_surf_map[iazim_surf1];
  iazim_surf_arr[0] = iazim_surf0; iazim_surf_arr[1] = iazim_surf1; 
  return 0;
}


int gravity_map_yinter_lon(double *yinter, int *order_interpo_map, double long_gc_corr, GRAVITY_T *Gravity, double y_lon0, double y_lon1, double y_lon2, double y_lon3){
    yinter[0] = y_lon0; yinter[1] = y_lon1; yinter[2] = y_lon2; yinter[3] = y_lon3;

    /* if (long_gc_corr*180/M_PI >= 360 - Gravity->dlon_map){ */
    /*   *order_interpo_map = 2; // otherwise the denominator in the interpolation formula is 0 */
    /* } */
    /* else{ */
    /*   *order_interpo_map = 4; */
    /* } */
    return 0;

}


int earth_pressure_map_yinter_azim_surf(double **yinter, int *order_interpo_map, double azim_surf_corr, GRAVITY_T *Gravity, double **y_az, int iazim_surf_arr[2]){
  int ii;
  for (ii = 0; ii < 3; ii++){
    yinter[ii][0] = y_az[iazim_surf_arr[0]][ii]; yinter[ii][1] = y_az[iazim_surf_arr[1]][ii]; 
  }
    return 0;

}

int gravity_map_yinter_lat(double *yinter, int *order_interpo_map, double lat_gc, GRAVITY_T *Gravity, double y_lat0, double y_lat1, double y_lat2, double y_lat3){
      if ((lat_gc*180/M_PI < Gravity->min_lat_map) || (lat_gc*180/M_PI >= Gravity->max_lat_map + Gravity->dlat_map)){ // all ilat the same
	yinter[0] = y_lat0;
	*order_interpo_map = 1;
      }
      else if ((lat_gc*180/M_PI < Gravity->min_lat_map + Gravity->dlat_map) && (lat_gc*180/M_PI >= Gravity->min_lat_map)){ // ilat1 = ilat0
	yinter[0] = y_lat0; yinter[1] = y_lat2; yinter[2] = y_lat3;
	*order_interpo_map = 3;
      }
      else if ((lat_gc*180/M_PI >= Gravity->max_lat_map) && (Gravity->max_lat_map < Gravity->max_lat_map + Gravity->dlat_map)){ //ilat3 = ilat2 = ilat1
	yinter[0] = y_lat0; yinter[1] = y_lat1;
	*order_interpo_map = 2;
      }
      else if ( (lat_gc*180/M_PI >= Gravity->max_lat_map-Gravity->dlat_map) && (lat_gc*180/M_PI < Gravity->max_lat_map)){ // ilat3 =ilat2
	yinter[0] = y_lat0; yinter[1] = y_lat1; yinter[2] = y_lat2;
	*order_interpo_map = 3;
    }
      else{
	yinter[0] = y_lat0; yinter[1] = y_lat1; yinter[2] = y_lat2; yinter[3] = y_lat3;
    	*order_interpo_map = 4;
      }

return 0;
}



int earth_pressure_map_yinter_elev_surf(double **yinter, int *order_interpo_map, double elev_surf, GRAVITY_T *Gravity, double **y_el){
    int ii;
  for (ii = 0; ii < 3; ii++){

    if (elev_surf*180/M_PI < Gravity->min_elev_surf_map){
	yinter[ii][0] = y_el[0][ii];
	*order_interpo_map = 1;
      }
    else if (elev_surf*180/M_PI >= Gravity->max_elev_surf_map){
	yinter[ii][0] = y_el[0][ii]; 
	*order_interpo_map = 1;
    }
      else{
	yinter[ii][0] = y_el[0][ii]; yinter[ii][1] = y_el[1][ii];
    	*order_interpo_map = 2;
      }
  }
return 0;
}


int earth_pressure_map_yinter_zenith(double **yinter, int *order_interpo_map, double zenith, GRAVITY_T *Gravity, double **y_zen){
      int ii;
  for (ii = 0; ii < 3; ii++){


    if (zenith*180/M_PI < Gravity->min_zenith_map){
	yinter[ii][0] = y_zen[0][ii]; 
	*order_interpo_map = 1;
      }
    else if (zenith*180/M_PI >= Gravity->max_zenith_map){
	yinter[ii][0] = y_zen[0][ii]; 
	*order_interpo_map = 1;
    }
      else{
	yinter[ii][0] = y_zen[0][ii]; yinter[ii][1] = y_zen[1][ii];
    	*order_interpo_map = 2;
      }
  }
return 0;
}


int gravity_map_yinter_radius_earth(double **yinter, int *order_interpo_map, double rmag, GRAVITY_T *Gravity, double **y_rad){
      int ii;
  for (ii = 0; ii < 3; ii++){
    if (rmag < Gravity->min_radius_map){
      yinter[ii][0] = y_rad[0][ii];
      *order_interpo_map = 1;
    }
    else if (rmag >= Gravity->max_radius_map){
      yinter[ii][0] = y_rad[0][ii]; 
      *order_interpo_map = 1;
    }
    else{
      yinter[ii][0] = y_rad[0][ii]; yinter[ii][1] = y_rad[1][ii]; 
      *order_interpo_map = 2;
    }
  }
  return 0;
}

int gravity_map_yinter_radius(double *yinter, int *order_interpo_map, double rmag, GRAVITY_T *Gravity, double y_radius0, double y_radius1, double y_radius2, double y_radius3){

      if ((rmag < Gravity->min_radius_map) || (rmag >= Gravity->max_radius_map + Gravity->dradius_map)){ // all iradius the same
      yinter[0] = y_radius0;
      *order_interpo_map = 1;

    }
    else if ((rmag < Gravity->min_radius_map+Gravity->dradius_map) && (rmag >= Gravity->min_radius_map)){ // iradius0 = iradius1
      yinter[0] = y_radius0; yinter[1] = y_radius2; yinter[2] = y_radius3;
      *order_interpo_map = 3;

    }
    else if ((rmag >= Gravity->max_radius_map) && (rmag < Gravity->max_radius_map + Gravity->dradius_map)){ // iradius2 and 3 = iradius 1
      yinter[0] = y_radius0; yinter[1] = y_radius1;
      *order_interpo_map = 2;

    }
    else if ((rmag >= Gravity->max_radius_map-Gravity->dradius_map) && (rmag < Gravity->max_radius_map)){ // iradius3 = iradius2
      yinter[0] = y_radius0; yinter[1] = y_radius1; yinter[2] = y_radius2;
      *order_interpo_map = 3;

    }
    else{
	    yinter[0] = y_radius0; yinter[1] = y_radius1; yinter[2] = y_radius2; yinter[3] = y_radius3;
	    *order_interpo_map = 4;
    }

  return 0;
}


int earth_pressure_map_xinter_elev_azim_surf(double *xinter_radius, double *xinter_zenith, double *xinter_elev_surf, double *xinter_azim_surf, double rmag, double zenith, double elev_surf, double azim_surf_corr, GRAVITY_T *Gravity, int *iradius_arr, int *izenith_arr, int *ielev_surf_arr, int *iazim_surf_arr)
	  {
  
  int iradius1, iradius2, iradius3, iradius0;
  int izenith0, izenith1, izenith2, izenith3;
  int ielev_surf0, ielev_surf1, ielev_surf2, ielev_surf3;
  int iazim_surf0, iazim_surf1, iazim_surf2, iazim_surf3;
  iradius0 = iradius_arr[0]; iradius1 = iradius_arr[1]; iradius2 = iradius_arr[2]; iradius3 = iradius_arr[3];
  izenith0 = izenith_arr[0]; izenith1 = izenith_arr[1]; izenith2 = izenith_arr[2]; izenith3 = izenith_arr[3];
  ielev_surf0 = ielev_surf_arr[0]; ielev_surf1 = ielev_surf_arr[1]; ielev_surf2 = ielev_surf_arr[2]; ielev_surf3 = ielev_surf_arr[3];
  iazim_surf0 = iazim_surf_arr[0]; iazim_surf1 = iazim_surf_arr[1]; iazim_surf2 = iazim_surf_arr[2]; iazim_surf3 = iazim_surf_arr[3];
  
    // xinter_azim_surf
      if (azim_surf_corr*180/M_PI <= Gravity->dazim_surf_map){ 
	xinter_azim_surf[0] = Gravity->azim_surf_map[Gravity->nazim_surf_map-2] - Gravity->azim_surf_map[Gravity->nazim_surf_map-1]; 
	xinter_azim_surf[1] = Gravity->azim_surf_map[iazim_surf1]; xinter_azim_surf[2] = Gravity->azim_surf_map[iazim_surf2]; xinter_azim_surf[3] = Gravity->azim_surf_map[iazim_surf3];
      }
      else if (azim_surf_corr*180/M_PI >= 360 - Gravity->dazim_surf_map){
	xinter_azim_surf[3] = Gravity->azim_surf_map[Gravity->nazim_surf_map-1] + (Gravity->azim_surf_map[1] - Gravity->azim_surf_map[0]);
	xinter_azim_surf[0] = Gravity->azim_surf_map[iazim_surf0]; xinter_azim_surf[1] = Gravity->azim_surf_map[iazim_surf1]; xinter_azim_surf[2] = Gravity->azim_surf_map[iazim_surf2];
      }
      else{
	xinter_azim_surf[0] = Gravity->azim_surf_map[iazim_surf0]; xinter_azim_surf[1] = Gravity->azim_surf_map[iazim_surf1]; xinter_azim_surf[2] = Gravity->azim_surf_map[iazim_surf2]; xinter_azim_surf[3] = Gravity->azim_surf_map[iazim_surf3];

      }

      // xinter_zenith
      if ((zenith*180/M_PI < Gravity->min_zenith_map + Gravity->dzenith_map) && (zenith*180/M_PI >= Gravity->min_zenith_map)){ // izenith1 = izenith0
	xinter_zenith[0] = Gravity->zenith_map[izenith0]; xinter_zenith[1] = Gravity->zenith_map[izenith2]; xinter_zenith[2] =  Gravity->zenith_map[izenith3];
      }
      else if ( (zenith*180/M_PI >= Gravity->max_zenith_map-Gravity->dzenith_map) && (zenith*180/M_PI < Gravity->max_zenith_map)){ // izenith3 =izenith2
	xinter_zenith[0] = Gravity->zenith_map[izenith0]; xinter_zenith[1] = Gravity->zenith_map[izenith1]; xinter_zenith[2] = Gravity->zenith_map[izenith2];
    }
      else{
	xinter_zenith[0] = Gravity->zenith_map[izenith0]; xinter_zenith[1] = Gravity->zenith_map[izenith1]; xinter_zenith[2] = Gravity->zenith_map[izenith2]; xinter_zenith[3] = Gravity->zenith_map[izenith3];
      }


      // xinter_elev_surf
      if ((elev_surf*180/M_PI < Gravity->min_elev_surf_map + Gravity->delev_surf_map) && (elev_surf*180/M_PI >= Gravity->min_elev_surf_map)){ // ielev_surf1 = ielev_surf0
	xinter_elev_surf[0] = Gravity->elev_surf_map[ielev_surf0]; xinter_elev_surf[1] = Gravity->elev_surf_map[ielev_surf2]; xinter_elev_surf[2] =  Gravity->elev_surf_map[ielev_surf3];
      }
      else if ( (elev_surf*180/M_PI >= Gravity->max_elev_surf_map-Gravity->delev_surf_map) && (elev_surf*180/M_PI <= Gravity->max_elev_surf_map)){ // ielev_surf3 =ielev_surf2
	xinter_elev_surf[0] = Gravity->elev_surf_map[ielev_surf0]; xinter_elev_surf[1] = Gravity->elev_surf_map[ielev_surf1]; xinter_elev_surf[2] = Gravity->elev_surf_map[ielev_surf2];
    }
      else{
	xinter_elev_surf[0] = Gravity->elev_surf_map[ielev_surf0]; xinter_elev_surf[1] = Gravity->elev_surf_map[ielev_surf1]; xinter_elev_surf[2] = Gravity->elev_surf_map[ielev_surf2]; xinter_elev_surf[3] = Gravity->elev_surf_map[ielev_surf3];
      }

      
      // xinter_radius
      if ((rmag < Gravity->min_radius_map) || (rmag >= Gravity->max_radius_map + Gravity->dradius_map)){ // all iradius the same
      xinter_radius[0] = Gravity->radius_map[iradius0];
    }
    else if ((rmag < Gravity->min_radius_map+Gravity->dradius_map) && (rmag >= Gravity->min_radius_map)){ // iradius0 = iradius1
      xinter_radius[0] = Gravity->radius_map[iradius0]; xinter_radius[1] =  Gravity->radius_map[iradius2]; xinter_radius[2] = Gravity->radius_map[iradius3];
    }
    else if ((rmag >= Gravity->max_radius_map) && (rmag < Gravity->max_radius_map + Gravity->dradius_map)){ // iradius2 and 3 = iradius 1
      xinter_radius[0] = Gravity->radius_map[iradius0]; xinter_radius[1] =  Gravity->radius_map[iradius1];
    }
    else if ((rmag >= Gravity->max_radius_map-Gravity->dradius_map) && (rmag < Gravity->max_radius_map)){ // iradius3 = iradius2
      xinter_radius[0] = Gravity->radius_map[iradius0]; xinter_radius[1] =  Gravity->radius_map[iradius1]; xinter_radius[2] = Gravity->radius_map[iradius2];
    }
    else{
          xinter_radius[0] = Gravity->radius_map[iradius0];   xinter_radius[1] = Gravity->radius_map[iradius1]; xinter_radius[2] = Gravity->radius_map[iradius2]; xinter_radius[3] = Gravity->radius_map[iradius3];
    }

  return 0;
}

int gravity_map_xinter(double *xinter_lon, double *xinter_lat, double *xinter_radius,  double long_gc_corr, double lat_gc, double rmag, GRAVITY_T *Gravity,
		       int *ilon_arr, int *ilat_arr, int *iradius_arr){
  
  int iradius1, ilat1, ilon1, iradius2, ilat2, ilon2, ilon3, ilat3, iradius3, ilon0, ilat0, iradius0;
  ilon0 = ilon_arr[0]; ilon1 = ilon_arr[1]; ilon2 = ilon_arr[2]; ilon3 = ilon_arr[3];
  ilat0 = ilat_arr[0]; ilat1 = ilat_arr[1]; ilat2 = ilat_arr[2]; ilat3 = ilat_arr[3];
  iradius0 = iradius_arr[0]; iradius1 = iradius_arr[1]; iradius2 = iradius_arr[2]; iradius3 = iradius_arr[3];
  
    // xinter_lon
      if (long_gc_corr*180/M_PI <= Gravity->dlon_map){ // ilon0 = Gravity->nlon_map-2 but ilon1 = 0 so avoid a difference of about 360 degrees between both lon by setting xinter_lon[0] to, for example, -1 (instead of 359)
	xinter_lon[0] = Gravity->lon_map[Gravity->nlon_map-2] - Gravity->lon_map[Gravity->nlon_map-1]; 
	xinter_lon[1] = Gravity->lon_map[ilon1]; xinter_lon[2] = Gravity->lon_map[ilon2]; xinter_lon[3] = Gravity->lon_map[ilon3];
      }
      else if (long_gc_corr*180/M_PI >= 360 - Gravity->dlon_map){
	xinter_lon[3] = Gravity->lon_map[Gravity->nlon_map-1] + (Gravity->lon_map[1] - Gravity->lon_map[0]);
	xinter_lon[0] = Gravity->lon_map[ilon0]; xinter_lon[1] = Gravity->lon_map[ilon1]; xinter_lon[2] = Gravity->lon_map[ilon2];
      }
      else{
	xinter_lon[0] = Gravity->lon_map[ilon0]; xinter_lon[1] = Gravity->lon_map[ilon1]; xinter_lon[2] = Gravity->lon_map[ilon2]; xinter_lon[3] = Gravity->lon_map[ilon3];

      }

      // xinter_lat
      if ((lat_gc*180/M_PI < Gravity->min_lat_map) || (lat_gc*180/M_PI >= Gravity->max_lat_map + Gravity->dlat_map)){ // all ilat the same
	xinter_lat[0] = Gravity->lat_map[ilat0];
      }
      else if ((lat_gc*180/M_PI < Gravity->min_lat_map + Gravity->dlat_map) && (lat_gc*180/M_PI >= Gravity->min_lat_map)){ // ilat1 = ilat0
	xinter_lat[0] = Gravity->lat_map[ilat0]; xinter_lat[1] = Gravity->lat_map[ilat2]; xinter_lat[2] =  Gravity->lat_map[ilat3];
      }
      else if ((lat_gc*180/M_PI >= Gravity->max_lat_map) && (Gravity->max_lat_map < Gravity->max_lat_map + Gravity->dlat_map)){ //ilat3 = ilat2 = ilat1
	xinter_lat[0] = Gravity->lat_map[ilat0]; xinter_lat[1] = Gravity->lat_map[ilat1];
      }
      else if ( (lat_gc*180/M_PI >= Gravity->max_lat_map-Gravity->dlat_map) && (lat_gc*180/M_PI < Gravity->max_lat_map)){ // ilat3 =ilat2
	xinter_lat[0] = Gravity->lat_map[ilat0]; xinter_lat[1] = Gravity->lat_map[ilat1]; xinter_lat[2] = Gravity->lat_map[ilat2];
    }
      else{
	xinter_lat[0] = Gravity->lat_map[ilat0]; xinter_lat[1] = Gravity->lat_map[ilat1]; xinter_lat[2] = Gravity->lat_map[ilat2]; xinter_lat[3] = Gravity->lat_map[ilat3];
      }

      // xinter_radius
      if ((rmag < Gravity->min_radius_map) || (rmag >= Gravity->max_radius_map + Gravity->dradius_map)){ // all iradius the same
      xinter_radius[0] = Gravity->radius_map[iradius0];
    }
    else if ((rmag < Gravity->min_radius_map+Gravity->dradius_map) && (rmag >= Gravity->min_radius_map)){ // iradius0 = iradius1
      xinter_radius[0] = Gravity->radius_map[iradius0]; xinter_radius[1] =  Gravity->radius_map[iradius2]; xinter_radius[2] = Gravity->radius_map[iradius3];
    }
    else if ((rmag >= Gravity->max_radius_map) && (rmag < Gravity->max_radius_map + Gravity->dradius_map)){ // iradius2 and 3 = iradius 1
      xinter_radius[0] = Gravity->radius_map[iradius0]; xinter_radius[1] =  Gravity->radius_map[iradius1];
    }
    else if ((rmag >= Gravity->max_radius_map-Gravity->dradius_map) && (rmag < Gravity->max_radius_map)){ // iradius3 = iradius2
      xinter_radius[0] = Gravity->radius_map[iradius0]; xinter_radius[1] =  Gravity->radius_map[iradius1]; xinter_radius[2] = Gravity->radius_map[iradius2];
    }
    else{
          xinter_radius[0] = Gravity->radius_map[iradius0];   xinter_radius[1] = Gravity->radius_map[iradius1]; xinter_radius[2] = Gravity->radius_map[iradius2]; xinter_radius[3] = Gravity->radius_map[iradius3];
    }

  return 0;
}



int compute_thrust(double athrust_i2cg_INRTL[3], // if section #THRUST exists in the main input file, this function computes the thrust acceleration in the ECI reference frame;
		   double r_i2cg_INRTL[3],
		   double v_i2cg_INRTL[3],
		   OPTIONS_T       *OPTIONS
		   ) {
  double T_lvlh_2_inrtl[3][3];
  double T_inrtl_2_lvlh[3][3];
    compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, r_i2cg_INRTL, v_i2cg_INRTL);
  m_trans(T_lvlh_2_inrtl, T_inrtl_2_lvlh);

  m_x_v(athrust_i2cg_INRTL, T_lvlh_2_inrtl, OPTIONS->thrust_accel_lvlh);

  //  v_print(athrust_i2cg_INRTL, "athrust_i2cg_INRTL");

  
  return 0;
}
