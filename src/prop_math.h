 #ifndef __PROP_MATH_H__
#define __PROP_MATH_H__
#include <math.h>
#include "mpi_fake/mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "propagator.h"
#define SMALL_NUM 1.0e-12
#define DEG2RAD 0.01745329251994
#define RAD2DEG 57.29577951308233
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))

#define  CONST_MAT          ( ConstSpiceDouble   (*) [3] )


// Matrix Printing Functions
int m_print(double m_to_print[3][3],
        char name[256]);

int m_print6(double m_to_print[6][6],
         char name[256]);

int m_print6_temp(double m_to_print[6][6],
         char name[256]);

int m_print9_temp(double m_to_print[9][9],
          char name[256]);

int m_print8_temp(double m_to_print[8][8],
          char name[256]);

// Vector Printing Functions
int v_norm_print(double v_to_print_norm[3],
         char name[256]);

int v_print( double v_to_print[3],
         char name[256]);

int v_print6( double v_to_print[6],
         char name[256]);
int v_print7( double v_to_print[7],
          char name[256]);

// Vector and Matrix math
// Compute the vector magnitude
int v_mag( double *v_mag,
           double v_in[3]);

// Compute the vector magnitude
int v_mag6( double *v_mag,
        double v_in[6]);

// Compute the normal vector (3 - dimension)
int v_norm( double u_out[3],
            double v_in[3]);

// Compute the normal vector (6 - dimension)
int v_norm6( double u_out[6],
         double v_in[6]);

// Compute the dot product
int v_dot(  double *dot,
            double v1[3],
            double v2[3]);

// Compute the cross product
int v_cross(    double v_cross[3],
                double v1[3],
                double v2[3]);

// Scale a vector
int v_scale(    double v_out[3],
                double v_in[3],
                double scale);

// Subtract vectors
int v_sub(  double v_out[3],
            double v1[3],
            double v2[3]);

// Add vectors
int v_add(  double v_out[3],
            double v1[3],
            double v2[3]);

// Copy vectors
int v_copy( double v_out[3],
            double v_in[3]);

// Copy matrices
int m_copy( double m_out[3][3],
            double m_in[3][3]);

// Copy matrices
int m_copy6( double m_out[6][6],
         double m_in[6][6]);

// Matrix x Vector
int m_x_v(  double v_out[3],
            double m_in[3][3],
            double v_in[3]);

// Matrix x Vector
int m_x_v6(  double v_out[6],
            double **m_in,
         double v_in[6]);

int m_x_v6bis(  double v_out[6],
            double m_in[6][6],
         double v_in[6]);

// Matrix x Vector
int m_x_v9(  double v_out[9],
            double **m_in,
         double v_in[9]);

// Matrix x Matrix (CBV 07/24/2015)
int m_x_m( double m_out[3][3],
       double m_in1[3][3],
       double m_in2[3][3] );

// Matrix x Matrix (CBV 10/04/2016)
int m_x_m6( double **m_out,
       double **m_in1,
        double **m_in2 );

int m_x_m6bis( double m_out[6][6],
       double m_in1[6][6],
              double m_in2[6][6] );

// Matrix Transpose
int m_trans(    double m_out[3][3],
                double m_in[3][3]);

// Matrix Transpose
int m_trans6(    double m_out[6][6],
         double m_in[6][6]);

// Coordinate Transformation Functions

// Equinoctial to Classical
int equin_to_class(double *sma, double *ecc, double *inc, double *raan, double *arg_per, double *true_ano, double *mean_ano, double af, double ag, double l, double n, double chi, double psi, double mu, double et, double fr); // chi looks like a x

// Compute derivative matrix elements equinoctial to classical (1)
int compute_T_deriv_equin_to_class(double T_equin_2_class[6][6], double af, double ag, double l, double n, double chi, double psi, double mu, double fr);

// Cartesian to Equinoctial (1)
int cart_to_equin( double *af, double *ag, double *l, double *n, double *chi, double *psi, double mu,  double fr, double rvec[3], double vvec[3]);

// Equinoctial to Cartesian (1)
int equin_to_cart(double rvec[3], double vvec[3],  double af, double ag, double l, double n, double chi, double psi, double mu, double fr);

// Compute derivative matrix elements equinoctial to cartesian (1)
int compute_T_deriv_equin_to_cart( double T_equin_to_cart[6][6], double af, double ag, double l, double n, double chi, double psi, double mu, double fr);

// Compute derivative matrix elements classical to cartesian (1)
//(page 28-31 of http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.112.5780&rep=rep1&type=pdf)
int compute_T_deriv_class_to_cart( double T_class_to_cart[6][6], double a, double e, double i, double o, double w, double nu, double mu );

// Inertial to NTW (CBV) (source: Vallado3 page 172)
int compute_T_inrtl_2_ntw( double T_inrtl_2_ntw[3][3],
               double r[3],
               double v[3]);

// compute transformation matrix from inertial to ntw but 6 by 6 dimensions. this is used to convert inertial covariance matricies *position, velocity) from inertial to ntw coordinates.
// important: contraritly to compute_T_inrtl_2_ntw, the outputs is ntw (twn in compute_T_inrtl_2_ntw). this is to compare to Vallado matlab's script covct2o2.m (in Code/cygnss/collision/vallado/matlab or email from Vallado o July 9, 2018)
// source: Vallado matlab's script covct2o2.m (in Code/cygnss/collision/vallado/matlab or email from Vallado o July 9, 2018)
int compute_T_inrtl_2_ntw_6by6( double T_inrtl_2_ntw_6by6[6][6],
               double r[3],
                double v[3]);

int q_copy(double q_out[4], double q_in[4]);

//  Inertial to LVLH http://degenerateconic.com/wp-content/uploads/2015/03/lvlh.pdf
int compute_T_inrtl_2_lvlh( double T_inrtl_2_lvlh[3][3],
                            double r[3],
                            double v[3]);

/* Inertial to Earth pressure frame. The Earth pressure frame is defined as follow: */
/*  Notations: */
/* - O center of the Earth */
/* - C position of the satellite */
/* - H projection of the satellite lcoation on the surface of the Earth (ie sub-satellite point) */
/* - S position of the Sun */
// - the z vector is the direction OH
// - the y vector is the vector in the plane OHS and in the direction of the Sun
// - the x vector completes the orthogonal basis
int compute_T_inrtl_2_earth_pres_frame( double T_inrtl_2_earth_pres_frame[3][3],
                    double r[3],
                    double et);

int compute_T_sc_to_lvlh(double T_sc_to_lvlh[3][3],
             double v_angle[3],
             int order_rotation[3],
             char attitude_profile[256],
             double *et,
             double r_i2cg_INRTL[3],
             double v_i2cg_INRTL[3],
             int file_is_quaternion,
             double quaternion[4],
             PARAMS_T *PARAMS
             );

//https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
int compute_T_enu_to_ecef( double T_enu_to_ecef[3][3], // out: rotation matrix ENU to ECEF
               double geodetic_latitude,  // in: geodetic latitude of the place at the surface of the Earth
               double longitude,   // in: longitude of the place at the surface of the Earth
               double flattening); // in: flattening parameter

int etprint( double et_to_print, char str_print[256] );
 
// FOOTNOTE:
// (1): (Vallado and Alfano 2015 "Updated Analytical Partials for Covariance Transformations and Optimization" (a few typos) and Danielson et al. 1995 "Semianalytic Satellite Theory" (that's how i found the typos in Alfano. Alfano numerical results are correct though but there are some typos in the expresssions of the equations)

#endif
