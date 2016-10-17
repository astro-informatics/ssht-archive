// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SSHT_ADJOINT
#define SSHT_ADJOINT

#include <complex.h>


void ssht_adjoint_mw_inverse_sov_sym(_Complex double *flm, 
				     _Complex double *f, 
				     int L, int spin, 
				     ssht_dl_method_t dl_method,
				     int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_real(_Complex double *flm, 
					  double *f, 
					  int L,
					  ssht_dl_method_t dl_method, 
					  int verbosity);
void ssht_adjoint_mw_forward_sov_sym(_Complex double *f, 
				     _Complex double *flm,
				     int L, int spin,
				     ssht_dl_method_t dl_method,
				     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_real(double *f, 
					  _Complex double *flm,
					  int L,
					  ssht_dl_method_t dl_method,
					  int verbosity);

void ssht_adjoint_mw_inverse_sov_sym_pole(_Complex double *flm, _Complex double *f,
					  _Complex double f_sp, double phi_sp,
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_real_pole(_Complex double *flm, 
					       double *f, 
					       double f_sp,
					       int L, 
					       ssht_dl_method_t dl_method,
					       int verbosity);
void ssht_adjoint_mw_forward_sov_sym_pole(_Complex double *f, 
					  _Complex double *f_sp, double *phi_sp,
					  _Complex double *flm, 
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_adjoint_mw_forward_sov_sym_real_pole(double *f, 
					       double *f_sp,
					       _Complex double *flm, 
					       int L, 
					       ssht_dl_method_t dl_method, 
					       int verbosity);


void ssht_adjoint_mw_inverse_sov_sym_ss(_Complex double *flm, _Complex double *f, 
					int L, int spin, 
					ssht_dl_method_t dl_method,
					int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_ss_real(_Complex double *flm, double *f, 
					     int L, 
					     ssht_dl_method_t dl_method, 
					     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss(_Complex double *f, _Complex double *flm,
					int L, int spin,
					ssht_dl_method_t dl_method,
					int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_real(double *f, 
					     _Complex double *flm,
					     int L,
					     ssht_dl_method_t dl_method,
					     int verbosity);

void ssht_adjoint_mw_inverse_sov_sym_ss_pole(_Complex double *flm, _Complex double *f,
					     _Complex double f_np, double phi_np,
					     _Complex double f_sp, double phi_sp,
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_ss_real_pole(_Complex double *flm, 
						  double *f, 
						  double f_np,
						  double f_sp,
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_pole(_Complex double *f, 
					     _Complex double *f_np, double *phi_np,
					     _Complex double *f_sp, double *phi_sp,
					     _Complex double *flm, 
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_real_pole(double *f, 
						  double *f_np,
						  double *f_sp,
						  _Complex double *flm, 
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity);


#endif
