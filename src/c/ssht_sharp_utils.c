#include "ssht_types.h"
#include "ssht_sharp_utils.h"
#include <string.h>
#include <fftw3.h>

typedef double complex dcmplx;

int ssht_use_libsharp_fwd(void)
  {
  static int res=-1;
  if (res==-1)
    {
#pragma omp critical (SSHT_SHARP_GETENV)
{
    char *env=getenv("SSHT_USE_LIBSHARP");
    if (!env)
      res=0;
    else
      res=(strcmp(env,"FWD")==0)||(strcmp(env,"BOTH")==0);
}
    }
  return res;
  }

int ssht_use_libsharp_inv(void)
  {
  static int res=-1;
  if (res==-1)
    {
#pragma omp critical (SSHT_SHARP_GETENV)
{
    char *env=getenv("SSHT_USE_LIBSHARP");
    if (!env)
      res=0;
    else
      res=(strcmp(env,"INV")==0)||(strcmp(env,"BOTH")==0);
}
    }
  return res;
  }

void ssht_flm2alm_r (const dcmplx *flm, int L0, int L, dcmplx ***alm,
  sharp_alm_info **ainfo)
  {
  sharp_make_triangular_alm_info(L-1,L-1,1,ainfo);
  ALLOC2D((*alm),dcmplx,1,(L*(L+1))/2);
  // rearrange a_lm into required format
  int l,m;
  for (m=0; m<L; ++m)
    for (l=m; l<L; ++l)
      (*alm)[0][sharp_alm_index(*ainfo,l,m)]=flm[l*l + l + m];
  for (m=0; m<L0; ++m)
    for (l=m; l<L0; ++l)
      (*alm)[0][sharp_alm_index(*ainfo,l,m)]=0.;
  }

void ssht_flm2alm_c (const dcmplx *flm, int L0, int L, int spin,
  dcmplx ***alm, sharp_alm_info **ainfo)
  {
  sharp_make_triangular_alm_info(L-1,L-1,1,ainfo);
  ALLOC2D((*alm),dcmplx,2,(L*(L+1))/2);
  // rearrange a_lm into required format
  double spinsign = (spin==0) ? 1. : -1.;
  double spinfac = (spin&1) ? -1. : 1.;
  for (int m=0; m<L; ++m)
    {
    double mfac=(m&1) ? -1.:1.;
    for (int l=m; l<L; ++l)
      {
      dcmplx pslm=flm[l*l+l+m], mslm=spinfac*mfac*conj(flm[l*l+l-m]);
      if (spin<0) SWAP(pslm,mslm,dcmplx);
      dcmplx E=spinsign*0.5*(pslm+spinfac*mslm);
      dcmplx B=-spinsign*0.5*_Complex_I*(pslm-spinfac*mslm);
      (*alm)[0][sharp_alm_index(*ainfo,l,m)]=E;
      (*alm)[1][sharp_alm_index(*ainfo,l,m)]=B;
      }
    }
  for (int m=0; m<L0; ++m)
    for (int l=m; l<L0; ++l)
      (*alm)[0][sharp_alm_index(*ainfo,l,m)]=
      (*alm)[1][sharp_alm_index(*ainfo,l,m)]= 0.;
  }

void ssht_alm2flm_c (dcmplx *flm, int L0, int L, int spin, dcmplx **alm,
  sharp_alm_info *ainfo)
  {
  double spinsign = (spin==0) ? 1. : -1.;
  double spinfac = (spin&1) ? -1. : 1.;
  for (int m=0; m<L; ++m)
    {
    double mfac=(m&1) ? -1:1;
    for (int l=m; l<L; ++l)
      {
      dcmplx E=alm[0][sharp_alm_index(ainfo,l,m)],
             B=alm[1][sharp_alm_index(ainfo,l,m)];
      if (spin>=0)
        {
        flm[l*l+l+m]=spinsign*(E+_Complex_I*B);
        flm[l*l+l-m]=spinsign*(mfac*(conj(E)+_Complex_I*conj(B)));
        }
      else
        {
        flm[l*l+l-m]=spinfac*mfac*conj(spinsign*(E+_Complex_I*B));
        flm[l*l+l+m]=spinfac*conj(spinsign*(conj(E)+_Complex_I*conj(B)));
        }
      }
    }
  for (int m=0; m<L0; ++m)
    for (int l=m; l<L0; ++l)
      flm[l*l+l+m] = flm[l*l+l-m] = 0.;
  }

void ssht_alm2flm_r (dcmplx *flm, int L0, int L, dcmplx **alm,
  sharp_alm_info *ainfo)
  {
  for (int m=0; m<L; ++m)
    for (int l=m; l<L; ++l)
      {
      flm[l*l + l + m]=alm[0][sharp_alm_index(ainfo,l,m)];
      flm[l*l + l - m]=conj(flm[l*l + l + m])*((m&1)? -1:1);
      }
  for (int m=0; m<L0; ++m)
    for (int l=m; l<L0; ++l)
      flm[l*l+l+m] = flm[l*l+l-m] = 0.;
  }

void ssht_sharp_mw_forward_complex(dcmplx *flm, const dcmplx *f,
  int L0, int L, int spin)
  {
  int nphi=2*L-1;
  int nm=nphi;
  int nth_mw=L;
  int nth_mwfull=2*L-1;
  int nth_hw=2*L;
  dcmplx **tmp1;
  ALLOC2D(tmp1,dcmplx,nth_hw,nphi);

  // FFT in phi direction
  {
  fftw_plan plan = fftw_plan_dft_1d(nphi,tmp1[0],tmp1[1],
    FFTW_FORWARD,FFTW_MEASURE|FFTW_UNALIGNED);
#pragma omp parallel for
  for (int ith=0; ith<nth_mw; ++ith)
    fftw_execute_dft(plan,(dcmplx *)(f+ith*nphi),tmp1[ith]);
  fftw_destroy_plan(plan);
  }

  {
  double norm=1./(nth_mwfull*nphi);
  double dtheta=-SSHT_PI/nth_mwfull;
  dcmplx *fact=RALLOC(dcmplx,nth_mw);
  for (int ith=0; ith<nth_mw; ++ith)
    fact[ith] = norm*cexp(_Complex_I*ith*dtheta);
  dcmplx *tmp=fftw_alloc_complex(nth_mwfull);
  fftw_plan plan1 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan plan2 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
  fftw_free(tmp);
  // loop over all m
#pragma omp parallel
{
  dcmplx *tmp=fftw_alloc_complex(nth_mwfull);
#pragma omp for
  for (int m=0; m<nm; ++m)
    {
    for (int ith=0; ith<nth_mw; ++ith)
      tmp[ith]= tmp1[ith][m];

    // theta extension
    int mreal= (m<L) ? m : m-(2*L-1);
    mreal+=spin;
    int sign = (mreal&1) ? -1. : 1.;
    for (int ith=nth_mw; ith<nth_mwfull; ++ith)
      tmp[ith]= sign*tmp[nth_mwfull-1-ith];

    fftw_execute_dft(plan1,tmp,tmp);

    // theta shift and normalisation
    tmp[0]*=norm;
    for (int ith=1; ith<nth_mw; ++ith)
      {
      tmp[ith] *= fact[ith];
      tmp[nth_mwfull-ith] *= conj(fact[ith]);
      }

    fftw_execute_dft(plan2,tmp,tmp);

    for (int ith=0; ith<nth_mw; ++ith)
      tmp1[2*ith][m]=tmp[ith];
    }
  fftw_free(tmp);
} // end of parallel region
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(fact);
  }

  // FFT in phi direction
  {
  dcmplx *ttt=RALLOC(dcmplx,nphi);
  fftw_plan plan = fftw_plan_dft_1d(nphi,ttt,ttt,
    FFTW_BACKWARD,FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
#pragma omp parallel for
  for (int ith=0; ith<nth_mw; ++ith)
    fftw_execute_dft(plan,tmp1[2*ith],tmp1[2*ith]);
  fftw_destroy_plan(plan);
  }
  // copy original map data
  for (int ith=0; ith<nth_mw; ++ith)
    for (int m=0; m<nphi; ++m)
      tmp1[2*ith+1][m]=f[ith*nphi+m];

  if (spin<0)
    for (int ith=0; ith<nth_hw; ++ith)
      for (int m=0; m<nphi; ++m)
        tmp1[ith][m]=conj(tmp1[ith][m]);

  sharp_geom_info *tinfo;
  sharp_make_cc_geom_info (nth_hw, nphi, 0., 2, 2*nphi, &tinfo);
  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(L-1,L-1,1,&alms);
  dcmplx **alm;
  ALLOC2D(alm,dcmplx,2,(L*(L+1))/2);
  double * fr=(double *)(tmp1[0]);
  double *frp[2];
  frp[0]=fr;
  frp[1]=fr+1;
  sharp_execute(SHARP_MAP2ALM,abs(spin),alm,&frp[0],tinfo,alms,(spin==0)?2:1,SHARP_DP,NULL,NULL);
  ssht_alm2flm_c(flm,L0,L,spin,alm,alms);
  DEALLOC2D(alm);
  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  DEALLOC2D(tmp1);
  }

void ssht_sharp_mw_forward_real(dcmplx *flm, const double *f, int L0, int L)
  {
  int nphi=2*L-1;
  int nm=L;
  int nth_mw=L;
  int nth_mwfull=2*L-1;
  int nth_hw=2*L;
  double **tmp1;
  ALLOC2D(tmp1,double,nth_hw,2*nm);

  // FFT in phi direction
  {
  fftw_plan plan = fftw_plan_dft_r2c_1d(nphi,tmp1[0],(dcmplx *)tmp1[1],
    FFTW_MEASURE|FFTW_UNALIGNED);
#pragma omp parallel for
  for (int ith=0; ith<nth_mw; ++ith)
    fftw_execute_dft_r2c(plan,(double *)(f+ith*nphi),(dcmplx *)tmp1[ith]);
  fftw_destroy_plan(plan);
  }

  {
  double norm=1./(nth_mwfull*nphi);
  double dtheta=-SSHT_PI/nth_mwfull;
  dcmplx *fact=RALLOC(dcmplx,nth_mw);
  for (int ith=0; ith<nth_mw; ++ith)
    fact[ith] = norm*cexp(_Complex_I*ith*dtheta);
  dcmplx *tmp=fftw_alloc_complex(nth_mwfull);
  fftw_plan plan1 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan plan2 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
  fftw_free(tmp);
  // loop over all m
#pragma omp parallel
{
  dcmplx *tmp=fftw_alloc_complex(nth_mwfull);
#pragma omp for
  for (int m=0; m<nm; ++m)
    {
    for (int ith=0; ith<nth_mw; ++ith)
      tmp[ith]= tmp1[ith][2*m]+_Complex_I*tmp1[ith][2*m+1];

    // theta extension
    int sign = (m&1) ? -1. : 1.;
    for (int ith=nth_mw; ith<nth_mwfull; ++ith)
      tmp[ith]= sign*tmp[nth_mwfull-1-ith];

    fftw_execute_dft(plan1,tmp,tmp);

    // theta shift and normalisation
    tmp[0]*=norm;
    for (int ith=1; ith<nth_mw; ++ith)
      {
      tmp[ith] *= fact[ith];
      tmp[nth_mwfull-ith] *= conj(fact[ith]);
      }

    fftw_execute_dft(plan2,tmp,tmp);

    for (int ith=0; ith<nth_mw; ++ith)
      {
      tmp1[2*ith][2*m]=creal(tmp[ith]);
      tmp1[2*ith][2*m+1]=cimag(tmp[ith]);
      }
    }
  fftw_free(tmp);
} // end of parallel region
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(fact);
  }
  // FFT in phi direction
  {
  dcmplx *ttt=RALLOC(dcmplx,nm);
  fftw_plan plan = fftw_plan_dft_c2r_1d(nphi,ttt,(double *)ttt,
    FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
#pragma omp parallel for
  for (int ith=0; ith<nth_mw; ++ith)
    fftw_execute_dft_c2r(plan,(dcmplx *)tmp1[2*ith],tmp1[2*ith]);
  fftw_destroy_plan(plan);
  }
  // copy original map data
  for (int ith=0; ith<nth_mw; ++ith)
    for (int m=0; m<nphi; ++m)
      tmp1[2*ith+1][m]=f[ith*nphi+m];

  sharp_geom_info *tinfo;
  sharp_make_cc_geom_info (nth_hw, nphi, 0., 1, 2*nm, &tinfo);
  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(L-1,L-1,1,&alms);
  dcmplx **alm;
  ALLOC2D(alm,dcmplx,1,(L*(L+1))/2);
  sharp_execute(SHARP_MAP2ALM,0,alm,tmp1,tinfo,alms,1,SHARP_DP,NULL,NULL);
  ssht_alm2flm_r(flm,L0,L,alm,alms);
  DEALLOC2D(alm);
  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  DEALLOC2D(tmp1);
  }

void ssht_sharp_mws_forward_complex(dcmplx *flm, const dcmplx *f, int L0, int L, int spin)
  {
  int nphi=2*L;
  int nth_mw=L+1;
  int nth_mwfull=2*L;
  int nth_hw=2*L+1;
  dcmplx **tmp1;
  ALLOC2D(tmp1,dcmplx,nth_hw,nphi);

  {
  double norm=1./nth_mwfull;
  double dtheta=SSHT_PI/nth_mwfull;
  dcmplx *fact=RALLOC(dcmplx,nth_mw);
  for (int ith=0; ith<nth_mw; ++ith)
    fact[ith] = norm*cexp(_Complex_I*ith*dtheta);
  dcmplx *tmp=fftw_alloc_complex(nth_mwfull);
  fftw_plan plan1 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan plan2 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
  fftw_free(tmp);
  // loop over all m
#pragma omp parallel
{
  dcmplx *tmp=fftw_alloc_complex(nth_mwfull);
#pragma omp for
  for (int m=0; m<nphi; ++m)
    {
    for (int ith=0; ith<nth_mw; ++ith)
      tmp[ith]= f[ith*nphi+m];

    double nsign = (spin&1) ? -1. : 1.;
    int m_opposite=(m+nphi/2)%nphi;
    for (int ith=nth_mw; ith<nth_mwfull; ++ith)
      tmp[ith]= nsign*f[nphi*(nth_mwfull-ith)+m_opposite];

    fftw_execute_dft(plan1,tmp,tmp);

    // theta shift and normalisation
    tmp[0]*=norm;
    for (int ith=1; ith<nth_mw-1; ++ith)
      {
      tmp[ith] *= fact[ith];
      tmp[nth_mwfull-ith] *= conj(fact[ith]);
      }
    tmp[nth_mw-1]*=norm;

    fftw_execute_dft(plan2,tmp,tmp);

    for (int ith=0; ith<nth_mw-1; ++ith)
      tmp1[2*ith+1][m]=tmp[ith];
    }
  fftw_free(tmp);
} // end of parallel region
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(fact);
  }

  // copy original map data
  for (int ith=0; ith<nth_mw; ++ith)
    for (int m=0; m<nphi; ++m)
      tmp1[2*ith][m]=f[ith*nphi+m];

  if (spin<0)
    for (int ith=0; ith<nth_hw; ++ith)
      for (int m=0; m<nphi; ++m)
        tmp1[ith][m]=conj(tmp1[ith][m]);

  sharp_geom_info *tinfo;
  sharp_make_cc_geom_info (nth_hw, nphi, 0., 2, 2*nphi, &tinfo);
  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(L-1,L-1,1,&alms);
  double * fr=(double *)(tmp1[0]);
  double *frp[2];
  frp[0]=fr;
  frp[1]=fr+1;
  dcmplx **alm;
  ALLOC2D(alm,dcmplx,2,(L*(L+1))/2);
  sharp_execute(SHARP_MAP2ALM,abs(spin),alm,&frp[0],tinfo,alms,(spin==0)?2:1,SHARP_DP,NULL,NULL);
  ssht_alm2flm_c(flm,L0,L,spin,alm,alms);
  DEALLOC2D(alm);
  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  DEALLOC2D(tmp1);
  }

void ssht_sharp_mws_forward_real(dcmplx *flm, const double *f, int L0, int L)
  {
  int nphi=2*L;
  int nth_mw=L+1;
  int nth_mwfull=2*L;
  int nth_hw=2*L+1;
  double **tmp1;
  ALLOC2D(tmp1,double,nth_hw,nphi);

  {
  double norm=1./nth_mwfull;
  double dtheta=SSHT_PI/nth_mwfull;
  dcmplx *fact=RALLOC(dcmplx,nth_mw);
  for (int ith=0; ith<nth_mw; ++ith)
    fact[ith] = norm*cexp(_Complex_I*ith*dtheta);
  dcmplx *tmp=fftw_alloc_complex(nth_mwfull);
  fftw_plan plan1 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan plan2 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
  fftw_free(tmp);
  // loop over all m
#pragma omp parallel
{
  dcmplx *tmp=fftw_alloc_complex(nth_mwfull);
#pragma omp for
  for (int m=0; m<nphi; ++m)
    {
    for (int ith=0; ith<nth_mw; ++ith)
      tmp[ith]= f[nphi*ith+m];
    int m_opposite=(m+nphi/2)%nphi;
    for (int ith=nth_mw; ith<nth_mwfull; ++ith)
      tmp[ith]= f[nphi*(nth_mwfull-ith)+m_opposite];

    fftw_execute_dft(plan1,tmp,tmp);

    // theta shift and normalisation
    tmp[0]*=norm;
    for (int ith=1; ith<nth_mw-1; ++ith)
      {
      tmp[ith] *= fact[ith];
      tmp[nth_mwfull-ith] *= conj(fact[ith]);
      }
    tmp[nth_mw-1]*=norm;

    fftw_execute_dft(plan2,tmp,tmp);

    for (int ith=0; ith<nth_mw-1; ++ith)
      tmp1[2*ith+1][m]=creal(tmp[ith]);
    }
  fftw_free(tmp);
} // end of parallel region
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(fact);
  }

  // copy original map data
  for (int ith=0; ith<nth_mw; ++ith)
    for (int m=0; m<nphi; ++m)
      tmp1[2*ith][m]=f[ith*nphi+m];

  sharp_geom_info *tinfo;
  sharp_make_cc_geom_info (nth_hw, nphi, 0., 1, nphi, &tinfo);
  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(L-1,L-1,1,&alms);
  dcmplx **alm;
  ALLOC2D(alm,dcmplx,1,(L*(L+1))/2);
  sharp_execute(SHARP_MAP2ALM,0,alm,tmp1,tinfo,alms,1,SHARP_DP,NULL,NULL);
  ssht_alm2flm_r(flm,L0,L,alm,alms);
  DEALLOC2D(alm);
  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  DEALLOC2D(tmp1);
  }
