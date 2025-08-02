/******************************************************************************************************/
/*                           HYREC: Hydrogen and Helium Recombination Code                            */
/*                      Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                      */
/*                                                                                                    */
/*         history.c: functions for recombination history                                             */
/*                                                                                                    */
/*         Version: October 2012                                                                      */
/*                                                                                                    */
/*         Revision history:                                                                          */
/*            - written November 2010                                                                 */
/*            - January 2011: changed various switches (notably for post-Saha expansions)             */
/*                             so that they remain valid for arbitrary cosmologies                    */
/*            - May 2012:   - added explicit dependence on fine structure constant and electron mass  */
/*                          - modified call of rec_build_history                                      */
/*                             and improved numerical radiative transfer equations                    */
/*                             so the Lyman-lines spectrum can be extracted                           */
/*                           - split some functions for more clarity                                  */
/*            - October 2012: added some wrapper functions for running CAMB with HyRec                */
/*            - March 2023: added module for DM&PBH energy injection as HYREC-2 doen not seem to work */
/*                            with camb                                                               */
/*                             (courtesy of Antony Lewis)                                             */
/******************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include "hyrectools.h"
#include "helium.h"
#include "hydrogen.h"
#include "history.h"
#include "Interpolate_EFF.h"
#include "DarkSide.h"

/*****************************************************************************
Setting derived cosmological parameters needed for recombination calculation
******************************************************************************/

void rec_set_derived_params(REC_COSMOPARAMS *param)
{
    /* Separated by AML to avoid duplication */

    double z, Pion, Tresc, RLya, four_betaB;

    param->nH0 = 11.223846333047 * param->obh2 * (1. - param->Y); /* number density of hydrogen today in m-3 */
    param->odmh2 = param->omh2 - param->obh2;
    param->fHe = param->Y / (1 - param->Y) / 3.97153; /* abundance of helium by number */
    /* these should depend on fsR and meR, strictly speaking; however, these are corrections to corrections */

    /* Total number of redshift steps */
    param->nz = (long)floor(2 + log((1. + ZSTART) / (1. + ZEND)) / DLNA);

    /* (Added May 2012)
       In order to save memory for the radiation field tables, compute right away
       the range of redshifts at which radiative transfer is actually followed,
       and the corresponding number of steps.
       Radiative transfer starts when Tr < TR_MAX and stops when the ionization probablity
       from n=2 (estimated with simple Peebles model) is less than PION_MAX. */

    param->izH0 = (long)floor(1 + log(kBoltz * param->T0 / square(param->fsR) / param->meR * (1. + ZSTART) / TR_MAX) / DLNA);
    param->zH0 = (1. + ZSTART) * exp(-param->izH0 * DLNA) - 1.;

    Pion = 1.;
    z = 900.;
    while (Pion > PION_MAX && z > ZEND)
    {
        Tresc = kBoltz * param->T0 * (1. + z) / param->fsR / param->fsR / param->meR; /* Tr rescaled for alpha, me */
        RLya = LYA_FACT(param->fsR, param->meR) * rec_HubbleConstant(param, z) / (1e-6 * param->nH0 * cube(1. + z));
        four_betaB = SAHA_FACT(param->fsR, param->meR) * Tresc * sqrt(Tresc) * exp(-0.25 * EI / Tresc) * alphaB_PPB(Tresc, param->fsR, param->meR);
        Pion = four_betaB / (3. * RLya + L2s_rescaled(param->fsR, param->meR) + four_betaB);
        z -= 10.;
    }

    param->nzrt = (long)floor(2 + log((1. + ZSTART) / (1. + z)) / DLNA) - param->izH0;
}

/*************************************************************************************
Hubble expansion rate in sec^-1.
Note: neutrinos are assumed to be massless.
*************************************************************************************/

#ifdef CAMB

extern double dtauda_(double *);

double rec_HubbleConstant(REC_COSMOPARAMS *param, double z)
{
    double a;

    a = 1. / (1. + z);
    /* conversion from d tau/ da in Mpc to H(z) in 1/s */
    return 1. / (a * a) / dtauda_(&a) / 3.085678e22 * 2.99792458e8;
}

HRATEEFF rate_table;
TWO_PHOTON_PARAMS twog_params;
REC_COSMOPARAMS param;
int firstTime = 0;
double *xe_output, *tm_output, **Dfnu_hist, *Dfminus_Ly_hist[3];
double logstart;

void hyrec_init()
{

    /* Build effective rate table */
    char *buffer = (char *)malloc(1024);
    getcwd(buffer, 1024);
    chdir(HYRECPATH);
    rate_table.logTR_tab = create_1D_array(NTR);
    rate_table.TM_TR_tab = create_1D_array(NTM);
    rate_table.logAlpha_tab[0] = create_2D_array(NTM, NTR);
    rate_table.logAlpha_tab[1] = create_2D_array(NTM, NTR);
    rate_table.logR2p2s_tab = create_1D_array(NTR);
    read_rates(&rate_table);

    /* Read two-photon rate tables */
    read_twog_params(&twog_params);
    chdir(buffer);
    free(buffer);
}

void rec_build_history_camb_(const double *OmegaC, const double *OmegaB, const double *OmegaN,
                             const double *Omegav, const double *h0inp, const double *tcmb, const double *yp, const double *num_nu,
                             const double *DM_Channel, const double *Mdm, const double *Pann, const double *Gamma, const double *PBH_Model,
                             const double *PBH_Distribution, const double *Mbh, const double *fbh, const double *PBH_Lognormal_Sigma,
                             const double *PBH_PWL_Mmax, const double *PBH_PWL_Gamma, const double *PBH_Spin)
{

    double h2 = *h0inp / 100.;
    h2 = h2 * h2;
    param.T0 = *tcmb;
    param.obh2 = *OmegaB * h2;
    param.omh2 = (*OmegaB + *OmegaC) * h2;
    param.okh2 = (1 - *OmegaC - *OmegaB - *Omegav - *OmegaN) * h2;
    param.odeh2 = *Omegav * h2;
    param.w0 = -1; /* not actually used */
    param.wa = 0;
    param.Y = *yp;
    param.Nnueff = *num_nu;
    param.fsR = param.meR = 1.; /*** Default: today's values ***/

    // Passing DM params from camb to HyRec
    param.DM_Channel = *DM_Channel;
    param.Mdm = *Mdm;
    param.Pann = *Pann;
    param.Gamma = *Gamma;
    param.PBH_Model = *PBH_Model;
    param.PBH_Distribution = *PBH_Distribution;
    param.Mbh = *Mbh;
    param.fbh = *fbh;
    param.PBH_Lognormal_Sigma = *PBH_Lognormal_Sigma;
    param.PBH_PWL_Mmax = *PBH_PWL_Mmax;
    param.PBH_PWL_Gamma = *PBH_PWL_Gamma;
    param.PBH_Spin = *PBH_Spin;
    param.Use_SSCK = 0; // Do this later for camb
    
    // Set default camb params for DM&PBH
    param.odmh2 = param.omh2 - param.obh2;

    rec_set_derived_params(&param);

    if (firstTime == 0)
    {
        hyrec_init();
        firstTime = 1;
        logstart = -log(1. + ZSTART);
        xe_output = create_1D_array(param.nz);
        tm_output = create_1D_array(param.nz);
    }
    Dfnu_hist = create_2D_array(NVIRT, param.nzrt);
    Dfminus_Ly_hist[0] = create_1D_array(param.nzrt); /* Ly-alpha */
    Dfminus_Ly_hist[1] = create_1D_array(param.nzrt); /* Ly-beta  */
    Dfminus_Ly_hist[2] = create_1D_array(param.nzrt); /* Ly-gamma */
    rec_build_history(&param, &rate_table, &twog_params, xe_output, tm_output, Dfnu_hist, Dfminus_Ly_hist);
    free_2D_array(Dfnu_hist, NVIRT);
    free(Dfminus_Ly_hist[0]);
    free(Dfminus_Ly_hist[1]);
    free(Dfminus_Ly_hist[2]);
}

double hyrec_xe_(double *a)
{
    double loga = log(*a);
    if (loga < logstart)
    {
        return xe_output[0];
    }
    return rec_interp1d(logstart, DLNA, xe_output, param.nz, loga);
}

double hyrec_tm_(double *a)
{
    double loga = log(*a);
    if (loga < logstart)
    {
        return tm_output[0];
    }
    return rec_interp1d(logstart, DLNA, tm_output, param.nz, loga);
}

#else

double rec_HubbleConstant(REC_COSMOPARAMS *param, double z)
{

    int j;
    double rho, rho_i[5], ainv;
    double ogh2;

    ainv = 1. + z; /* inverse scale factor */

    /* Matter */
    rho_i[0] = param->omh2 * ainv * ainv * ainv;

    /* Curvature */
    rho_i[1] = param->okh2 * ainv * ainv;

    /* Dark energy */
    rho_i[2] = param->odeh2 * pow(ainv, 3 * (1 + param->w0)) * exp(3 * param->wa * (log(ainv) - 1. + 1. / ainv));

    /* Get radiation density.
     * Coefficient based on 1 AU = 1.49597870691e11 m (JPL SSD) */
    ogh2 = 4.48162687719e-7 * param->T0 * param->T0 * param->T0 * param->T0;
    rho_i[3] = ogh2 * ainv * ainv * ainv * ainv;

    /* Neutrino density -- scaled from photons assuming lower temperature */
    rho_i[4] = 0.227107317660239 * rho_i[3] * param->Nnueff;

    /* Total density, including curvature contribution */
    rho = 0.;
    for (j = 0; j < 5; j++)
        rho += rho_i[j];

    /* Conversion to Hubble rate */
    return (3.2407792896393e-18 * sqrt(rho));
}
#endif

/*************************************************************************************************
Cosmological parameters Input/Output
*************************************************************************************************/

void rec_get_cosmoparam(FILE *fin, FILE *fout, REC_COSMOPARAMS *param)
{
    char Param_Name[100];

    /*-------- Dark Matter params --------*/
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%d", &(param->Use_SSCK));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->DM_Channel));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->Mdm));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->Pann));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->Gamma));

    /*-------- Primordial Black Hole params --------*/
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->PBH_Model));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->PBH_Distribution));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->Mbh));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->fbh));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->PBH_Lognormal_Sigma));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->PBH_PWL_Mmax));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->PBH_PWL_Gamma));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->PBH_Spin));

    /*-------- Cosmological params --------*/

    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->T0));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->obh2));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->omh2));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->okh2));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->odeh2));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg %lg", &(param->w0), &(param->wa));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->Y));
    fscanf(fin, "%s", Param_Name);
    fscanf(fin, "%lg", &(param->Nnueff));

    // Ensure all params are correctly entered
    /*
    printf("Use_SSCK = %d\n", param->Use_SSCK);
    printf("DM_Channel = %f\n", param->DM_Channel);
    printf("Mdm = %f\n", param->Mdm);
    printf("Pann = %E\n", param->Pann);
    printf("Gamma = %E\n", param->Gamma);
    printf("PBH_Model = %f\n", param->PBH_Model);
    printf("PBH_Distribution = %f\n", param->PBH_Distribution);
    printf("Mbh = %E\n", param->Mbh);
    printf("fbh = %E\n", param->fbh);
    printf("PBH_Lognormal_Sigma = %f\n", param->PBH_Lognormal_Sigma);
    printf("PBH_PWL_Mmax = %E\n", param->PBH_PWL_Mmax);
    printf("PBH_PWL_Gamma = %f\n", param->PBH_PWL_Gamma);
    printf("PBH_Spin = %f\n", param->PBH_Spin);

    printf("Tcmb = %f\n", param->T0);
    printf("obh2 = %f\n", param->obh2);
    printf("omh2 = %f\n", param->omh2);
    printf("okh2 = %f\n", param->okh2);
    printf("odeh2 = %f\n", param->odeh2);
    printf("w0 = %f\n", param->w0);
    printf("wa = %f\n", param->wa);
    printf("YHe = %f\n", param->Y);
    printf("Neff = %f\n", param->Nnueff);
    */
    
    /****** Added May 2012: explicit dependence on fine-structure constant and electron mass ******/
    /** fsR = alpha_fs(rec) / alpha_fs(today), meR = me(rec) / me(today) **/

    param->fsR = param->meR = 1.0; /*** Default: today's values ***/

    /**** UNCOMMENT IF WANT TO USE DIFFERENT VALUES OF alpha_fs OR me *****/

    /* if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter ratio of fine structure constant at last scattering to today's value, alpha(rec)/alpha(now): "); */
    /* fscanf(fin, "%lg", &(param->fsR)); */
    /* if (fout!=NULL && PROMPT==1) fprintf(fout, "Enter ratio of electron mass at last scattering to today's value, me(rec)/me(now): "); */
    /* fscanf(fin, "%lg", &(param->meR));  */

    rec_set_derived_params(param);

    if (fout != NULL && PROMPT == 1)
        fprintf(fout, "\n");
}

/*****************************************************************************************
Matter temperature -- 1st order steady state, from Hirata 2008.
The input and output temperatures are in KELVIN.
******************************************************************************************/

double rec_Tmss(double xe, double Tr, double H, double fHe, double fsR, double meR, double *DarkArray)
{
    double coeff_inv, nH, dEdVdt_Heat, DM_Term, coeff;
    dEdVdt_Heat = DarkArray[2];
    nH = DarkArray[3];
    coeff_inv = H / (fsR * fsR / cube(meR) * 4.91466895548409e-22) / pow(Tr, 4.0) * (1. + xe + fHe) / xe;
    coeff = 1.0 / coeff_inv;
    DM_Term = dEdVdt_Heat / kBoltz / (1.5 * nH * (1. + xe + fHe)) / H / (1. + coeff);

    // return Tr / (1. + H / (fsR * fsR / cube(meR) * 4.91466895548409e-22) / pow(Tr, 4.0) * (1. + xe + fHe) / xe);
    return Tr / (1. + coeff_inv) + DM_Term;

    /* Coefficient = 8 sigma_T a_r / (3 m_e c) */
    /* Here Tr, Tm are the actual (not rescaled) temperatures */
}

/******************************************************************************************
Matter temperature evolution derivative. Input and output temperatures are in KELVIN.
Added May 2012: when Tm = Tr, return -Tr (needed by CLASS)
Modified by JSC: allow Tm > Tr
******************************************************************************************/

double rec_dTmdlna(double xe, double Tm, double Tr, double H, double fHe, double fsR, double meR, double *DarkArray)
{
    double Q_adia, Q_compt, Q_dm, z, nH, dEdVdt_Heat;
    Q_adia = -2. * Tm;
    Q_compt = fsR * fsR / meR / meR / meR * 4.91466895548409e-22 * Tr * Tr * Tr * Tr * xe / (1. + xe + fHe) * (Tr - Tm) / H;

    dEdVdt_Heat = DarkArray[2];
    nH = DarkArray[3];
    Q_dm = dEdVdt_Heat / kBoltz / (1.5 * nH * (1. + xe + fHe)) / H;

    // return (Tr / Tm - 1. < 1e-10 ? -Tr : Q_adia + Q_compt + Q_dm);
    return Q_adia + Q_compt + Q_dm;
    /* Here Tr, Tm are the actual (not rescaled) temperatures */
}

/**********************************************************************************************
Second-order integrator for HeII abundance during HeII->HeI recombination.
Assumes hydrogen ionization is in Saha equilibrium (free electrons being provided by both H and He).
If xHeII is close enough to Saha equilibrium do a post-Saha expansion.
May 2012: removed unused z_prev, z_prev2 variables
***********************************************************************************************/

void rec_get_xe_next1_He(REC_COSMOPARAMS *param, double z_in, double *xHeII,
                         double *dxHeIIdlna_prev, double *dxHeIIdlna_prev2, int *post_saha)
{

    double H, xH1s, xH1s_p, xH1s_m, xHeIISaha, dxHeIISaha_dlna, DdxHeIIdlna_Dxe, dxHeIIdlna, z_out, Dxe;

    H = rec_HubbleConstant(param, z_in);
    xH1s = rec_saha_xH1s(*xHeII, param->nH0, param->T0, z_in, param->fsR, param->meR);
    dxHeIIdlna = rec_helium_dxHeIIdlna(xH1s, *xHeII, param->nH0, param->T0, param->fHe, H, z_in, param->fsR, param->meR);

    /* Post-Saha approximation during the early phase of HeII->HeI recombination */
    if (*post_saha == 1)
    {
        z_out = (1. + z_in) * exp(-DLNA) - 1.;
        H = rec_HubbleConstant(param, z_out);
        xHeIISaha = rec_saha_xHeII(param->nH0, param->T0, param->fHe, z_out, param->fsR, param->meR);

        dxHeIISaha_dlna = (1. + z_out) * (rec_saha_xHeII(param->nH0, param->T0, param->fHe, z_out - 0.5, param->fsR, param->meR) - rec_saha_xHeII(param->nH0, param->T0, param->fHe, z_out + 0.5, param->fsR, param->meR));

        Dxe = 0.01 * (param->fHe - xHeIISaha);
        xH1s_p = rec_saha_xH1s(xHeIISaha + Dxe, param->nH0, param->T0, z_out, param->fsR, param->meR);
        xH1s_m = rec_saha_xH1s(xHeIISaha - Dxe, param->nH0, param->T0, z_out, param->fsR, param->meR);

        DdxHeIIdlna_Dxe = (rec_helium_dxHeIIdlna(xH1s_p, xHeIISaha + Dxe, param->nH0, param->T0, param->fHe, H, z_out, param->fsR, param->meR) - rec_helium_dxHeIIdlna(xH1s_m, xHeIISaha - Dxe, param->nH0, param->T0, param->fHe, H, z_out, param->fsR, param->meR)) / 2. / Dxe;

        *xHeII = xHeIISaha + dxHeIISaha_dlna / DdxHeIIdlna_Dxe;

        /* Check that the post-Saha approximation is still valid. If not, switch it off for future iterations */
        if (fabs(*xHeII - xHeIISaha) > DXHEII_MAX)
            *post_saha = 0;
    }

    /* Otherwise integrate ODE */
    else
        *xHeII += DLNA * (1.25 * dxHeIIdlna - 0.25 * (*dxHeIIdlna_prev2));

    /* Update stored derivatives */
    *dxHeIIdlna_prev2 = *dxHeIIdlna_prev;
    *dxHeIIdlna_prev = dxHeIIdlna;
}

/****************************************************************************************************
Post-Saha value of neutral hydrogen abundance at the redshift z_out.
Tm = Tr is assumed.
xHeII_out is the given value of Helium abundance at z_out (may be non-zero if Helium is still recombining)
iz_out is the index corresponding to z_out.
****************************************************************************************************/

double rec_xH1s_postSaha(REC_COSMOPARAMS *param, unsigned iz_out, double z_out, double xHeII_out,
                         HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                         double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist, int *post_saha, double *DarkArray)
{

    double ainv, xH1sSaha, xHIISaha, dxH1sSaha_dlna, dxH1sdlna_Saha, DdxH1sdlna_DxH1s, H, T, nH, Dxe, xH1s;

    xH1sSaha = rec_saha_xH1s(xHeII_out, param->nH0, param->T0, z_out, param->fsR, param->meR);
    xHIISaha = 1. - xH1sSaha;
    H = rec_HubbleConstant(param, z_out);
    T = kBoltz * param->T0 * (ainv = 1. + z_out); /* Convert to eV for hydrogen rec functions */
    nH = 1e-6 * param->nH0 * ainv * ainv * ainv;  /* Convert to cm^-3 for hydrogen rec functions */

    dxH1sSaha_dlna = (1. + z_out) * (rec_saha_xH1s(xHeII_out, param->nH0, param->T0, z_out - 0.5, param->fsR, param->meR) - rec_saha_xH1s(xHeII_out, param->nH0, param->T0, z_out + 0.5, param->fsR, param->meR));
    /* (partial xHII)/(partial lna). Use xH1s = 1-xHII for better accuracy. */
    if (xHeII_out != 0.)
    {
        Dxe = 0.01 * (param->fHe - xHeII_out);
        dxH1sSaha_dlna += (rec_saha_xH1s(xHeII_out + Dxe, param->nH0, param->T0, z_out, param->fsR, param->meR) - rec_saha_xH1s(xHeII_out - Dxe, param->nH0, param->T0, z_out, param->fsR, param->meR)) / 2. / Dxe * rec_helium_dxHeIIdlna(xH1sSaha, xHeII_out, param->nH0, param->T0, param->fHe, H, z_out, param->fsR, param->meR);
        /* (partial xHII)/(partial xHeII).dxHeII/dlna */
    }
    dxH1sdlna_Saha = -rec_dxHIIdlna(MODEL, xHIISaha + xHeII_out, xHIISaha, nH, H, T, T, rate_table, twog_params,
                                    Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz_out - param->izH0, z_out, param->fsR, param->meR, DarkArray);
    Dxe = 0.01 * xH1sSaha;
    DdxH1sdlna_DxH1s = (rec_dxHIIdlna(MODEL, xHIISaha + Dxe + xHeII_out, xHIISaha + Dxe, nH, H, T, T, rate_table, twog_params,
                                      Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz_out - param->izH0, z_out, param->fsR, param->meR, DarkArray) -
                        rec_dxHIIdlna(MODEL, xHIISaha - Dxe + xHeII_out, xHIISaha - Dxe, nH, H, T, T, rate_table, twog_params,
                                      Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz_out - param->izH0, z_out, param->fsR, param->meR, DarkArray)) /
                       2. / Dxe;

    xH1s = xH1sSaha + (dxH1sSaha_dlna - dxH1sdlna_Saha) / DdxH1sdlna_DxH1s;

    /* Check that we are still close enough to Saha equilibrium. If not, switch post-saha expansion off */
    if (fabs(xH1s - xH1sSaha) > DXHII_MAX)
        *post_saha = 0;

    return xH1s;
}

/*****************************************************************************************************
Second-order integrator used to evolve simultaneously H(1s) and HeII.
When hydrogen is close to Saha equilibrium but there is still a significant amount of HeII,
use a post-Saha expansion for hydrogen. The other way around is not treated (simply integrate H and He until
there is almost no HeII left, then integrate H only)
******************************************************************************************************/

void get_rec_next2_HHe(REC_COSMOPARAMS *param, unsigned iz_in, double z_in, double Tm_in, double *xH1s, double *xHeII,
                       HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist,
                       double *dxHIIdlna_prev, double *dxHeIIdlna_prev, double *dxHIIdlna_prev2, double *dxHeIIdlna_prev2, int *post_saha, double *DarkArray)
{

    double H, dxHeIIdlna, dxHIIdlna, TR, TM, nH, zout, ainv, xH1s_in, xHeII_in, xe_in;

    xH1s_in = *xH1s;
    xHeII_in = *xHeII;
    xe_in = xHeII_in + (1. - xH1s_in);

    /* Evolve HeII by solving ODE */
    H = rec_HubbleConstant(param, z_in);
    dxHeIIdlna = rec_helium_dxHeIIdlna(xH1s_in, xHeII_in, param->nH0, param->T0, param->fHe, H, z_in, param->fsR, param->meR);
    *xHeII += DLNA * (1.25 * dxHeIIdlna - 0.25 * (*dxHeIIdlna_prev2));

    /* Compute Hydrogen derivative at input time. Even if using the post-Saha expansion, needed for getting the correct radiation field at z_in */
    TR = kBoltz * param->T0 * (ainv = 1. + z_in);
    TM = kBoltz * Tm_in;
    nH = 1e-6 * param->nH0 * ainv * ainv * ainv;
    dxHIIdlna = rec_dxHIIdlna(MODEL, xe_in, (1. - xH1s_in), nH, H, TM, TR, rate_table, twog_params,
                              Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz_in - param->izH0, z_in, param->fsR, param->meR, DarkArray);

    /* If Hydrogen is still close to Saha equilibrium do a post-Saha expansion for Hydrogen */
    if (*post_saha == 1)
    {
        zout = (1. + z_in) * exp(-DLNA) - 1.; /* Redshift for the output */
        *xH1s = rec_xH1s_postSaha(param, iz_in + 1, zout, *xHeII, rate_table, twog_params,
                                  Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, post_saha, DarkArray);
    }

    /* Otherwise solve HII ODE */
    else
        *xH1s -= DLNA * (1.25 * dxHIIdlna - 0.25 * (*dxHIIdlna_prev2));

    /* Update derivatives */
    *dxHIIdlna_prev2 = *dxHIIdlna_prev;
    *dxHIIdlna_prev = dxHIIdlna;
    *dxHeIIdlna_prev2 = *dxHeIIdlna_prev;
    *dxHeIIdlna_prev = dxHeIIdlna;
}

/*********************************************************************************************************
Second-order integrator to evolve hydrogen only, assuming helium has entirely recombined.
Tm is given as an input (to avoid computing it twice) and fixed to quasi-equilibrium value with Tr.
**********************************************************************************************************/

void rec_get_xe_next1_H(REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out,
                        HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params,
                        double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist,
                        double *dxedlna_prev, double *dxedlna_prev2, int *post_saha, double *DarkArray)
{

    double dxedlna, TR, nH, ainv, H, TM, zout;
    int model;

    TR = kBoltz * param->T0 * (ainv = 1. + z_in);
    nH = 1e-6 * param->nH0 * ainv * ainv * ainv;
    H = rec_HubbleConstant(param, z_in);
    TM = kBoltz * Tm_in;

    /* Switch off radiative transfer calculation if needed (param->nzrt and izH0 are evaluated in rec_get_cosmoparam) */
    model = (iz - param->izH0 < param->nzrt || MODEL != FULL) ? MODEL : EMLA2s2p;

    dxedlna = rec_dxHIIdlna(model, xe_in, xe_in, nH, H, TM, TR, rate_table, twog_params,
                            Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz - param->izH0, z_in, param->fsR, param->meR, DarkArray);

    /* If close to Saha equilibrium (with xHeII = 0), do a post-Saha expansion */
    if (*post_saha == 1)
    {
        zout = (1. + z_in) * exp(-DLNA) - 1.; /* Redshift for the output */
        *xe_out = 1. - rec_xH1s_postSaha(param, iz + 1, zout, 0., rate_table, twog_params,
                                         Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, post_saha, DarkArray);
    }

    /* Otherwise evolve ODE */
    else
        *xe_out = xe_in + DLNA * (1.25 * dxedlna - 0.25 * (*dxedlna_prev2));

    /* Update previous derivatives */
    *dxedlna_prev2 = *dxedlna_prev;
    *dxedlna_prev = dxedlna;
}

/**********************************************************************************************
Second-order integrator for evolving xe and Tm simultaneously
Used for Hydrogen recombination only
May 2012: added a switch so Peebles model can be used at low redshift.
***********************************************************************************************/

void rec_get_xe_next2_HTm(int func_select, REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                          HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[],
                          double **Dfnu_hist, double *dxedlna_prev, double *dTmdlna_prev, double *dxedlna_prev2, double *dTmdlna_prev2, double *DarkArray)
{

    double dxedlna, dTmdlna, TR, nH, ainv, H, TM;
    int model;

    TR = kBoltz * param->T0 * (ainv = 1. + z_in);
    nH = 1e-6 * param->nH0 * ainv * ainv * ainv;
    H = rec_HubbleConstant(param, z_in);
    TM = kBoltz * Tm_in;

    /* Switch off radiative transfer calculation if needed */
    model = (iz - param->izH0 < param->nzrt || func_select != FULL) ? func_select : EMLA2s2p;

    dxedlna = rec_dxHIIdlna(model, xe_in, xe_in, nH, H, TM, TR, rate_table, twog_params,
                            Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz - param->izH0, z_in, param->fsR, param->meR, DarkArray);

    dTmdlna = rec_dTmdlna(xe_in, Tm_in, TR / kBoltz, H, param->fHe, param->fsR, param->meR, DarkArray);

    *xe_out = xe_in + DLNA * (1.25 * dxedlna - 0.25 * (*dxedlna_prev2));
    *Tm_out = Tm_in + DLNA * (1.25 * dTmdlna - 0.25 * (*dTmdlna_prev2));

    *dxedlna_prev2 = *dxedlna_prev;
    *dTmdlna_prev2 = *dTmdlna_prev;
    *dxedlna_prev = dxedlna;
    *dTmdlna_prev = dTmdlna;
}

/****************************************************************************************************
Builds a recombination history
Added May 2012: The radiation field was added as an input so it can be extracted if desired
****************************************************************************************************/

void rec_build_history(REC_COSMOPARAMS *param, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                       double *xe_output, double *Tm_output, double **Dfnu_hist, double *Dfminus_Ly_hist[3])
{
    // jsc

    long iz;
    double z, dxHIIdlna_prev, dxHIIdlna_prev2, dTmdlna_prev, dTmdlna_prev2, dxHeIIdlna_prev, dxHeIIdlna_prev2;
    double Delta_xe, xHeII, xH1s, f1, f3, f4;
    double **Dfminus_hist;
    int post_saha;
    // Dark Matter related quantities, MUST be of the form:
    // DarkArray = [dEdVdt_HIon, dEdVdt_LyA, dEdVdt_Heat, .... ]
    Validate_Inputs(param);
    double DarkArray[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};

    Dfminus_hist = create_2D_array(NVIRT, param->nzrt);

    /* Make sure the input spectrum is initialized at zero */
    for (iz = 0; iz < param->nzrt; iz++)
        Dfminus_Ly_hist[0][iz] = Dfminus_Ly_hist[1][iz] = Dfminus_Ly_hist[2][iz] = 0;

    z = ZSTART;

    /********* He III -> II Saha phase. Tm = Tr. Stop when xHeIII = 1e-8 *********/
    Delta_xe = param->fHe; /* Delta_xe = xHeIII here */

    for (iz = 0; iz < param->nz && Delta_xe > 1e-8; iz++)
    {
        Update_DarkArray(z, param, DarkArray);
        if (param->Use_SSCK == 1)
        {
            if (iz > 0)
            {
                f1 = SSCK_EFF(xe_output[iz - 1], 1);
                f3 = SSCK_EFF(xe_output[iz - 1], 3);
                f4 = SSCK_EFF(xe_output[iz - 1], 4);
            }
            else
            {
                f1 = 0.0;
                f3 = 0.0;
                f4 = 0.0;
            }
            DarkArray[0] = DarkArray[2] * f1;
            DarkArray[1] = DarkArray[2] * f3;
            DarkArray[2] = DarkArray[2] * f4;
        }

        z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
        xe_output[iz] = rec_xesaha_HeII_III(param->nH0, param->T0, param->fHe, z, &Delta_xe, param->fsR, param->meR);
        Tm_output[iz] = param->T0 * (1. + z);
        Check_Error(xe_output, Tm_output, z, iz, param);
        if (debug_mode)
        {
            printf("Stage_1: z = %f, xe = %f, Tm = %f\n", z, xe_output[iz], Tm_output[iz]);
        }
    }

    /******** He II -> I recombination.
              Hydrogen in Saha equilibrium with the free electrons.
              Tm fixed to steady state.
              Integrate until TR is low enough that can start integrating hydrogen recombination
              (this occurs at index izH0 computed in rec_get_cosmoparam).
              Start with post-Saha expansion.
     ********/

    dxHeIIdlna_prev2 = (xe_output[iz - 2] - xe_output[iz - 4]) / 2. / DLNA;
    dxHeIIdlna_prev = (xe_output[iz - 1] - xe_output[iz - 3]) / 2. / DLNA;

    xHeII = rec_saha_xHeII(param->nH0, param->T0, param->fHe, z, param->fsR, param->meR);
    post_saha = 1; /* Start with post-saha expansion */

    for (; iz < param->izH0 + 1; iz++)
    {
        Update_DarkArray(z, param, DarkArray);
        if (param->Use_SSCK == 1)
        {
            f1 = SSCK_EFF(xe_output[iz - 1], 1);
            f3 = SSCK_EFF(xe_output[iz - 1], 3);
            f4 = SSCK_EFF(xe_output[iz - 1], 4);
            DarkArray[0] = DarkArray[2] * f1;
            DarkArray[1] = DarkArray[2] * f3;
            DarkArray[2] = DarkArray[2] * f4;
        }

        rec_get_xe_next1_He(param, z, &xHeII, &dxHeIIdlna_prev, &dxHeIIdlna_prev2, &post_saha);
        z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
        xH1s = rec_saha_xH1s(xHeII, param->nH0, param->T0, z, param->fsR, param->meR);
        xe_output[iz] = (1. - xH1s) + xHeII;
        Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0 * (1. + z), rec_HubbleConstant(param, z), param->fHe, param->fsR, param->meR, DarkArray);
        Check_Error(xe_output, Tm_output, z, iz, param);
        if (debug_mode)
        {
            printf("Stage_2: z = %f, xe = %f, Tm = %f\n", z, xe_output[iz], Tm_output[iz]);
        }
    }

    /******** H II -> I and He II -> I simultaneous recombination (rarely needed but just in case)
              Tm fixed to steady state.
              Integrate H and He simultaneously until xHeII < XHEII_MIN
              Start with post-saha expansion for hydrogen
     ********/

    dxHIIdlna_prev2 = (xe_output[iz - 2] - xe_output[iz - 4]) / 2. / DLNA - dxHeIIdlna_prev2;
    dxHIIdlna_prev = (xe_output[iz - 1] - xe_output[iz - 3]) / 2. / DLNA - dxHeIIdlna_prev;
    post_saha = 1;

    for (; iz < param->nz && xHeII > XHEII_MIN; iz++)
    {
        Update_DarkArray(z, param, DarkArray);
        if (param->Use_SSCK == 1)
        {
            f1 = SSCK_EFF(xe_output[iz - 1], 1);
            f3 = SSCK_EFF(xe_output[iz - 1], 3);
            f4 = SSCK_EFF(xe_output[iz - 1], 4);
            DarkArray[0] = DarkArray[2] * f1;
            DarkArray[1] = DarkArray[2] * f3;
            DarkArray[2] = DarkArray[2] * f4;
        }

        get_rec_next2_HHe(param, iz - 1, z, Tm_output[iz - 1], &xH1s, &xHeII, rate_table, twog_params, Dfminus_hist, Dfminus_Ly_hist,
                          Dfnu_hist, &dxHIIdlna_prev, &dxHeIIdlna_prev, &dxHIIdlna_prev2, &dxHeIIdlna_prev2, &post_saha, DarkArray);
        xe_output[iz] = (1. - xH1s) + xHeII;
        z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
        Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0 * (1. + z), rec_HubbleConstant(param, z), param->fHe, param->fsR, param->meR, DarkArray);
        Check_Error(xe_output, Tm_output, z, iz, param);

        if (debug_mode)
        {
            printf("Stage_3: z = %f, xe = %f, Tm = %f\n", z, xe_output[iz], Tm_output[iz]);
        }
    }

    /******** H recombination. Helium assumed entirely neutral.
              Tm fixed to steady-state until its relative difference from Tr is DLNT_MAX
    ********/

    for (; iz < param->nz && 1. - Tm_output[iz - 1] / param->T0 / (1. + z) < DLNT_MAX; iz++)
    {
        Update_DarkArray(z, param, DarkArray);
        if (param->Use_SSCK == 1)
        {
            f1 = SSCK_EFF(xe_output[iz - 1], 1);
            f3 = SSCK_EFF(xe_output[iz - 1], 3);
            f4 = SSCK_EFF(xe_output[iz - 1], 4);
            DarkArray[0] = DarkArray[2] * f1;
            DarkArray[1] = DarkArray[2] * f3;
            DarkArray[2] = DarkArray[2] * f4;
        }

        rec_get_xe_next1_H(param, z, xe_output[iz - 1], Tm_output[iz - 1], xe_output + iz, rate_table, iz - 1, twog_params,
                           Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, &dxHIIdlna_prev, &dxHIIdlna_prev2, &post_saha, DarkArray);
        z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
        Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0 * (1. + z), rec_HubbleConstant(param, z), param->fHe, param->fsR, param->meR, DarkArray);
        Check_Error(xe_output, Tm_output, z, iz, param);

        if (debug_mode)
        {
            printf("Stage_4: z = %f, xe = %f, Tm = %f\n", z, xe_output[iz], Tm_output[iz]);
        }
    }

    /******** Evolve xe and Tm simultaneously until the lower bounds of integration tables are reached.
              Note that the radiative transfer calculation is switched off automatically in the functions
              rec_get_xe_next1_H and rec_get_xe_next2_HTm when it is no longer relevant.
    ********/

    dTmdlna_prev2 = (Tm_output[iz - 2] - Tm_output[iz - 4]) / 2. / DLNA;
    dTmdlna_prev = (Tm_output[iz - 1] - Tm_output[iz - 3]) / 2. / DLNA;

    for (; iz < param->nz && kBoltz * param->T0 * (1. + z) / param->fsR / param->fsR / param->meR > TR_MIN && Tm_output[iz - 1] / param->T0 / (1. + z) > TM_TR_MIN; iz++)
    {
        Update_DarkArray(z, param, DarkArray);
        if (param->Use_SSCK == 1)
        {
            f1 = SSCK_EFF(xe_output[iz - 1], 1);
            f3 = SSCK_EFF(xe_output[iz - 1], 3);
            f4 = SSCK_EFF(xe_output[iz - 1], 4);
            DarkArray[0] = DarkArray[2] * f1;
            DarkArray[1] = DarkArray[2] * f3;
            DarkArray[2] = DarkArray[2] * f4;
        }

        rec_get_xe_next2_HTm(MODEL, param, z, xe_output[iz - 1], Tm_output[iz - 1], xe_output + iz, Tm_output + iz,
                             rate_table, iz - 1, twog_params, Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist,
                             &dxHIIdlna_prev, &dTmdlna_prev, &dxHIIdlna_prev2, &dTmdlna_prev2, DarkArray);
        z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
        Check_Error(xe_output, Tm_output, z, iz, param);

        if (debug_mode)
        {
            printf("Stage_5: z = %f, xe = %f, Tm = %f\n", z, xe_output[iz], Tm_output[iz]);
        }
    }

    /***** For low redshifts (z < 20 or so) use Peeble's model (Tm is evolved with xe).
            The precise model does not metter much here as
            1) the free electron fraction is basically zero (~1e-4) in any case and
            2) the universe is going to be reionized around that epoch
     *****/

    for (; iz < param->nz; iz++)
    {
        Update_DarkArray(z, param, DarkArray);
        if (param->Use_SSCK == 1)
        {
            f1 = SSCK_EFF(xe_output[iz - 1], 1);
            f3 = SSCK_EFF(xe_output[iz - 1], 3);
            f4 = SSCK_EFF(xe_output[iz - 1], 4);
            DarkArray[0] = DarkArray[2] * f1;
            DarkArray[1] = DarkArray[2] * f3;
            DarkArray[2] = DarkArray[2] * f4;
        }

        rec_get_xe_next2_HTm(PEEBLES, param, z, xe_output[iz - 1], Tm_output[iz - 1], xe_output + iz, Tm_output + iz,
                             rate_table, iz - 1, twog_params, Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist,
                             &dxHIIdlna_prev, &dTmdlna_prev, &dxHIIdlna_prev2, &dTmdlna_prev2, DarkArray);
        z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
        Check_Error(xe_output, Tm_output, z, iz, param);

        if (debug_mode)
        {
            printf("Stage_6: z = %f, xe = %f, Tm = %f\n", z, xe_output[iz], Tm_output[iz]);
        }
    }

    /* Cleanup */
    free_2D_array(Dfminus_hist, NVIRT);
}
/******************************************************************************************************************************/
