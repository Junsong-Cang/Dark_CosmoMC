/* Dark Matter Energy injection module
main reference : 2108.13256
*/

// Rho_cr * c^2 /h^2, in eV/cm^3
#define Rhocr_C2_no_h2 1.054041699E4
#define Nearly_Zero 1.0E-50
#define Zmax_Energy_Injection 2000

/* Set number of timesteps in integrating PBH deposition rate
100 takes 8 seconds and should be good enough
*/
#define PBH_Integration_Size 100
#define debug_mode 0
// Only include this after definning Rhocr_C2_no_h2 otherwise there will be redefinition error
#include "Accreting_PBH.h"

void Validate_Inputs(REC_COSMOPARAMS *params)
{
    /* Check input params
    1. Particle Channel
    2. PBH Spin
    3. PBH Model: Currently only support [1 2 3]
    */
    int Particle_Channel, PBH_Model, PBH_Distribution;
    Particle_Channel = Convert_to_Int(params->DM_Channel);
    PBH_Model = Convert_to_Int(params->PBH_Model);
    PBH_Distribution = Convert_to_Int(params->PBH_Distribution);

    if ((Particle_Channel < 1) || (Particle_Channel > 12))
    {
        printf("Error from Validate_Inputs@HyRec: Unknown particle channel selected\n");
        exit(1);
    }
    if ((params->PBH_Spin > 0.99999) || (params->PBH_Spin < -1.0E-5))
    {
        printf("Error from Validate_Inputs@HyRec: wrong PBH_Spin.\n");
        exit(1);
    }
    if ((PBH_Model < 1) || (PBH_Model > 3))
    {
        printf("Error from Validate_Inputs@HyRec: Wrong choice of PBH_Model\n");
        exit(1);
    }
    if (fabs(params->PBH_PWL_Gamma) < 0.001)
    {
        printf("Error: PBH Power-Law index cannot be 0.\n");
        exit(1);
    }
    if (PBH_Distribution == 3)
    {
        if ((params->Mbh) > (params->PBH_PWL_Mmax))
        {
            printf("Error: Mmin > Mmax in power-law.\n");
            exit(1);
        }
    }
}

// ALL dEdVdt are in ev/cm^3/s unit

double dEdVdt_decay_inj(double z, REC_COSMOPARAMS *params)
{
    double Gamma, Omch2, r;
    Gamma = params->Gamma;
    Omch2 = params->odmh2;
    if (Gamma > 0)
    {
        r = Gamma * Omch2 * pow(1 + z, 3.0) * Rhocr_C2_no_h2;
    }
    else
    {
        r = 0.0;
    }
    return r;
}

double dEdVdt_ann_inj(double z, REC_COSMOPARAMS *params)
{
    double Omch2, r;
    Omch2 = params->odmh2;
    // 1E9 converts GeV to eV
    r = params->Pann * square(Rhocr_C2_no_h2 * Omch2 * cube(1. + z)) * 1.0E-9;
    return r;
}

double dEdVdt_Hawking_inj(double z, double Mbh, REC_COSMOPARAMS *params)
{
    // Hawking Radiation injection, Normalised to mbh>10^17 g
    // See Eq (3.10) of arxiv 2108.13256
    double r;
    r = (5.626976744186047e+29 / cube(Mbh) * params->fbh * params->odmh2 * cube(1. + z));
    return r;
}

double dEdVdt_Accretion_inj(double z, double Mbh, REC_COSMOPARAMS *params)
{
    /* Energy injeciton rate from accreting PBH
    ---- inputs ----
    Mbh : PBH mass in msun
    */
    double r, dEdt, nbh, rho_dm;
    rho_dm = params->odmh2 * 1.879E-26 * pow(1 + z, 3);
    nbh = params->fbh * rho_dm / (Mbh * Solar_Mass) / 1E6;
    dEdt = PBH_Accretion_Luminisoty(Mbh, z);
    r = dEdt * nbh;

    return r;
}

double EFF_SSCK(double z, int dep_channel)
{
    double xe, r;
    xe = LCDM_Xe(z);
    if (dep_channel == 1)
    {
        r = (1 - xe) / 3;
    }
    else if (dep_channel == 3)
    {
        r = (1 - xe) / 3;
    }
    else if (dep_channel == 4)
    {
        r = (1 + 2 * xe) / 3;
    }
    else
    {
        printf("Wrong choice of dep_channel.\n");
        exit(1);
    }
    // 1-xe might be <0
    r = fmax(r, 0);
    return r;
}

double dEdVdt_decay_dep(double z, REC_COSMOPARAMS *params, int dep_channel)
{
    double inj, r, EFF, Mdm;
    int DM_Channel;
    DM_Channel = Convert_to_Int(params->DM_Channel);
    inj = dEdVdt_decay_inj(z, params);
    EFF = Interp_EFF_DM_Decay(params->Mdm, z, dep_channel, DM_Channel);
    r = EFF * inj;
    return r;
}

double dEdVdt_Hawking_Mono_dep(double z, double Mbh, REC_COSMOPARAMS *params, int dep_channel)
{
    // Hawking Radiation monochromatic deosition rate
    double inj, r, EFF, Mdm;
    inj = dEdVdt_Hawking_inj(z, Mbh, params);
    EFF = Interp_EFF_Hawking(params->Mbh, z, params->PBH_Spin, params->PBH_Model, dep_channel);
    r = EFF * inj;
    return r;
}

double dEdVdt_Accretion_Mono_dep(double z, double Mbh, REC_COSMOPARAMS *params, int dep_channel)
{
    double inj, r, EFF;
    inj = dEdVdt_Accretion_inj(z, Mbh, params);
    EFF = EFF_SSCK(z, dep_channel);
    r = EFF * inj;
    return r;
}

double PBH_Mass_Function(double m, REC_COSMOPARAMS *params)
{
    // return PBH mass function, result has dimension of m^-1, unit of m is determined by PBH_Model
    double r, mc, sbh, mmin, mmax, y, pi, x;
    int PBH_Distribution;
    mc = params->Mbh;
    sbh = params->PBH_Lognormal_Sigma;
    mmin = mc;
    mmax = params->PBH_PWL_Mmax;
    y = params->PBH_PWL_Gamma;
    pi = 3.141592653589793;
    PBH_Distribution = Convert_to_Int(params->PBH_Distribution);

    if (PBH_Distribution == 2)
    {
        // log-normal
        r = pow(log(m / mc), 2);
        r = -r / (2 * sbh * sbh);
        r = exp(r) / (m * sbh * sqrt(2 * pi));
    }
    else if (PBH_Distribution == 3)
    {
        // Power-law
        if ((mmin <= m) && (m <= mmax))
        {
            r = y * pow(m, y - 1);
            r = r / (pow(mmax, y) - pow(mmin, y));
        }
        else
        {
            r = 0.0;
        }
    }
    else if (PBH_Distribution == 4)
    {
        x = pow(m / mc, 2.85);
        r = x * exp(-x);
        r = 3.2 * r / mc;
    }

    r = 0;
    return r;
}

double Integrate_PBH_dEdVdt(double z, REC_COSMOPARAMS *params, int dep_channel)
{
    // Find dEdVdt_dep for general PBH distributions, applicable to both Hawking Radiation and Accreting PBH
    int PBH_Distribution, PBH_Model, idx;
    double Mmin, Mmax, Tmin, Tmax, r, x1, x2, x1_, x2_, xc, sbh, dx, x, f1, f2, f, m;
    PBH_Distribution = Convert_to_Int(params->PBH_Distribution);
    PBH_Model = Convert_to_Int(params->PBH_Model);
    sbh = params->PBH_Lognormal_Sigma;

    if ((PBH_Model == 2) || (PBH_Model == 3))
    {
        Tmax = 0.9999 * Kerr_PBH_Temperature_Axis[Kerr_PBH_Temperature_Size - 1];
        Tmin = 1.0001 * Kerr_PBH_Temperature_Axis[0];
        Mmin = 1.06E22 / Tmax;
        Mmax = 1.06E22 / Tmin;
    }
    else
    {
        printf("Error: Wrong choice of PBH_Model.\n");
        exit(1);
    }

    // Doing integration in log space, x defined as log(m)
    x1 = log(Mmin);
    x2 = log(Mmax);
    // Optimise integration boundary
    if (PBH_Distribution == 2)
    {
        // log-normal: +- 5 sigma
        xc = log(params->Mbh);
        x1_ = xc - 5 * sbh;
        x2_ = xc + 5 * sbh;
    }
    else if (PBH_Distribution == 3)
    {
        // power-law: [Mmin, Mmax]
        x1_ = log(params->Mbh);
        x2_ = log(params->PBH_PWL_Mmax);
    }
    // updating
    x1 = fmax(x1, x1_);
    x2 = fmin(x2, x2_);

    dx = (x2 - x1) / ((double)PBH_Integration_Size - 1.0);
    r = 0.0;
    x = x1;

    for (idx = 0; idx < PBH_Integration_Size; idx++)
    {
        m = exp(x);
        if ((PBH_Model == 2) || (PBH_Model == 3))
        {
            f1 = dEdVdt_Hawking_Mono_dep(z, m, params, dep_channel);
        }
        else
        {
            printf("Accreting PBH not ready.\n");
            exit(1);
        }
        f2 = PBH_Mass_Function(m, params);
        f = f1 * f2 * m;
        r = r + f * dx;
        x = x + dx;
    }

    return r;
}

double dEdVdt_Hawking_dep(double z, REC_COSMOPARAMS *params, int dep_channel)
{
    // Hawking Radiation for general mass distributions
    int PBH_Distribution;
    double r;

    PBH_Distribution = Convert_to_Int(params->PBH_Distribution);
    if (PBH_Distribution == 1)
    {
        r = dEdVdt_Hawking_Mono_dep(z, params->Mbh, params, dep_channel);
    }
    else
    {
        r = Integrate_PBH_dEdVdt(z, params, dep_channel);
    }

    return r;
}

double dEdVdt_ann_dep(double z, REC_COSMOPARAMS *params, int dep_channel)
{
    double inj, r, EFF, Mdm;
    int DM_Channel;
    DM_Channel = (int)round(params->DM_Channel);
    inj = dEdVdt_ann_inj(z, params);
    EFF = Interp_EFF_DM_Annihilation(params->Mdm, z, dep_channel, DM_Channel);
    r = EFF * inj;
    return r;
}

double dEdVdt_deposited(double z, REC_COSMOPARAMS *params, int dep_channel)
{
    /* Energy Deposition Rate in ev/cm^3/s
     -- inputs --
     dep_channel = 1: HIon
                   3: LyA
                   4: Heating
    */

    double r_dec, r_ann, r_bh, r;
    int PBH_Model, PBH_Distribution;
    PBH_Model = Convert_to_Int(params->PBH_Model);
    PBH_Distribution = Convert_to_Int(params->PBH_Distribution);

    // Check params
    if ((dep_channel == 2) || (dep_channel == 5))
    {
        printf("HeIon and Continnum dep channels not allowed, exitting\n");
        exit(1);
    }
    if (z > Zmax_Energy_Injection)
    {
        return 0.0;
    }
    // decay
    if (params->Gamma < Nearly_Zero)
    {
        r_dec = 0.0;
    }
    else
    {
        r_dec = dEdVdt_decay_dep(z, params, dep_channel);
    }
    // annihilation
    if (params->Pann < Nearly_Zero)
    {
        r_ann = 0.0;
    }
    else
    {
        r_ann = dEdVdt_ann_dep(z, params, dep_channel);
    }
    // pbh
    if (params->fbh < Nearly_Zero)
    {
        r_bh = 0.0;
    }
    else
    {
        if (PBH_Model == 1)
        {
            if (PBH_Distribution != 1)
            {
                printf("Module not ready for accreting PBH with extended distribution");
                exit(1);
            }
            else
            {
                r_bh = dEdVdt_Accretion_Mono_dep(z, params->Mbh, params, dep_channel);
            }
        }
        else if (PBH_Model == 2)
        {
            r_bh = dEdVdt_Hawking_dep(z, params, dep_channel);
        }
        else
        {
            printf("Unknown PBH Model.\n");
            exit(1);
        }
    }

    r = r_dec + r_ann + r_bh;

    return r;
}

void Update_DarkArray(double z, REC_COSMOPARAMS *params, double *DarkArray)
{
    double dEdVdt_HIon, dEdVdt_LyA, dEdVdt_Heat, nH;
    dEdVdt_HIon = dEdVdt_deposited(z, params, 1);
    dEdVdt_LyA = dEdVdt_deposited(z, params, 3);
    dEdVdt_Heat = dEdVdt_deposited(z, params, 4);
    // Let's fill dEdVdt first and see what we can do with it
    DarkArray[0] = dEdVdt_HIon;
    DarkArray[1] = dEdVdt_LyA;
    DarkArray[2] = dEdVdt_Heat;
    // H number density in cm^3
    nH = params->nH0 * cube(1 + z) * 1.E-6;
    DarkArray[3] = nH;
}

void Check_Error(double xe, double T)
{
    // Check for inifnity and NaN in Xe and T
    if (isfinite(xe) == 0)
    {
        printf("Error from Check_Error@HyRec: xe is NaN or infinite\n");
        exit(1);
    }
    if (isfinite(T) == 0)
    {
        printf("Error from Check_Error@HyRec: T is NaN or infinite\n");
        exit(1);
    }
}