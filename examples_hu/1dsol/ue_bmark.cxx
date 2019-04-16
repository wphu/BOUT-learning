/*******************************************************************************
 * UEDGE benchmark case
 *
 * Solves equations for
 *  density Ni
 *  parallel ion velocity Vi
 *  electron and ion temperatures Te, Ti
 *
 * Intended to be run for NZ=1 (i.e. X and Y only) for comparison with UEDGE
 *
 *******************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <derivs.hxx>

#include <cmath>


// 2D initial profiles
//Field2D Ni0;

// 3D evolving fields
Field2D  Ni;
Field2D  Vi;
Field2D  Pii;
Field2D  Ti;
Field2D  Te;
Field2D  E;




//3D intermidia fields for boundary appling


// parameters
BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;


// max number of grid file nx
int nx,ny;

//logical paramter
bool Test_Ni_Eq;

double nt  = 1.0e19;
double cst = 600000;
double Liz   = 4.0;

double Sp = nt * cst / Liz;
double e  = 1.602e-19;
double K  = 1.3806e-23;
double mi = 2.0 * 1.67262158e-27;
double Qeq = 0.0;
double QEi = 0.0;
double QEe = 0.0;
double QR  = 0.0;
double k0e = 2000.0;
double k0i = 60.0;
double t_pv= 1.0e-6;

int physics_init(bool restarting)
{
    output.write("Solving transport equations for Ni, Vi, Ti, Te\n");

    /////////////// LOAD DATA FROM GRID FILE //////////////

    mesh->get(mesh->dx,   "dx");
    mesh->get(nx,"NX");
    // Load normalisation values
    GRID_LOAD(Te_x);
    GRID_LOAD(bmag);
    output.write("Matrics coefficient : rho_snx %d  \n",nx);

    bmag *= 1.0e4;

    /////////////// READ OPTIONS //////////////////////////

    // Read some parameters
    Options *globalOptions = Options::getRoot();
    Options *options = globalOptions->getSection("uedge");
    OPTION(options, AA, 2.0);
    OPTION(options, ZZ, 1.0);



    ////////////// CALCULATE PARAMETERS ///////////////////

    output.write("Diffusion coefficients: rho_s %e AA %e Te_x %e bmag %e ZZ %e \n",
                 rho_s,AA,Te_x,bmag,ZZ);

    rho_s = 1.02*sqrt(AA*Te_x)/ZZ/bmag;


    ///////////// NORMALISE QUANTITIES ////////////////////

    output.write("\tNormalising to rho_s = %e\n", rho_s);
    mesh->dx = 0.01; ///= rho_s;

    // Normalise magnetic field
    output.write("Diffusion coefficients: rho_s %e AA %e Te_x %e bmag %e ZZ %e \n",
                  rho_s,AA,Te_x,bmag,ZZ);


    Ni.setBoundary("Ni");
    Vi.setBoundary("Vi");
    Pii.setBoundary("Pii");
    Ti.setBoundary("Ti");
    Te.setBoundary("Te");
    E.setBoundary("E");

    SOLVE_FOR(Ni);
    SOLVE_FOR(Vi);
    SOLVE_FOR(Pii);
    SOLVE_FOR(Ti);
    SOLVE_FOR(Te);
    SOLVE_FOR(E);

    for (int jx=0;jx<mesh->ngx;jx++)
    {
        for (int jy=0;jy<mesh->ngy;jy++)
        {
            Ni[jx][jy]=1.0e17;
            Vi[jx][jy]=1.0;
            Pii[jx][jy]=1.0;
            Ti[jx][jy]=100.0;
            Te[jx][jy]=100.0;
            E[jx][jy]=10.0;
        }
    }

    return(0);
}




int physics_run(BoutReal t)
{
    Ni.applyBoundary();
    Vi.applyBoundary();
    Pii.applyBoundary();
    Ti.applyBoundary();
    Te.applyBoundary();
    E.applyBoundary();

    ddt(Ni)     = -DDX(Ni * Vi) + Sp;

    ddt(Vi)     = -DDX(mi * Ni * Vi * Vi + Ni * K * Ti + Ni * K * Te + Pii) / (mi * Ni)
                  - 0.0;

    ddt(Pii)    = -Pii - (4.0/9.0) * Ni * K * Ti * t_pv * DDX(Vi);

    ddt(Ti)     = -DDX(2.5 * Ni * K * Ti * Vi + 0.5 * mi * Ni * Vi * Vi * Vi + Pii * Vi)
                  + k0i * Ti^2.5 * D2DX2(Ti) + e * Ni * Vi * E + Qeq + QEi;

    ddt(Te)     = -DDX(2.5 * Ni * K * Te * Vi) + k0e * Te^2.5 * D2DX2(Te)
                  - e * Ni * Vi * E - Qeq + QR + QEe;

    ddt(E)      = -E - 0.71 * K * DDX(Te) / e - DDX(Ni * K * Te) / (e * Ni);


    return(0);
}
