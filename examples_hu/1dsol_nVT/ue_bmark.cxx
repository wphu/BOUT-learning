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
#include <field_factory.hxx>

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


Field2D  Sp;



//3D intermidia fields for boundary appling


// parameters
BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;


// max number of grid file nx
int nx,ny;

//logical paramter
bool Test_Ni_Eq;

double nt  = 1.0e16;
double cst = 600000;
double Liz   = 4.0;

double e  = 1.602e-19;
double K  = 1.3806e-23;
double mi = 2.0 * 1.67262158e-27;
double Qeq = 0.0;
double QEi = 0.0;
double QEe = 0.0;
double QR  = 0.0;
double k0e = 2000.0;
double k0i = 600.0;
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
    mesh->dx = 0.1; ///= rho_s;

    // Normalise magnetic field
    output.write("Diffusion coefficients: rho_s %e AA %e Te_x %e bmag %e ZZ %e \n",
                  rho_s,AA,Te_x,bmag,ZZ);


    Ni.setBoundary("Ni");
    Vi.setBoundary("Vi");
    Ti.setBoundary("Ti");
    Sp.setBoundary("Sp");

    SOLVE_FOR(Ni);
    SOLVE_FOR(Vi);
    SOLVE_FOR(Ti);
    SOLVE_FOR(Sp);

    int ls = mesh->ngx * 0.45;
    int length = mesh->ngx * 0.5;

    for (int jx=0;jx<mesh->ngx;jx++)
    {
        for (int jy=0;jy<mesh->ngy;jy++)
        {
            Ni[jx][jy]=1.0e17;
            Vi[jx][jy]=0.0;
            Ti[jx][jy]=1000.0;

	    if(jx > length - ls && jx < length + ls)
	    {
		Sp[jx][jy] = nt * cst / Liz;
	    }
	    else
	    {
	        Sp[jx][jy] = 0.0;
	    }
        }
    }

    return(0);
}




int physics_run(BoutReal t)
{

    mesh->communicate(Ni);
    mesh->communicate(Vi);
    mesh->communicate(Ti);

    ddt(Ni)     = -Ni * DDX(Vi) - VDDX(Vi, Ni) + Sp;

    ddt(Vi)     = -2.0 * VDDX(Vi, mi * Vi) / (mi) -  VDDX(Vi * Vi, mi * Ni) / (mi * Ni)
		  - DDX(2.0 *  K * Ti) / (mi) - Ti * DDX(2.0 * Ni * K) / (mi * Ni);


    ddt(Ti)     = -VDDX(Vi, 5.0 * Ni * K * Ti) - 5.0 * Ni * K * Ti * DDX(Vi)
		  - VDDX(Vi * Vi * Vi, 0.5 * mi * Ni) - 0.5 * mi * Ni * DDX(Vi * Vi * Vi)
                  + k0i * Ti^2.5 * D2DX2(Ti) + Qeq + QEi;

    ddt(Sp)	= 0.0;

    return(0);
}
