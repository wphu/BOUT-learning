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





//3D intermidia fields for boundary appling


// parameters
BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;


// max number of grid file nx
int nx,ny;

//logical paramter
bool Test_Ni_Eq;

int physics_init(bool restarting)
{
  output.write("Solving transport equations for Ni, Vi, Ti, Te\n");

  /////////////// LOAD DATA FROM GRID FILE //////////////

  mesh->get(mesh->dx,   "dx");

  mesh->get(nx,"NX");

 output.write("Matrics coefficient : rho_snx %d  \n",nx);

  // Load normalisation values
  GRID_LOAD(Te_x);

  GRID_LOAD(bmag);


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


;  mesh->dx /= rho_s;

  // Normalise magnetic field

  
  output.write("Diffusion coefficients: rho_s %e AA %e Te_x %e bmag %e ZZ %e \n",
	       rho_s,AA,Te_x,bmag,ZZ);

  /////////////// CALCULATE METRICS /////////////////

 
   Ni.setBoundary("Ni");

  SOLVE_FOR(Ni);

  for (int jx=0;jx<mesh->ngx;jx++)
     {
      for (int jy=0;jy<mesh->ngy;jy++)
	{
// 	for (int jz=0;jz<mesh->ngz;jz++)
//          {
//	   Ni[jx][jy][jz]=(BoutReal)jx/(BoutReal)(mesh->ngx-1);
//	   Ni[jx][jy]=(BoutReal)jx*(BoutReal)jy/(BoutReal)(mesh->ngx-1);
	   Ni[jx][jy]=(BoutReal)jx/(BoutReal)(mesh->ngx-1);
//           }
        }
     }



  return(0);
}




int physics_run(BoutReal t)
{


  Ni.applyBoundary();
 

  // ddt(Ni) = 1000.*D2DX2(Nit)+1.;
  ddt(Ni) =1.*D2DX2(Ni);



  return(0);
}























