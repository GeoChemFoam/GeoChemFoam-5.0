/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the GeoChemFoam Group at Heriot-Watt University. Please visit our
    website for more information <https://github.com/GeoChemFoam>.

Application
    processheatTransfer 

Description
    calculate heat transfer coefficient between solid and fluid

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "fvCFD.H"
#include "argList.H"
#include "primitivePatchInterpolation.H"
#include "timeSelector.H"


using namespace Foam;


int main(int argc, char *argv[])
{
	#   include "setRootCase.H"
	#   include "createTime.H"


	instantList timeList = timeSelector::select0(runTime, args);

        #   include "createNamedMesh.H"

        volScalarField eps0
        (
            IOobject
            (
                "eps",
                 runTime.timeName(),
                 mesh,
                 IOobject::READ_IF_PRESENT,
                 IOobject::NO_WRITE
             ),
             mesh,
             dimensionedScalar("esp",dimless,1.0)
        );

        volScalarField tinv0
        (
            IOobject
            (
                "tinv",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
             ),
             mesh,
             dimensionedScalar("tinv",dimless,1.0)
        );



        string filename = "ST.csv";
	std::ofstream csvfile(filename.c_str());
	csvfile << "time Tf Tm coeff\n";

        
       forAll(timeList, timeStep)
       {	
	    	runTime.setTime(timeList[timeStep], timeStep);

		Info<< endl<<timeStep<< "    Time = " << runTime.timeName() << endl;

                #   include "createNamedMesh.H"

                IOdictionary postProcessDict
                (
                    IOobject
                    (
                        "postProcessDict",
                        "system",
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    )
                );

                scalar z1
                (
                    postProcessDict.lookupOrDefault<scalar>("z1",-1e30)
                );


                scalar z2
                (
                    postProcessDict.lookupOrDefault<scalar>("z2",1e30)
                );

                scalar y1
                (
                    postProcessDict.lookupOrDefault<scalar>("y1",-1e30)
                );


                scalar y2
                (
                    postProcessDict.lookupOrDefault<scalar>("y2",1e30)
                );


                scalar x1
                (
                    postProcessDict.lookupOrDefault<scalar>("x1",1e-30)
                );


                scalar x2
                (
                    postProcessDict.lookupOrDefault<scalar>("x2",1e30)
                );


                volScalarField clip
                (
                        IOobject
                        (
                            "clip",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh,
                        dimensionedScalar("clip",dimVolume,0.0)
                );

                volScalarField coordz=mesh.C().component(2);
                volScalarField coordy=mesh.C().component(1);
                volScalarField coordx=mesh.C().component(0);


               forAll(mesh.cells(),j)
               {
                    scalar zj=coordz[j];
                    scalar yj=coordy[j];
                    scalar xj=coordx[j];
                    if (zj>=z1 && zj<=z2 && yj>=y1 && yj<=y2 && xj>=x1 && xj<=x2)
                    {
                        clip[j]=mesh.V()[j];
                    }
                }

               Info<< "Reading transportProperties\n" << endl;

               IOdictionary transportProperties
               (
                   IOobject
                   (
                       "transportProperties",
                       runTime.constant(),
                       mesh,
                       IOobject::MUST_READ,
                       IOobject::NO_WRITE
                   )
               );


	       dimensionedScalar DT("DT", dimViscosity, transportProperties);

		volScalarField eps 
		(
			IOobject
			(
				"eps",
				runTime.timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			),
			eps0
		);

		volScalarField eps2 = 1-eps;

                volScalarField tinv
                (
                    IOobject
                    (
                        "tinv",
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    tinv0
                );


                volScalarField T
                (
                        IOobject
                        (
                                "T",
                                runTime.timeName(),
                                mesh,
                                IOobject::MUST_READ,
                                IOobject::NO_WRITE
                        ),
                        mesh
                );


                scalar Tf = T.weightedAverage(eps*clip).value();

                scalar Tm = T.weightedAverage(eps2*clip).value();

                surfaceScalarField D  = fvc::interpolate(eps*tinv*DT);;
            
                surfaceScalarField phiD = D*fvc::snGrad(T)*mesh.magSf();

                surfaceScalarField epsf = fvc::interpolate(eps);
 
                volScalarField F = fvc::div(phiD*epsf)-eps*fvc::div(phiD);

                scalar ST = F.weightedAverage(clip).value();

                volScalarField gEps = mag(fvc::grad(eps));

                scalar A = gEps.weightedAverage(clip).value();

		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << Tf << " " << Tm << " "  << -ST/A/(Tm-Tf+1e-9)  << "\n";
		}

        }

        csvfile.close();

	return 0;
}


// ************************************************************************* //
