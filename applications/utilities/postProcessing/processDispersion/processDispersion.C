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
    processDispersion 

Description
    calculate dispersion coefficient using closure variable B

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


        string filename = "disp.csv";
	std::ofstream csvfile(filename.c_str());
	csvfile << "time Dxx Dyy Dzz Dxy Dxz Dyz\n";

        
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
                    postProcessDict.lookupOrDefault<scalar>("x1",-1e30)
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


              Info<< "Reading molecular diffusionvity DT\n" << endl;

              dimensionedScalar DT("DT",transportProperties);

              dimensionedScalar gamma("gamma",transportProperties);

              volVectorField B
              (
                  IOobject
                  (
                    "B",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                  ),
                  mesh
                );

                volVectorField U
                (
                        IOobject
                        (
                                "U",
                                runTime.timeName(),
                                mesh,
                                IOobject::MUST_READ,
                                IOobject::NO_WRITE
                        ),
                        mesh
                );

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
                    mesh,
                    dimensionedScalar("eps",dimless,1.0)
                );


                volScalarField tinv
                (
                    IOobject
                    (
                        "tinv",
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("tinv",dimless,1.0)
                );

                volScalarField beta
                (
                    IOobject
                    (
                        "beta",
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("beta",dimless,0.0)
                );

                volScalarField Lpore
                (
                    IOobject
                    (
                         "Lpore",
                         runTime.timeName(),
                         mesh,
                         IOobject::READ_IF_PRESENT,
                         IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("Lpore",dimLength,0.0)
                );

                //U = fvc::average(fvc::interpolate(U));
                volScalarField Pe = mag(U)/0.429463*Lpore/DT;
                volTensorField Deff = DT*tinv*tensor(I);
	        Deff.replace(tensor::XX,DT*tinv*(1+0*beta*pow(Pe,gamma)));
                //Deff.replace(tensor::YY,DT*tinv*(1+beta*pow(Pe,gamma)));


                //Deff.replace(tensor::XY,DT*tinv*(beta/0.15*pow(Pe,0.57)));
                //Deff.replace(tensor::YX,DT*tinv*(beta/0.15*pow(Pe,0.57)));


                dimensionedScalar poro ("poro",eps.weightedAverage(clip));

                dimensionedVector Umoy ("Umoy", U.weightedAverage(clip));

                Info<< "Umoy = " << Umoy.value() << nl << endl;

                Info << "poro = " << poro.value() << nl << endl;

                volVectorField b = B -sum(eps*B*clip)/sum(eps*clip);//fvc::domainIntegrate(B)/sum(mesh.V());
                //volVectorField b = B -fvc::domainIntegrate(B)/sum(mesh.V());


                dimensionedTensor D ("D", dimensionSet(0,2,-1,0,0,0,0), tensor::zero);

                volTensorField Utildab ("Utildab",((U-eps/poro*Umoy)*b));

                volTensorField gradb ("gradb", Deff & fvc::grad(b));

                D = sum(Deff*eps*clip)/poro/sum(clip) - sum(Utildab*clip)/poro/sum(clip) + sum(gradb*eps*clip)/poro/sum(clip);


                dimensionedTensor Dsym = D;//0.5*(D+D.T());
                scalar Dxx = Dsym.component(tensor::XX).value();
                scalar Dyy = Dsym.component(tensor::YY).value();
                scalar Dzz = Dsym.component(tensor::ZZ).value();
                scalar Dxy = Dsym.component(tensor::XY).value();
                scalar Dyx = Dsym.component(tensor::YX).value();
                scalar Dxz = Dsym.component(tensor::XZ).value();
                scalar Dzx = Dsym.component(tensor::ZX).value();
                scalar Dyz = Dsym.component(tensor::YZ).value();
                scalar Dzy = Dsym.component(tensor::ZY).value();


		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << Dxx << " " << Dxy << " " << Dxz << " " << Dyx << " " << Dyy << " " << Dyz << " " << Dzx << " " << Dzy << " " << Dzz  << "\n";
		}

        }

        csvfile.close();

	return 0;
}


// ************************************************************************* //
