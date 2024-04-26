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
    processFF

Description
    calculate porosity and formation factor

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
        #   include "createNamedMesh.H"

	instantList timeList = timeSelector::select0(runTime, args);
 
        scalar poro = 0.0;
        scalar Deff = 0.0;
        scalar FF   = 0.0;

	std::ofstream csvfile("FF.csv");
	csvfile << "time poro FF\n";

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
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
		);

		IOdictionary transportProperties 
		(
			IOobject
			(
				"transportProperties",
				"constant",
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
		);

		dimensionedScalar DT("DT",dimViscosity,transportProperties);


		label direction
		(
		    readLabel(postProcessDict.lookup("direction"))
		);

                scalar z1
                (
                    readScalar(postProcessDict.lookup("z1"))
                );


                scalar z2
                (
                    readScalar(postProcessDict.lookup("z2"))
                );

                scalar y1
                (
                    readScalar(postProcessDict.lookup("y1"))
                );


                scalar y2
                (
                    readScalar(postProcessDict.lookup("y2"))
                );


                scalar x1
                (
                    readScalar(postProcessDict.lookup("x1"))
                );


                scalar x2
                (
                   readScalar(postProcessDict.lookup("x2"))
                );

                word gradient
                (
                    postProcessDict.lookup("gradient")
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

                scalar vol = gSum(clip);

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

                volScalarField T
                (
                        IOobject
                        (
                                "T",
                                runTime.timeName(),
                                mesh,
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE
                        ),
                        mesh
                );

                scalar L=0.0;

                if (direction==0) L=x2-x1;
                if (direction==1) L=y2-y1;
                if (direction==2) L=z2-z1;

                poro = eps.weightedAverage(clip).value()*vol/(x2-x1)/(y2-y1)/(z2-z1);

                volVectorField J = eps*DT*fvc::grad(T);

                volScalarField Jd = J.component(direction);
                scalar Javg =  Jd.weightedAverage(clip).value();
                scalar JD   = Javg*vol/(x2-x1)/(y2-y1)/(z2-z1);


                volScalarField gradRate= fvc::grad(T) & J;
                scalar gradR = gradRate.weightedAverage(clip).value();

                scalar dTdz= -gradR/Javg;

                if (gradient=="direct")
                {
                    Info<< "gadient is calculated directly"
                       << endl;

                    word inletName
                    (
                        postProcessDict.lookup("inletName")
                    );

                    word outletName
                    (
                        postProcessDict.lookup("outletName")
                    );

                    label inletPatchID  = mesh.boundaryMesh().findPatchID(inletName);
                    label outletPatchID = mesh.boundaryMesh().findPatchID(outletName);

                    const fvPatchScalarField& TInlet  = T.boundaryField()[inletPatchID];
                    const fvPatchScalarField& TOutlet = T.boundaryField()[outletPatchID];
                    scalarField magSfInlet  = TInlet.patch().magSf();
                    scalarField magSfOutlet = TOutlet.patch().magSf();

                    dTdz = mag(gSum(TInlet*magSfInlet)/gSum(magSfInlet)-gSum(TOutlet*magSfOutlet)/gSum(magSfOutlet))/L+1e-13;
                }
    
                else if (gradient=="volumeAveraged")
                {
                    Info<< "gradient is calculated by voolume averaged"
                      << endl;
                }

               else
               {
                   Info<< "gradient should be equal to direct or volumeAveraged"
                       << endl
                       << abort(FatalError);
               }


                Deff = -JD/dTdz;
                FF = DT.value()/Deff;

		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << poro << " " << FF << "\n";
		}
	}

	csvfile.close();

	return 0;
}


// ************************************************************************* //
