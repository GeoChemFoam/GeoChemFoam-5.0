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
    processPoroPerm

Description
    calculate porosity and permeability from flow field for each time-step

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
 
        scalar poro = 0.0;
        scalar perm = 0.0;

	std::ofstream csvfile("poroPerm.csv");
	csvfile << "time poro perm Lpore Re UD\n";

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

		dimensionedScalar nu("nu",dimViscosity,transportProperties);

        word velocity
        (
            postProcessDict.lookup("velocity")
        );
 
        word pressureDrop 
        (
            postProcessDict.lookup("pressureDrop")
        );


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
			mesh,
		    dimensionedScalar("eps",dimless,1.0)
		);

                volScalarField Kinv
                (
                    IOobject
                    (
                        "Kinv",
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("Kinv",dimless/dimArea,0.0)
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

                volScalarField p
                (
                        IOobject
                        (
                                "p",
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

                volScalarField Ud = U.component(direction);
                scalar Uavg =  Ud.weightedAverage(clip).value();
                scalar UD=Uavg*vol/(x2-x1)/(y2-y1)/(z2-z1);


               if (velocity=="inlet")
               {
		 Info<< "Velocity is calcluated at the inlet"
                       << endl;

                   word inletName
                   (
                       postProcessDict.lookup("inletName")
                   );

                   label inletPatchID  = mesh.boundaryMesh().findPatchID(inletName);

                   const fvPatchVectorField& Uf = U.boundaryField()[inletPatchID];
                   scalarField Ufd = Uf.component(direction);

                   const fvPatchScalarField& pInlet  = p.boundaryField()[inletPatchID];
                   scalarField magSfInlet  = pInlet.patch().magSf();

                   scalar Q  = gSum(Ufd*magSfInlet);
                   scalar A  = (x2-x1)*(y2-y1)*(z2-z1)/L;
                   UD = Q/A;

               }
               else if (velocity=="volumeAveraged")
               {
		Info<< "Velocity is volume averaged"
                       << endl;

                  //UD=Uavg*vol/(x2-x1)/(y2-y1)/(z2-z1);
               }
               else
               { 
                   Info<< "velocity should be equal to inlet or volumeAveraged"
                       << endl
                       << abort(FatalError);
               }

                volScalarField gradRate= fvc::grad(p) & U;
                scalar gradP = gradRate.weightedAverage(clip).value();

                scalar dpdz= -gradP/(Uavg+1e-30);


                if (pressureDrop=="direct")
                {
			Info<< "pressure drop is calculated directly"
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

                    const fvPatchScalarField& pInlet  = p.boundaryField()[inletPatchID];
                    const fvPatchScalarField& pOutlet = p.boundaryField()[outletPatchID];
                    scalarField magSfInlet  = pInlet.patch().magSf();
                    scalarField magSfOutlet = pOutlet.patch().magSf();

                    dpdz = mag(gSum(pInlet*magSfInlet)/gSum(magSfInlet)-gSum(pOutlet*magSfOutlet)/gSum(magSfOutlet))/L+1e-13;

                }
                else if (pressureDrop=="velocityWeightedAverageGradient")
                {
			Info<< "Pressure drop is calculated by a velocity weighted average gradient"
                       << endl;

                    //dpdz= -gradP/(Uavg+1e-30);
                }
                else if (pressureDrop=="velocityWeightedAverageViscousForce")
                {
			Info<< "Pressure drop is calculated by a velocity weighted average viscous force"
                       << endl;

                    volScalarField viscPRate=(nu/eps)*(fvc::laplacian(U) & U)-nu*Kinv*(U & U);

                    scalar viscP = viscPRate.weightedAverage(clip).value();

                    dpdz=-viscP/(Uavg+1e-30);
                }
                else if (pressureDrop=="momentumSource")
                {
                    //dimensionedVector momentumSource("momentumSource",dimAcceleration,transportProperties);
		    Info<< "Reading pressure gradient\n" << endl;
                    IOdictionary fvOptions
                    (
                        IOobject
                        (
                            "fvOptions",
                            runTime.constant(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        )
                    );
		    const dictionary& gradPDict =fvOptions.subDict("momentumSource").subDict("injectionRateSuSp");
                    dimensionedVector momentumSource = Tuple2<vector, scalar>(gradPDict.lookup("U")).first();
		    dpdz = momentumSource.component(direction).value();
                }
               else
               {
                   Info<< "pressureDrop should be equal to direct or velocityWeightedAverageGradient or velocityWeightedAverageViscousForce or momentumSource"
                       << endl
                       << abort(FatalError);
               }



		poro = eps.weightedAverage(clip).value()*vol/(x2-x1)/(y2-y1)/(z2-z1);

                perm = UD*nu.value()/(dpdz+1e-30);

                if (perm<0)
                {
                   Info<< "negative permeability - did the velocity field converge?"
                       << endl
                       << abort(FatalError);
                
                }
                //scalar perm = -Uavg*vol/(x2-x1)/(y2-y1)/(z2-z1)*nu.value()/(dpdz+1e-30); 
		scalar Lpore = std::sqrt(8*perm/poro);
		scalar Re = UD/poro*Lpore/nu.value();


		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << poro << " " << perm << " " << Lpore << " " << Re << " " << UD << "\n";
		}
	}

	csvfile.close();

        std::ofstream SPdict("system/SPDict");
        SPdict << "FoamFile\n";
        SPdict << "{\n";
        SPdict << "version     4.0;\n";
        SPdict << "format      ascii;\n";
        SPdict << "class       dictionary;\n";
        SPdict << "location    \"system\";\n";
        SPdict << "object      postProcessDict;\n";
        SPdict << "}\n";
        SPdict << "\n";
        SPdict << "poro " << poro  << " ;\n";
        SPdict << "perm " << perm << " ;\n";


	return 0;
}


// ************************************************************************* //
