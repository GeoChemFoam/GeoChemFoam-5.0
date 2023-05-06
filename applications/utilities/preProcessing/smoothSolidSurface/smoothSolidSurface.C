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
    smoothSolidSurface

Description
    smooth solid surface volume fraction eps 

\*---------------------------------------------------------------------------*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "fvCFD.H"
#include "argList.H"
#include "primitivePatchInterpolation.H"
#include "timeSelector.H"
#include "simpleControl.H"


using namespace Foam;

int main(int argc, char *argv[])
{
	#   include "setRootCase.H"
	#   include "createTime.H"

	instantList timeList = timeSelector::select0(runTime, args);

	forAll(timeList, timeStep)
	{
		runTime.setTime(timeList[timeStep], timeStep);

		Info<< endl<<timeStep<< "    Time = " << runTime.timeName() << endl;

		#   include "createNamedMesh.H"
		
                simpleControl simple(mesh);

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

		volScalarField eps 
		(
			IOobject
			(
				"eps",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);

		dimensionedScalar kf(dimensionedScalar::getOrDefault("kf",transportProperties,dimless/sqr(dimLength),Zero));

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
                    kf*pow(1-eps+1e-13,2)/pow(eps,3)
                );

		volScalarField kfvol = Kinv*pow(eps,3)/pow(1-eps+1e-13,2);

		
		scalar cS = readScalar(simple.dict().lookup("cS"));
		label nS = readLabel(simple.dict().lookup("nS"));
		for(int i=0;i<nS;i++)
		{
			eps = cS*fvc::average(fvc::interpolate(eps))+(1-cS)*eps;
		}
    
		Kinv = kfvol*pow(1-eps+1e-13,2)/pow(eps,3);

		eps.write();
		Kinv.write();

	}

	return 0;
}


// ************************************************************************* //
