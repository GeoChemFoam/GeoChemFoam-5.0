if (adjustTimeStep)
{
    scalar maxDeltaTFact = maxMeshCo/(meshCoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    scalar newTime = runTime.value()+runTime.deltaTValue();
    if (newTime>runTime.endTime().value()) runTime.setDeltaT(runTime.endTime().value()-runTime.value(),false); 
    else if (newTime+0.5*runTime.deltaTValue() > runTime.endTime().value()) runTime.setDeltaT((runTime.endTime().value()-runTime.value())/2.0,false); 
    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
