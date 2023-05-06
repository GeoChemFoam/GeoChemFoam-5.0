
if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0) && (deltaEpsMax > maxDeltaEps))
    {
        runTime.setDeltaT
        (
            min
            (
                maxDeltaEps*runTime.deltaT().value()/deltaEpsMax,
                maxDeltaT
            )
        );
    }
}

// ************************************************************************* //
