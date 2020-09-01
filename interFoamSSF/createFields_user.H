
    volVectorField fc
    (
        IOobject
        (
            "fc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("fc", dimMass/(dimLength*dimLength*dimTime*dimTime), vector(0,0,0))
    );

