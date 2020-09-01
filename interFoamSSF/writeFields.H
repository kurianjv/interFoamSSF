Info<<"Volume integral of bubble: " << fvc::domainIntegrate(alpha2) << endl;
Info<<"Volume integral of water: " << fvc::domainIntegrate(alpha1) << endl;

Info  << "   Uavg (m/s) = " << mag(average(U)).value() << endl;
Info  << "   Umax (m/s) = " << max(mag(U)).value() << endl; 

Info  << "   p_rgh in bubble (kg/ms2) = " << fvc::domainIntegrate(alpha2*p_rgh)/fvc::domainIntegrate(alpha2) << endl; 
Info  << "   p in bubble (kg/ms2) = " << fvc::domainIntegrate(alpha2*p)/fvc::domainIntegrate(alpha2) << endl; 

Info << "Total water flow at out of domain"   <<  fvc::domainIntegrate(fvc::div(phi)) << endl;


volScalarField Ux = U.component(0);
volScalarField Uy = U.component(1);
volScalarField Uz = U.component(2);

Info  << "   Bubble avg. Ux (m/s) = " << fvc::domainIntegrate(alpha2*Ux)/fvc::domainIntegrate(alpha2) << endl;
Info  << "   Bubble avg. Uy (m/s) = " << fvc::domainIntegrate(alpha2*Uy)/fvc::domainIntegrate(alpha2) << endl;
Info  << "   Bubble avg. Uz (m/s) = " << fvc::domainIntegrate(alpha2*Uz)/fvc::domainIntegrate(alpha2) << endl;
Info  << "   Bubble avg. U (m/s) = " << fvc::domainIntegrate( alpha2*mag(U) )/fvc::domainIntegrate(alpha2) << endl;


