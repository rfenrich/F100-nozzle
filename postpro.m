meshGMF = ReadGMF('axinoz.mesh');
meshGMF = meshPrepro(meshGMF);
meshSU2 = convertGMFtoSU2(meshGMF);

SolSU2 = ReadSU2Sol('restart_flow.dat');

nozzle = nozzleCFDPostPro(meshSU2, SolSU2, nozzle, fluid, freestream);