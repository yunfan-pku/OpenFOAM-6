/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    sprayFoam

Description
    Transient solver for compressible, turbulent flow with a spray particle
    cloud.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "basicSprayCloud.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include <fstream>
#include <sstream>
#include "smallthermo.C"
#include "csvfile.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
#include "postProcess.H"

#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createControl.H"
#include "createTimeControls.H"
#include "createFields.H"
#include "createFieldRefs.H"
#include "compressibleCourantNo.H"
#include "setInitialDeltaT.H"
#include "initContinuityErrs.H"


    turbulence->validate();
    int timei = 0;
    int loop1 = 0;
    int loop2 = 0;
    //Time clockTime_vle(runTime);
    double thermotime = 0;
    forAll(xcoord, i)
    {
        xcoord[i] = mesh.C()[i].x();
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n"
        << endl;

    while (runTime.run())
    {

        timei++;
#include "readTimeControls.H"
#include "compressibleCourantNo.H"
#include "setDeltaT.H"

        runTime++;

        std::stringstream sstream1;
        sstream1 << timei << "_ptest.csv";
        std::string filename;
        sstream1 >> filename;
        //csvwrite(filename, xcoord, "x");

        Info << "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        if (!pimple.frozenFlow())
        {
#include "rhoEqn.H"
            loop1 = 0;

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {

                loop1++;
                /*
                std::stringstream sstream1;
                sstream1<<timei<<"_"<<loop1<<"_"<<loop2<<"Us.csv";
                std::string s;
                sstream1>>s;
                std::ofstream foutl1(s);
                */

#include "UEqn.H"
#include "YEqn.H"
#include "EEqn.H"

                //foutl1.close();

                loop2 = 0;

                // --- Pressure corrector loop
                while (pimple.correct())
                {

                    loop2++;
                    std::stringstream sstream1;
                    sstream1 << "_" << timei << "_" << loop1 << "_" << loop2 << "_";
                    std::string s;
                    sstream1 >> s;



                    /*
                                        std::ofstream fout(s);

                                        forAll(p,i)
                                        {
                                        fout<<p[i]<<",";
                                        }
                                        fout<<std::endl;
                                        */
#include "pEqn.H"


                                        /*
                                        forAll(p,i)
                                        {
                                        fout<<p[i]<<",";
                                        }
                                        fout<<std::endl;
                                        fout.close();
                                        */
                }
                thermo.correct();
                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            rho = thermo.rho();

            if (runTime.write())
            {
                combustion->Qdot()().write();
            }
        }
        else
        {
            if (runTime.writeTime())
            {
                parcels.write();
            }
        }
        Info << "VLE time =" << thermotime << endl;
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End\n"
        << endl;

    return 0;
}

// ************************************************************************* //
