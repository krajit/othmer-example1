/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield                | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O perati       on     |
    \\  /    A nd                  | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipul       ation  |
---------------------       ----------------------------------------------------------
License       
    This file is part        of OpenFOAM.
       
    OpenFOAM is free        software: you can redistribute it and/or modify it
    under the terms o       f the GNU General Public License as published by
    the Free Software        Foundation, either version 3 of the License, or
    (at your option)        any later version.
       
    OpenFOAM is distr       ibuted in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; wit       hout even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PAR       TICULAR PURPOSE.  See the GNU General Public License
    for more details.       

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

template <class Type>
void zeroCells(
    GeometricField<Type, fvPatchField, volMesh> &vf,
    const labelList &cells)
{
    forAll(cells, i)
    {
        vf[cells[i]] = Zero;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "postProcess.H"

#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createControl.H"
#include "createFields.H"
#include "initContinuityErrs.H"
#include "initAdjointContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Cost function value
    scalar J = 0;
    scalar Jold = 0;
    scalar Jk = 0;

    // Compute cost function value
#include "costFunctionValue.H"

    std::ofstream file("results.csv");
    file << 0 << "," << J << nl;
    file.close();

    Info << "\nStarting time loop\n"
         << endl;

    while (simple.loop() && (fabs(J - Jold) > tol) && (gamma > tol))
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // save old cost value
        Jold = J;

#include "stateEquation.H"
#include "adjointEquation.H"

        // Save current control
        alphak = alpha;

// calculate current cost
#include "costFunctionValue.H"
        Jk = J;

        bool gammaFound = false;

        // calculate derivative^2 integrate(U . Ua dv). Why??
        scalar phip0 = gSum(volField * (Ua.internalField() & U.internalField()));

        while ((!gammaFound) && (gamma > tol))
        {
            alpha = alpha - gamma * (Ua & U);

            // truncate u for constrained control set
            forAll(alpha, i)
            {
                alpha[i] = min(alpha[i], alphaMax[i]);
                alpha[i] = max(alpha[i], alphaMin[i]);
            }
            alpha.correctBoundaryConditions();

            // get new u
            #include "adjointEquation.H"
            // get new cost
            #include "costFunctionValue.H"

            // backtracking step to find alpha
            if (J <= Jk - c1 * gamma * phip0)
            {
                Info << "gamma found, gamma = " << gamma << ", J = " << J << ", phip0" << phip0 << endl;
                gammaFound = true;
            }
            else
            {
                Info << "gamma NOT found, gamma = " << gamma << ", J = " << J << ", phip0" << phip0 << endl;
                gamma = c2 * gamma;
            }
        }

        file.open("results.csv", std::ios::app);
        file << runTime.value() << "," << J << nl;
        file.close();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
