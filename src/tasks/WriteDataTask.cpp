/**
 * @file   WriteDataTask.cpp
 *
 * @date   Oct. 18, 2024
 * @author Nikolay Golubev
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "tasks/WriteDataTask.h"
/* Include Serenity Internal Headers */
#include "io/FormattedOutputStream.h"  // Filtered output streams.
#include "io/HDF5.h"                   // Write HDF5 files.
#include "system/SystemController.h"   // Access system name.
/* Include Std and External Headers */
#include <iomanip> // Format ASCII output files.

#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/transformer/Ao2MoTransformer.h"
#include "io/FormattedOutput.h"
#include "math/RegularRankFourTensor.h"
#include "system/SystemController.h"

namespace Serenity {

WriteDataTask::WriteDataTask(std::shared_ptr<SystemController> system) : 
_system(system),
_basis(system->getBasisController()),
_orbitals(system->getActiveOrbitalController<RSCF>()),
_elstruct(system->getElectronicStructure<RSCF>()),
_coeffs(std::make_unique<Eigen::MatrixXd>(_system->getActiveOrbitalController<RSCF>()->getCoefficients())),
_eris(nullptr){
}

void WriteDataTask::run() {
  printSubSectionTitle("Saving SCF data to hdf5 file");
  OutputControl::n.printf(
    "  This module will save orbital energies, coefficients, and molecular integrals to hdf5 file.\n\n");

  if(settings.activeOrbs.size() == 0){//default active space - all orbitals

  }
  else if(settings.activeOrbs.size() == 2){//initial and final orbitals are specified

  }
  else{
    throw SerenityError("SAVEDATA task: Wrong syntax in activeOrbs keyword.");
  }
  

  //some info
  //unsigned int OrbitalIn=

  //extract information about the system
  const unsigned int nBasisFunc = _basis->getNBasisFunctions();
  const unsigned int nOccOrb = _system->getNOccupiedOrbitals<RSCF>();
  const double energy = _elstruct->getEnergy();

  //const std::string baseName = "data";
  //const std::string fileName = _system->getSystemName() + "." + baseName + ".hdf5";
  const std::string fileName = settings.fname;
  OutputControl::nOut << "  Writing data to file " << fileName << "\n" << std::endl;

  //create hdf5 file stream
  HDF5::H5File file(fileName.c_str(), H5F_ACC_TRUNC);

  //write orbital energies
  OutputControl::nOut << "  Saving orbital energies..." << std::endl;
  const auto& orbitalEnergies = _orbitals->getEigenvalues();
  HDF5::save(file, "energies", orbitalEnergies);

  //write orbital irreps
  OutputControl::nOut << "  Saving orbital symmetries..." << std::endl;
  Eigen::VectorXi irreps(nBasisFunc,1);//symmetry is not implemented in serenity, thus fill with ones
  HDF5::save(file, "irreps", irreps);

  //write orbital coefficients
  OutputControl::nOut << "  Saving orbital coefficients..." << std::endl;
  Eigen::VectorXd coeffs;
  this->prepCoefficients(coeffs);
  HDF5::save(file, "coefficients", coeffs);

  //prepare and save MO integrals and their indices
  OutputControl::nOut << "  Saving integrals..." << std::endl;
  Eigen::VectorXd integrals;
  Eigen::VectorXd indices;
  this->prepERIS(integrals,indices);
  HDF5::save(file, "integrals", integrals);
  HDF5::save(file, "indices", indices);
  const unsigned int nIntegrals = integrals.size();

  OutputControl::nOut << "\n" << std::endl;
  OutputControl::nOut << "  Size of the orbital basis: "<<nBasisFunc<< std::endl;
  OutputControl::nOut << "  Number of occupied orbitals: "<<nOccOrb<< std::endl;
  OutputControl::nOut << "  Total energy: "<<energy<< std::endl;
  OutputControl::nOut << "  Number of non-zero integrals: "<<nIntegrals<< std::endl;

  HDF5::save_scalar_attribute(file, "nBasisFunc", nBasisFunc);
  HDF5::save_scalar_attribute(file, "nOccOrb", nOccOrb);
  HDF5::save_scalar_attribute(file, "energy", energy);
  HDF5::save_scalar_attribute(file, "nIntegrals", nIntegrals);


  file.close();

  OutputControl::nOut << std::endl;
  OutputControl::nOut << "  Done!" << std::endl;
  OutputControl::nOut << std::endl;

  return;
}

void WriteDataTask::prepCoefficients(Eigen::VectorXd& coeffs){
  coeffs = *_coeffs;
  return;
}

void WriteDataTask::prepERIS(Eigen::VectorXd& integrals, Eigen::VectorXd& indices){
  const unsigned int nBasisFunc = _basis->getNBasisFunctions();
  if (!_eris) {
    _eris = std::unique_ptr<RegularRankFourTensor<double>>(new RegularRankFourTensor<double>(nBasisFunc, 0.0));
  }
  auto& eris = *_eris;

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, _basis, 1E-10);

  auto const storeERIS = [&eris](const unsigned int& a, const unsigned int& b, const unsigned int& i,
                                 const unsigned int& j, const Eigen::VectorXd& integral, const unsigned int threadId) {
    (void)threadId;
    eris(b, a, i, j) = integral(0);
    eris(b, a, j, i) = integral(0);
    eris(a, b, j, i) = integral(0);
    eris(a, b, i, j) = integral(0);
    eris(i, j, b, a) = integral(0);
    eris(i, j, a, b) = integral(0);
    eris(j, i, b, a) = integral(0);
    eris(j, i, a, b) = integral(0);
  };

  looper.loop(storeERIS, _coeffs->lpNorm<Eigen::Infinity>());

  Ao2MoTransformer ao2mo(_basis);

  ao2mo.transformTwoElectronIntegrals(eris, eris, *_coeffs, nBasisFunc);

  std::vector<double> ints;
  std::vector<double> indx;

  for (unsigned int i = 0; i < nBasisFunc; ++i) {
    for (unsigned int j = i; j < nBasisFunc; ++j) {
      for (unsigned int k = i; k < nBasisFunc; ++k) {
        for (unsigned int l = k; l < nBasisFunc; ++l) {
          double val=eris(i,j,k,l);
          if(fabs(val) >= settings.eriThreshold){
            ints.push_back(val);
            indx.push_back(i);
            indx.push_back(j);
            indx.push_back(k);
            indx.push_back(l);
          }
        }
      }
    }
  }

  integrals = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ints.data(), ints.size());
  indices   = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(indx.data(), indx.size());

  return;
}

} /* namespace Serenity */
