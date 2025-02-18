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
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutput.h"
#include "math/RegularRankFourTensor.h"
#include "system/SystemController.h"
#include "geometry/Geometry.h"
#include "parameters/Constants.h"

#include "basis/CartesianToSphericalTransformer.h"

namespace Serenity {

WriteDataTask::WriteDataTask(std::shared_ptr<SystemController> system) : 
_system(system),
_basis(system->getBasisController()),
_orbitals(system->getActiveOrbitalController<RSCF>()),
_elstruct(system->getElectronicStructure<RSCF>()),
_geometry(system->getGeometry()),
_coeffs(std::make_unique<Eigen::MatrixXd>(_system->getActiveOrbitalController<RSCF>()->getCoefficients())),
_eris(nullptr){
}

void WriteDataTask::run() {
  printSubSectionTitle("Saving SCF data to hdf5 file");
  OutputControl::n.printf(
    "  This module will save orbital energies, coefficients, and molecular integrals to hdf5 file.\n\n");

  _nBasisFunc = _basis->getNBasisFunctions();
  const unsigned int nOccOrb = _system->getNOccupiedOrbitals<RSCF>();

  //default active space settings
  _nActOrbs = _nBasisFunc;
  _iActOrb = 0;
  _fActOrb = _nBasisFunc;//-1;
  _nOccActOrbs = nOccOrb;

  if(!settings.activeOrbs.empty()){
    if(settings.activeOrbs.size() == 2){//initial and final orbitals are specified
      _iActOrb = settings.activeOrbs[0]-1;
      _fActOrb = settings.activeOrbs[1];//-1;
      if(_iActOrb>=_fActOrb){
        throw SerenityError(
          "SAVEDATA task: Wrong active space is requested.\n"
          "               Index of the initial orbital must be lower than that of the final orbital!");
      }
      if(_iActOrb>=_nBasisFunc || _fActOrb>_nBasisFunc){
          throw SerenityError(
            "SAVEDATA task: Wrong active space is requested.\n"
            "               Index of the initial or final orbitals must be lower than the size of the basis!");
      }
      if(_iActOrb<=0 || _iActOrb>=nOccOrb){
          throw SerenityError(
            "SAVEDATA task: Wrong active space is requested.\n"
            "               Index of the initial orbital must be in the occupied space!");
      }
      if(_fActOrb<=nOccOrb || _fActOrb>_nBasisFunc){
          throw SerenityError(
            "SAVEDATA task: Wrong active space is requested.\n"
            "               Index of the final orbital must be in the virtual space!");
      }    
      _nActOrbs = _fActOrb - _iActOrb;
      _nOccActOrbs = nOccOrb - _iActOrb;
    }
    else{
      throw SerenityError("SAVEDATA task: Wrong syntax in activeOrbs keyword. Use {i f}.");
    }  
  }

  //const std::string baseName = "data";
  //const std::string fileName = _system->getSystemName() + "." + baseName + ".hdf5";
  const std::string fileName = settings.fname;
  OutputControl::nOut << "  Writing data to file " << fileName << "\n" << std::endl;

  //create hdf5 file stream
  HDF5::H5File file(fileName.c_str(), H5F_ACC_TRUNC);

  //extract information about the system  
  const double energy = _elstruct->getEnergy();

  //write orbital energies
  OutputControl::nOut << "  Saving orbital energies..." << std::endl;
  const auto& orbitalEnergies = _orbitals->getEigenvalues();
  //bug fixed: syntax of .segment(i,n): Block containing n elements, starting at position i 
  HDF5::save(file, "energies", orbitalEnergies.segment(_iActOrb,_nActOrbs));

  //write orbital occupations
  OutputControl::nOut << "  Saving orbital occupations..." << std::endl;
  const auto& occs = _elstruct->getDensityMatrixController()->getOccupations();
  HDF5::save(file, "occs", occs.segment(_iActOrb,_nActOrbs));

  //write orbital irreps
  OutputControl::nOut << "  Saving orbital symmetries..." << std::endl;
  Eigen::VectorXi irreps(_nActOrbs);
  irreps.setOnes();//symmetry is not implemented in serenity, thus fill with ones
  HDF5::save(file, "irreps", irreps);

  //write orbital coefficients
  OutputControl::nOut << "  Saving orbital coefficients..." << std::endl;
  Eigen::VectorXd coeffs(_nActOrbs*_nBasisFunc);
  this->prepCoefficients(coeffs);
  HDF5::save(file, "coefficients", coeffs);

  //prepare and save MO integrals and their indices
  OutputControl::nOut << "  Saving integrals..." << std::endl;
  Eigen::VectorXd integrals;
  Eigen::VectorXi indices;
  this->prepERIS(integrals,indices);
  HDF5::save(file, "integrals", integrals);
  HDF5::save(file, "indices", indices);
  const unsigned int nIntegrals = integrals.size();

  //prepare and save AO basis
  OutputControl::nOut << "  Saving AO basis..." << std::endl;  
  unsigned int nAO;
  unsigned int maxNmbCC;
  Eigen::VectorXi ncc;
  Eigen::VectorXd cc;
  Eigen::VectorXd alpha;
  Eigen::VectorXi center;
  Eigen::VectorXi polynomial;
  prepAOBasis(nAO,maxNmbCC,ncc,cc,alpha,center,polynomial);
  HDF5::save_scalar_attribute(file, "nAO", nAO);
  HDF5::save_scalar_attribute(file, "maxNmbCC", maxNmbCC);
  HDF5::save(file, "ncc", ncc);  
  HDF5::save(file, "cc", cc);
  HDF5::save(file, "alpha", alpha);
  HDF5::save(file, "center", center);
  HDF5::save(file, "polynomial", polynomial);

  //prepare and save AO overlap matrix
  Eigen::MatrixXd sAO = _system->getOneElectronIntegralController()->getOverlapIntegrals();
  std::cout<<"Size: "<<sAO.cols()<<" "<<sAO.rows()<<std::endl;
  Eigen::VectorXd sAO_data(Eigen::Map<Eigen::VectorXd>(sAO.data(), sAO.cols()*sAO.rows()));
  HDF5::save(file, "sAO", sAO_data);

  //prepare and save molecular geometry
  OutputControl::nOut << "  Saving geometry..." << std::endl;
  unsigned int nAtoms = 0;
  Eigen::VectorXd coords;
  Eigen::VectorXd Zs;
  double Enuc = 0.;
  prepGeometry(nAtoms,coords,Zs,Enuc);
  HDF5::save_scalar_attribute(file, "nAtoms", nAtoms);
  HDF5::save_scalar_attribute(file, "Enuc", Enuc);
  HDF5::save(file, "geometry", coords);  
  HDF5::save(file, "Zs", Zs);

  //saving attributes
  OutputControl::nOut << "  Saving attributes..." << std::endl;    
  HDF5::save_scalar_attribute(file, "nSym", 1);                         //number of irreps (1 in serenity)
  HDF5::save_scalar_attribute(file, "nCenters", _system->getNAtoms());  //number of atoms
  HDF5::save_scalar_attribute(file, "nBasisFunc", _nBasisFunc);
  HDF5::save_scalar_attribute(file, "nActOrbs", _nActOrbs);
  HDF5::save_scalar_attribute(file, "nOccOrb", _nOccActOrbs);
  HDF5::save_scalar_attribute(file, "energy", energy);
  HDF5::save_scalar_attribute(file, "nIntegrals", nIntegrals);

  OutputControl::nOut << "\n" << std::endl;
  OutputControl::nOut << "  Total energy: "<<energy<< std::endl;  
  OutputControl::nOut << "  Active orbitals: "<<_iActOrb+1<<"..."<<_fActOrb<< std::endl;
  OutputControl::nOut << "  Number of active orbitals: "<<_nActOrbs<< std::endl;
  OutputControl::nOut << "  Number of occupied orbitals: "<<_nOccActOrbs<< std::endl;
  OutputControl::nOut << "  Number of non-zero integrals: "<<nIntegrals<< std::endl;

  file.close();

  OutputControl::nOut << std::endl;
  OutputControl::nOut << "  Done!" << std::endl;
  OutputControl::nOut << std::endl;

  return;
}

void WriteDataTask::prepCoefficients(Eigen::VectorXd& coeffs){
  for (unsigned int i = _iActOrb; i < _fActOrb; ++i) {
    for (unsigned int j = 0; j < _nBasisFunc; ++j) {
      unsigned int i_ = i - _iActOrb; //shifted orbital index to account for active space selection
      coeffs(i_ * _nBasisFunc + j) = (*_coeffs)(j,i);   //we reverse the indices to match theADCcode convention
    }
  }
  return;
}

void WriteDataTask::prepERIS(Eigen::VectorXd& integrals, Eigen::VectorXi& indices){
  if (!_eris) {
    _eris = std::unique_ptr<RegularRankFourTensor<double>>(new RegularRankFourTensor<double>(_nBasisFunc, 0.0));
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

  //TODO: it makes sense to rewrite the AO2MO transformation routine to perform it only for active orbitals
  ao2mo.transformTwoElectronIntegrals(eris, eris, *_coeffs, _nBasisFunc);

  std::vector<double> ints;
  std::vector<int> indx;

  for (unsigned int i = _iActOrb; i < _fActOrb; ++i) {
    for (unsigned int j = i; j < _fActOrb; ++j) {
      for (unsigned int k = i; k < _fActOrb; ++k) {
        for (unsigned int l = k; l < _fActOrb; ++l) {
          double val=eris(i,j,k,l);
          if(fabs(val) >= settings.eriThreshold){
            ints.push_back(val);
            indx.push_back(i-_iActOrb);
            indx.push_back(j-_iActOrb);
            indx.push_back(k-_iActOrb);
            indx.push_back(l-_iActOrb);
          }
        }
      }
    }
  }

  integrals = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ints.data(), ints.size());
  indices   = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(indx.data(), indx.size());

  return;
}

//find to what atom a given shell belongs
unsigned int WriteDataTask::find_atom_indx(const double& x,
                                           const double& y,
                                           const double& z){

  unsigned int nAtoms = _geometry->getNAtoms();

  const auto atoms = _geometry->getAtoms();

  unsigned int ia = 0;
  for(ia = 0; ia<nAtoms; ia++){
    const auto& atom = atoms[ia];

    const double xa = atom->getX();
    const double ya = atom->getY();
    const double za = atom->getZ();

    double dist = sqrt((x-xa)*(x-xa) + (y-ya)*(y-ya) + (z-za)*(z-za));

    if(dist<1.E-5)
      break;
  }

  return ia;
}

void WriteDataTask::prepAOBasis(
  unsigned int& nAO, unsigned int& maxNmbCC,
  Eigen::VectorXi& ncc_,
  Eigen::VectorXd& cc_,
  Eigen::VectorXd& alpha_,
  Eigen::VectorXi& center_,
  Eigen::VectorXi& polynomial_)
{
  size_t nShells = _basis->getReducedNBasisFunctions(); //number of atomic orbitals
  maxNmbCC = _basis->getMaxNumberOfPrimitives(); //maximal number of contractions

  const auto shells = _basis->getBasis();

  //calculate number of AO
  nAO = 0;
  for(size_t i=0; i<nShells; i++){
      auto shell = shells[i];
      auto l = shell->getAngularMomentum();
      nAO += N_SHELL_CART[l];
  }

  std::vector<int> ncc(nAO);
  std::vector<double> cc(nAO*maxNmbCC,0.0);
  std::vector<double> alpha(nAO*maxNmbCC,0.0);
  std::vector<int> center(nAO);
  std::vector<int> polynomial(nAO*3);

  size_t indx = 0;
  for (unsigned int i = 0; i < nShells; ++i) {
    auto shell = shells[i];

    auto ncnts = shell->getNContracted();
    auto nexps = shell->getNPrimitives();
    auto l = shell->getAngularMomentum();

    //extract exponents and contraction coefficients for the given shell
    const auto exponents = shell->alpha;
    const auto contractions = shell->contr[0].coeff;

    //coordinates of the shell
    const double x = shell->getX();
    const double y = shell->getY();
    const double z = shell->getZ();

    //find to which atom the shell belongs
    auto atom_indx = find_atom_indx(x,y,z);

    //generate various l projections
    for (int aX = l; aX >= 0; --aX) {
      for (int aY = l - aX; aY >= 0; --aY) {
        int aZ = l - aY - aX;

        //assign values of the polynomial
        polynomial[3*indx    ] = aX;
        polynomial[3*indx + 1] = aY;
        polynomial[3*indx + 2] = aZ;

        //assign number of contractions and atom index
        ncc[indx] = nexps;
        center[indx] = atom_indx + 1;

        //normalize the AO
        long double norm = 
            std::pow(
                (long double)(double_factorial(2 * aX - 1)) 
              * (long double)(double_factorial(2 * aY - 1)) 
              * (long double)(double_factorial(2 * aZ - 1)),-0.5)
          * std::pow(double_factorial(2 * l - 1), 0.5);

        //assign exponentials and contraction coefficients
        for(size_t icc=0; icc<nexps; icc++){
            cc   [indx*maxNmbCC + icc] = norm * contractions[icc];
            alpha[indx*maxNmbCC + icc] = exponents[icc];
        }

        indx++;
      }
    }
  }

  ncc_    = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(ncc.data(), ncc.size());
  cc_     = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(cc.data(), cc.size());
  alpha_  = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(alpha.data(), alpha.size());
  center_ = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(center.data(), center.size());
  polynomial_ = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(polynomial.data(), polynomial.size());

  return;
}

void WriteDataTask::prepGeometry(
  unsigned int& nAtoms_,
  Eigen::VectorXd& coords_,
  Eigen::VectorXd& Zs_,
  double& Enuc_)
{
  Enuc_ = _geometry->getCoreCoreRepulsion();
  nAtoms_ = _geometry->getNAtoms();

  const auto atoms = _geometry->getAtoms();
  
  std::vector<double> coords(3*nAtoms_,0.0);
  std::vector<double> Zs(nAtoms_,0.0);

  for(size_t ia = 0; ia<nAtoms_; ia++){
    const auto& atom = atoms[ia];

    const int Z = atom->getNuclearCharge();
    const double x = atom->getX();
    const double y = atom->getY();
    const double z = atom->getZ();

    coords[3 * ia ]    = x;
    coords[3 * ia + 1] = y;
    coords[3 * ia + 2] = z;
    Zs[ia] = (double)Z;
  }

  coords_ = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(coords.data(), coords.size());
  Zs_     = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Zs.data(), Zs.size());

  return;
}

} /* namespace Serenity */
