/**
 * @file   WriteDataTask.h
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

#ifndef SERENITY_WRITEDATATASK_H
#define SERENITY_WRITEDATATASK_H

#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

#include "data/OrbitalController.h"
#include "data/ElectronicStructure.h"

#include "math/Matrix.h"
#include "math/RegularRankFourTensor.h"

namespace Serenity {

#define RSCF Options::SCF_MODES::RESTRICTED

/* Forward declarations */
class SystemController;
class BasisController;

struct WriteDataTaskSettings {
  WriteDataTaskSettings() : 
  fname("data.hdf5"),
  activeOrbs({}),
  eriThreshold(-1){};
  REFLECTABLE(
    (std::string)fname,
    (std::vector<unsigned int>)activeOrbs,
    (double)eriThreshold
  )
};

/**
 * @class WriteDataTask WriteDataTask.h
 * @brief A task which write HF data to file.
 */
class WriteDataTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param system The system controller.
   */
  WriteDataTask(std::shared_ptr<SystemController> system);
  /**
   * @brief Default destructor.
   */
  virtual ~WriteDataTask() = default;
  /**
   * @brief Execute the task.
   */
  void run();
  /**
   * @brief The settings.
   */
  WriteDataTaskSettings settings;

 private:
  std::shared_ptr<SystemController> _system;
  std::shared_ptr<BasisController> _basis;
  std::shared_ptr<OrbitalController<RSCF>> _orbitals;
  std::shared_ptr<ElectronicStructure<RSCF>> _elstruct;

  std::unique_ptr<Eigen::MatrixXd> _coeffs;
  std::unique_ptr<RegularRankFourTensor<double>> _eris;

  unsigned int _nBasisFunc;   //number of basis functions
  unsigned int _nActOrbs;     //number of active orbitals to save
  unsigned int _nOccActOrbs;  //number of occupied active orbitals
  unsigned int _iActOrb;      //index of initial active orbital
  unsigned int _fActOrb;      //index of final active orbital

  void prepCoefficients(Eigen::VectorXd& coeffs);
  void prepERIS(Eigen::VectorXd& integrals, Eigen::VectorXd& indices);
};

} /* namespace Serenity */

#endif // SERENITY_WRITEDATATASK_H
