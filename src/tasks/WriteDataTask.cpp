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

namespace Serenity {

WriteDataTask::WriteDataTask(std::shared_ptr<SystemController> system) : _system(system) {
}

void WriteDataTask::run() {
  printSubSectionTitle("Saving Hartree-Fock data to hdf5 file");
  OutputControl::n.printf(
      "    This module will save HF energies, vectors, and molecular integrals to hdf5 file.\n\n");

  return;
}

} /* namespace Serenity */
