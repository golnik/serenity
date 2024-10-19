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

namespace Serenity {
/* Forward declarations */
class SystemController;

struct WriteDataTaskSettings {
  WriteDataTaskSettings() : dummy(false){};
  REFLECTABLE((bool)dummy)
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
};

} /* namespace Serenity */

#endif // SERENITY_WRITEDATATASK_H
