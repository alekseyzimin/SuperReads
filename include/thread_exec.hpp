/*  This file is part of Jflib.

    Jflib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jflib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jflib.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __JFLIB_THREAD_EXEC_HPP__
#define __JFLIB_THREAD_EXEC_HPP__

#include <pthread.h>
#include <vector>
#include <exception>
#include <stdexcept>
#include <string>
#include <string.h>
#include <err.hpp>

class thread_exec {
  struct thread_info {
    int          id;
    pthread_t    thid;
    thread_exec *self;
  };
  static void *start_routine(void *);
  std::vector<struct thread_info> infos;

public:
  define_error_class(Error);
  thread_exec() {}
  virtual ~thread_exec() {}
  virtual void start(int id) = 0;
  void exec(int nb_threads);
  void join();
  void exec_join(int nb_threads) {
    exec(nb_threads);
    join();
  }
};

#endif // __THREAD_EXEC_HPP__
