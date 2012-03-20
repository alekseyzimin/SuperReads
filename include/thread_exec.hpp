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

/** Thread_exec is a helper class to spawn threads. To use, create a
 * class inheriting from thread_exec and implement the virtual method
 * start, which is the routine that you want to run in parallel. Upon
 * calling `exec(n)`, the `start(int)` method will be called `n`
 * times, each in a different thread. The parameter `id` to start will
 * have value `0` to `n-1`.
 * 
 * The `exec_join()` method does exec and then join for all the
 * threads. I.e., unlike `exec()` which returns immediately,
 * `exec_join()` returns when all the threads are done.
 *
 *  Example:
 *  ~~~{.cc}
 * class MyThread : public thread_exec {
 * public:
 *   virtual void start(int id) {
 *     std::cout << "I am thread " << id << " " << pthread_self() << "\n";
 *   }
 * };
 *
 * MyThread t;
 * t.exec_join(5);
 * ~~~
 *
 * The start method is called 5 times and 0 to 4 will be printed in
 * some random order with a different thread number. The output is
 * likely to be mixed between the different threads. See
 * [o_multiplexer](@ref jflib::o_multiplexer) to avoid this.
 */
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
  /** Method to execute in sub-thread. This method must be implement
   * in the sub-class.
   *
   * @param id Thread id. `id` will be in `[0;n-1]` if `n` threads are
   * started.
   */
  virtual void start(int id) = 0;
  /** Start threads and return immediately.
   * 
   * @param nb_threads Number of threads to start.
   */
  void exec(int nb_threads);
  /** Wait for previously started threads. */
  void join();
  /** Start threads and wait for them. 
   * @param nb_threads  Number of threads to start.
   */
  void exec_join(int nb_threads) {
    exec(nb_threads);
    join();
  }
};

#endif // __THREAD_EXEC_HPP__
