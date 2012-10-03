/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _THREAD_POOL_H_
#define _THREAD_POOL_H_

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>

#include <stdexcept>

// Modifications
// - Changed data_available to a bool, because it really is a boolean
// - Moved the struct arguments inside the class, where it belongs
// - Use work_finished and pthread_cond_broadcast to signal that there is no more work
// - Keep track whether a thread was actually started
// - Call release_workers from the destructor
// - deleted default copy constructor

//#define DEBUG 1
#undef DEBUG

template <typename T, typename U>
class thread_pool {
  typedef U (*function_ptr)(T);
  struct arguments {
    T *f_arg;
    U *retVal;
  };
  struct thread_info {
    pthread_t id;
    bool      started;
  };

  int              num_threads;
  thread_info     *thread_id;
  arguments        arg; 
  volatile bool    data_available;
  volatile bool    work_finished;
  volatile int     done;
  struct timespec  t;
  function_ptr     F;

  pthread_mutex_t lock1, lock2, lock3;
  pthread_cond_t  cond_avl, cond_acq;

  static void* __startThread(void* self) {
    ((thread_pool*)self)->worker();
    return 0;
  }

  void worker() {
    U *retValue;
    T  argumentLocal;

#ifdef DEBUG
    printf("Thread %ld is alive!!! %p\n", pthread_self(), this);
#endif
    while(true){
#ifdef DEBUG
      printf("Worker waiting for data %d\n",data_available);
#endif
      pthread_mutex_lock(&lock1);
      if(work_finished)
        goto finished;
      while(!data_available) {
        pthread_cond_wait(&cond_avl,&lock1);
        if(work_finished)
          goto finished;
      }
#ifdef DEBUG
      printf("Worker copying data\n");
#endif
      argumentLocal  = *arg.f_arg;
      retValue       = arg.retVal;
      pthread_mutex_lock(&lock3);
      data_available = false;
      pthread_cond_signal(&cond_acq);
      pthread_mutex_unlock(&lock3);
      pthread_mutex_unlock(&lock1);

      /* call the function */
      *retValue = (*F)(argumentLocal);

      pthread_mutex_lock(&lock2);
      done++;
#ifdef DEBUG
      printf("Thread %d: done = %d\n",(int)pthread_self(),done);
#endif
      pthread_mutex_unlock(&lock2);
    }

  finished:
    pthread_mutex_unlock(&lock1);
  }

public:

  thread_pool(int num_threads_, function_ptr F_) :
    num_threads(num_threads_),
    thread_id(new thread_info[num_threads]),
    arg(),
    data_available(false),
    work_finished(false),
    done(0),
    t({0, 5000}),
    F(F_)
  {
#ifdef DEBUG
    printf("Pool %p\n", this);
#endif
    pthread_mutex_init(&lock1,NULL);
    pthread_mutex_init(&lock2,NULL);
    pthread_mutex_init(&lock3,NULL);
    pthread_cond_init(&cond_avl,NULL);
    pthread_cond_init(&cond_acq,NULL);

    for(int i = 0; i < num_threads; ++i)
      thread_id[i].started = false;

    for(int i=0;i<num_threads;i++) {
      if(pthread_create(&thread_id[i].id, NULL, __startThread, this))
        throw std::runtime_error("Failed to start a thread in the pool");
      thread_id[i].started = true;
    }
  }

  thread_pool(const thread_pool& rhs) = delete;
  thread_pool& operator=(const thread_pool& rhs) = delete;

  ~thread_pool(void){
    release_workers();

    delete [] thread_id;
    pthread_mutex_destroy(&lock1);
    pthread_mutex_destroy(&lock2);
    pthread_mutex_destroy(&lock3);
    pthread_cond_destroy(&cond_acq);
    pthread_cond_destroy(&cond_avl);
  }

  int jobs_done(void){ return done; }

  void release_workers(void){
    pthread_mutex_lock(&lock1);
    work_finished = true;
    pthread_cond_broadcast(&cond_avl);
    pthread_mutex_unlock(&lock1);

    for(int i=0;i<num_threads;i++) {
      if(thread_id[i].started) {
        pthread_join(thread_id[i].id, NULL);
        thread_id[i].started = false;
      }
    }
  }

  void submit_job(T* arg_, U* ret){
#ifdef DEBUG
    printf("Waking worker %d\n",data_available);
#endif
    pthread_mutex_lock(&lock1);
    arg.f_arg      = arg_;
    arg.retVal     = ret;
    data_available = true;
    pthread_cond_signal(&cond_avl);
    pthread_mutex_unlock(&lock1);
#ifdef DEBUG
    printf("Signaled worker %d\n",data_available);
#endif
    /*wait until thread receives the data*/
    pthread_mutex_lock(&lock3);
    while(data_available)
      pthread_cond_wait(&cond_acq,&lock3);
    pthread_mutex_unlock(&lock3);
#ifdef DEBUG
    printf("Worker received data %d\n",data_available);
#endif

  }
};

#endif /* _THREAD_POOL_H_ */
