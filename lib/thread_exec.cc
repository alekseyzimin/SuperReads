/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <thread_exec.hpp>

void thread_exec::exec(int nb_threads) {
  struct thread_info info;
  infos.reserve(nb_threads);

  for(int i = 0; i < nb_threads; i++) {
    info.id   = i;
    info.self = this;
    infos.push_back(info);
    int err = pthread_create(&infos[i].thid, NULL, start_routine, &infos[i]);
    if(err)
      eraise(Error) << "Can't create thread " << i << err::str(i);
  }
}

void thread_exec::join() {
  typedef std::vector<std::pair<int,int> > err_info;
  err_info errors;
  for(unsigned int i = 0; i < infos.size(); i++) {
    int err = pthread_join(infos[i].thid, NULL);
    if(err) 
      errors.push_back(std::make_pair(i, err));
  }
  if(!errors.empty()) {
    std::ostringstream err;
    for(err_info::const_iterator it = errors.begin(); it != errors.end(); ++it)
      err << "Can't join thread " << it->first << err::str(it->second);
    eraise(Error) << "The following errors occured:\n" 
                  << err.str();
  }
}

void *thread_exec::start_routine(void *_info) {
  struct thread_info *info = (struct thread_info *)_info;
  info->self->start(info->id);
  return 0;
}
