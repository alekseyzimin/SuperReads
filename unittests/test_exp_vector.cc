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


#include <gtest/gtest.h>
#include <exp_vector.hpp>

TEST(ExpVector, Subscript) {
  exp_vector<int> ev;
  exp_vector<int, 5> evd;

  ev[10]  = 5;
  evd[10] = 10;
  for(int i = 0; i < 10; ++i) {
    EXPECT_EQ(0, ev[i]);
    EXPECT_EQ(5, evd[i]);
  }
  EXPECT_EQ(5, ev[10]);
  EXPECT_EQ(10, evd[10]);
}
