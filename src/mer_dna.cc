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


#include <jellyfish/mer_dna.hpp>
#include <iostream>

const char mer_dna_ns::dna_codes::rev_codes[4] = { 'A', 'C', 'G', 'T' };

#define R -1
#define I -2
#define O -3
#define A 0
#define C 1
#define G 2
#define T 3

const int mer_dna_ns::dna_codes::codes[256] = {
  O, O, O, O, O, O, O, O, O, O, I, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, R, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O, 
  O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O, 
  O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O, 
  O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O
};

#undef R
#undef I
#undef O
#undef A
#undef C
#undef G
#undef T

const char* const mer_dna_ns::error_different_k = "Length of k-mers are different";
const char* const mer_dna_ns::error_short_string = "Input string is to short";
