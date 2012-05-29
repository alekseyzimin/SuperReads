# SuperRead pipeline
# Copyright (C) 2012  Genome group at University of Maryland.
# 
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



#!/usr/bin/perl
# This exec expects 1 arg: the threshold number of occurrences which should
# be the maximum on the list; an arg of 5 would list all super-reads with
# a count of 5 or less
# Input on STDIN and should be one line for each super-read:
#  The number of times it occurs followed by the name of the super-read
$maxCount = $ARGV[0];
while ($line = <STDIN>) {
    chomp ($line);
    ($curCount, $superRead) = ($line =~ /^(\d+)\s+(\S.*\S)\s*$/);
    if ($curCount <= $maxCount) {
	print "$superRead\n"; }
}
