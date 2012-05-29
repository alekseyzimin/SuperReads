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
$shootingMethodSuperReadGroupFile = $ARGV[0];
$origSuperReadGroupFile = $ARGV[1];

open (FILE, $shootingMethodSuperReadGroupFile);
while ($line = <FILE>) {
    print $line;
    ($readName) = ($line =~ /^readName = (\d+)/);
    $alreadyOutput{$readName} = 1;
    $line = <FILE>;
    print $line;
}
close (FILE);

open (FILE, $origSuperReadGroupFile);
while ($line = <FILE>) {
    ($readName) = ($line =~ /^readName = (\d+)/);
    print $line unless ($alreadyOutput{$readName});
    $line = <FILE>;
    print $line unless ($alreadyOutput{$readName});
}
close (FILE);

