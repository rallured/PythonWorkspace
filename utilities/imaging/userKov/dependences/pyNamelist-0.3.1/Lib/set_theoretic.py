#
# Copyright (C) 2004 Mike Makowski (Makowski@fusion.gat.com)
#
# This program is part of the pyNamelist package
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
def intersect( *args ):

  res = []
  for x in args[0]:
    for other in args[1:]:
      if x not in other: break
      else:
	res.append(x)

  return res



def union( *args ):

  res = []
  for seq in args:
    for x in seq:
      if not x in res:
	res.append(x)

  return res



def complement( universe, subset ):

  """
  Elements of universe not found in a_set (the complement of subset relative to
  universe)
  """

  comp = []
  for element in universe:
    if element not in subset:
      comp.append(element)

  return comp



def unique( set ):

  """
  Returns a list of unique elements in the supplied set
  """

  unique_elements = []
  for element in set:
    if element not in unique_elements:
      unique_elements.append(element)

  return unique_elements
