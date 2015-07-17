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
#------------------------------------------------------------------------------
# Program: scan_dictionary
#
#   Routine that recursively scans a dictionary and prints either the 
#   dictionary names and optionally thier associated values. Based on an 
#   analogous IDL (struct_show) routine by Schachter which does the same for 
#   IDL structures.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   dictionary_entry = a python dictionary
#   print_values = switch indicating whether to print values of dictionary
#     entries, or only dictionary entry names (default is NOT to print values)
#   offset = string parameter used by the routine to indent names in proportion
#     to their depth in the dictionary
#
#------------------------------------------------------------------------------
# Outputs:
#
#   A listing of the name of the dictionary entry and optionally it's 
#   associated value
#
#------------------------------------------------------------------------------
# History (dd-mm-yy)
# 
#   05-10-00  Was remnant of the string_to_number module but obsolesced. 
#             Extracted basic recursive routine and expanded options. (mam)
#   13-10-00  Added sorting to dictionary keys so they wouldn't come out in
#             random order.
#
#------------------------------------------------------------------------------

def scan_dictionary( dictionary_entry, print_values = 'false', offset = '' ):

    if print_values != 'false':
	offset_increment = ''
    else:
	offset_increment = '    '

    try:
        keys = dictionary_entry.keys()
    except AttributeError:
	if print_values != 'false':
	    print dictionary_entry
	    print ' '
    else:
	keys.sort(case_independent_sort)
	for key in keys:
	    print offset, key
	    scan_dictionary( dictionary_entry[key], print_values, 
	        offset = offset + offset_increment )

    return


#------------------------------------------------------------------------------
# Program: case_independent_sort
#
#   Used together with the built-in function, sort, to generate a case
#   independent sort. Program is from Learning Python, by M Lutz and D Ascher,
#   p. 246
#
#------------------------------------------------------------------------------
# Inputs:
#
#   something = list to be sorted
#
#------------------------------------------------------------------------------
# Outputs:
#
#   other = input list sorted in a case independent way
#
#------------------------------------------------------------------------------
# History:
#
#   13-10-00  Entered program
#
#------------------------------------------------------------------------------
# Notes:
#
#   It seems that you have to know how the sort function work to understand the
#   interaction of this routine and the sort routine itself. Still don't quite
#   get it.
#
#------------------------------------------------------------------------------

def case_independent_sort( something, other ):

    import string

    something = string.lower(something)
    other = string.lower(other)

    return cmp( something, other )
