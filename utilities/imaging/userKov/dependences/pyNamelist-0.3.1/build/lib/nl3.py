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
import set_theoretic
import string
import string_to_number
import types, gzip

#-----------------------------------------------------------------------------
# Module: nl3
#
#   Contains routines to read and write namelists.
#
#       read, write
#
#------------------------------------------------------------------------------
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#------------------------------------------------------------------------------
# Program: read
#
#   Routine to read parse a FORTRAN style namelist and save them in a PYTHON
#   dictionary. The source of the namelist can either be a standard ascii
#   namelist file or a block of text form a source other than a file 
#   containing the elements of a namelist
#
#------------------------------------------------------------------------------
# Inputs:
#
#   namelist_name = string keyword containing the name of FORTRAN style 
#       namelist 
#   text_block = string keyword containing a block of text known to hold a
#       namelist
#
#------------------------------------------------------------------------------
# Outputs:
#
#   Python dictionary containing the contents of the namelist
#
#------------------------------------------------------------------------------
# History (dd-mm-yy)
# 
#   07-09-00  Makowski's first Python Program (that can really do something
#             useful)
#   11-09-00  Finished basic code. Sure there are some gotcha's out there, but
#             this works in the non-trivial test cases I've tried.
#   26-09-00  Generalized code to handle namelist parsing of both files and 
#             strings. Converted namelist_name from a normal argument to a
#             keyword argument and added a second keyword argument holding a
#             string containing a namelist. The change was instagated to be
#             able to handle g-files which are comprised of both undeclared
#             quantities followed by several namelists. (mam)
#   03-10-00  After a day of frustration decided that the best way to parse 
#             the number strings was to do so as part of this routine, rather 
#             than after all the strings were read. The biggest problem was
#             with the replicator strings.
#   05-10-00  Tidied up a few things:
#             - singletons saved as values rather than as a list of length 1
#   29-10-01  Added tag, _nl_sequence, so that files containing multiple 
#             namelist can write them in the same order that they are read.
#   10-07-02  Added coding to detect F90 namelist delimiters: & and /
#   12-08-03  Slight modifications made to handle Lahey compiler namelists
#
#------------------------------------------------------------------------------

def read( namelist_name = '', text_block = '' ):

  """
  Reads a namelist and converts into a dictionary which can then be modified
  by the user. Parameters:

      namelist_name = string containing file name to parse (default)
      text_block = string (block of text) to parse as a namelist

  The options are mutually exclusive; specifying both can lead to unexpected
  results. Usage:

      import nl3
      k1 = nl3.read('k101555.02750')

  """
    
  # Determine the souce of input: file or block of text. If both file and 
  # text_block are specified default to the file and issue a warning. If a
  # file, confirm that it is there.

  namelist_found = 'false'
  source = ''
  if namelist_name != '' and text_block == '':
    source = 'file'
    if namelist_name[-3:].lower() == '.gz':
      try:
	namelist = gzip.open( namelist_name, 'r' )
	namelist_found = 'true'
      except IOError:
	print 'IOError in nl3: No such file or directory: %s' % namelist_name
    else:
      try:
	namelist = open( namelist_name, 'r' )
	namelist_found = 'true'
      except IOError:
	print 'IOError in nl3: No such file or directory: %s' % namelist_name
  elif text_block != '' and namelist_name == '':
    source = 'text'
  elif text_block == '' and namelist_name == '':
    print 'Warning: No namelist source specified'
    print 'Returning empty dictionary.'
    print ''
  else:
    print 'Warning: Two namelist sources specified.'
    print 'Defaulting to file input.'
    print ''
    source = 'file'

  namelist_dictionary = {}
  if namelist_found == 'true' or source == 'text':

    # Stuff namelist_string from the appropriate source

    if namelist_found == 'true':
      namelist_string = namelist.read()
      namelist.close()
    else:
      namelist_string = text_block

    # The approach is to divide and conquer. Each namelist is assumed to be of
    # the form ( header, namelist ) where header lines are delimited by ";" 
    # characters and the namelist is bounded by "$" characters or "&" and '/' 
    # pairs. Finally, there may be a block of text following the last namelist 
    # in the file. This is saved in ( tailer ). The file is first parsed into 
    # these blocks and saved in a dictionary.

    namelist_count = 0

    if ( string.find( namelist_string, '$' ) != -1 ):
      begin_delimiter = '$'
      end_delimiter = '$'
      namelist_dictionary['_delimiter'] = '$'
    else:
      begin_delimiter = '&'
      end_delimiter = '/\012'
      namelist_dictionary['_delimiter'] = '&'

    while namelist_string != '':
	    
      # First peal off the header

      delimiter_index_1 = string.find( namelist_string, begin_delimiter )
      header_text = string.strip( namelist_string[ : delimiter_index_1 ] )
      namelist_string = namelist_string[ delimiter_index_1 : ]

      # Extract the name of the namelist.

      nl_name = string.split( namelist_string )[0]
      nl_name_len = len( nl_name )
      nl_name = string.upper( nl_name[ 1 : ] )
      namelist_string = namelist_string[ nl_name_len : ]

      # Next, peal off the body of the namelist. Remove line feeds and commas, 
      # and make sure equal signs are distinct from both variable names and 
      # values.

      delimiter_index_2 = string.find( namelist_string, end_delimiter )
      namelist_body = namelist_string[ : delimiter_index_2 ]
      namelist_body = string.join( string.split( namelist_body, '\012' ), ' ' )
      namelist_body = string.join( string.split( namelist_body, ',' ), ' ' )
      namelist_body = string.join( string.split( namelist_body, '=' ), ' = ' )
      namelist_string = namelist_string[ delimiter_index_2 : ]

      # Next, snip off the '$' or '$end' or '/' demarking the end of the namelist.

      if ( end_delimiter == '/\012' ):
 	namelist_string = namelist_string[3:]
      elif ( string.find( namelist_string, '$END' ) != -1 ):
	first_char = string.find( namelist_string, '$END' )
	namelist_string = namelist_string[ first_char + 5 : ]
      elif ( string.find( namelist_string, '$end' ) != -1 ):
	first_char = string.find( namelist_string, '$end' )
	namelist_string = namelist_string[first_char + 5: ]
      elif ( string.find( namelist_string, '$' ) != -1 ):
	first_char = string.find( namelist_string, '$' )
	namelist_string = namelist_string[ first_char + 2 : ]
      
      # Now stuff the pieces into the dictionary

      namelist_dictionary[nl_name] = {}
      namelist_dictionary[nl_name]['_nl_header'] = header_text
      namelist_dictionary[nl_name]['_nl_sequence'] = namelist_count
      namelist_dictionary[nl_name]['variables'] = string.join( 
	string.split( namelist_body ) )

      # Finally, check if there are any remaining '$' in the file. If not, 
      # what is left is the tailer.

      if string.find( namelist_string, begin_delimiter ) == -1:
	namelist_dictionary['tailer'] = namelist_string
	namelist_string = ''

      namelist_count = namelist_count + 1

    # Next extract the variable and their values and put them in the 
    # dictionary. To do so continue the divide and conquer strategy. Here, 
    # "=" signs form the delimiters between variable names and their values.
    # Note that this method of parsing the files uniquely identifies 
    # variables, so that only one copy (the last occurance) of a variable is
    #  saved. We also clean up and convert the string values of the variables
    # to python values

    for nl_name in set_theoretic.complement( namelist_dictionary.keys(), ['_delimiter'] ):
      if nl_name == 'tailer': continue

      # Convert the string into a list and initialize a dictionary for the 
      # namelist and its variables.

      variables_and_values = string.split( 
	namelist_dictionary[nl_name]['variables'] )
      del namelist_dictionary[nl_name]['variables']

      # Make a list of the indices of the equal signs which form the  
      # delimiters betweent the variable names and their values.

      equal_indices = []
      for i in range( len(variables_and_values) ):
	if variables_and_values[i] == '=': equal_indices.append(i)
      equal_indices.append( len(variables_and_values) + 1 )

      # Operate on the variables_and_values list:

      for j in range( len(equal_indices) - 1 ):

	# 1) Extract a variable/value pair

       	variable_name = variables_and_values[ equal_indices[j] - 1 ]
	variable_value = variables_and_values[ equal_indices[j] + 1 :
	  equal_indices[ j + 1 ] - 1 ]

	# 2) Treat special case resulting from strings containing whitespace

	if variable_value[0][0] == "'":
	  variable_value = [ "%s" % string.join(variable_value) ]

	# 3) Translate the values contained in the strings into Python numbers

	expanded_value = []
	for value in variable_value:
	  parsed_value = string_to_number.parse_number_string(value)
	  expanded_value = expanded_value + [ parsed_value['Value'] 
	    ] * parsed_value['Replications']

	# 4) Save the list under the variable name. If the value is a scalar,
        #    save the value itself rather than as a list containing one item

	if len( expanded_value ) == 1:
	  expanded_value = expanded_value[0]
	namelist_dictionary[nl_name][variable_name] = expanded_value

  return namelist_dictionary

#------------------------------------------------------------------------------
# Routine: write
#
#   Given a namelist dictionary, writes a namelist file
#
#------------------------------------------------------------------------------
# Author:
#
#   Mike Makowski
#
#------------------------------------------------------------------------------
# Inputs:
#
#   nl_dict = namelist dictionary as generated by read_namelist
#   nl_file_name = name of file to write namelist to
#   max_line_length = maximum line length of output file (default = 80)
#
#------------------------------------------------------------------------------
# Outputs:
#
#   File named, nl_file_name, with contents on nl_dict
#
#------------------------------------------------------------------------------
# History (dd-mm-yy):
#
#   17-10-01  Adopted from write method in k_file. This is a more general 
#             routine.
#
#------------------------------------------------------------------------------
# Notes:
#
#   _nl_header and _nl_sequence are internal parameters which are maintained
#   as  dictionary entries along with the namelist variables. Exception coding
#   treats these as needed.
#
#------------------------------------------------------------------------------

def write( nl_dict, nl_file_name, max_line_length=80 ):

  """
  Writes a namelist dictionary to a file. Parameters:

     nl_dict = namelist dictionary to write to file
     nl_file_name = string containing the name of the file to write to
     max_line_length = optional parameter spedicifying the line length
                       (default = 80)

  Usage:

     import nl3
     nl3.write( k1, 'k1_out.nl', max_line_length = 128 )
  """

  #----------------------------------------------------------------------------

  def fort_bool(value):

    if value[0] == 0 and value[1] == 'Boolean':
      return 'f'
    else:
      return 't'
  
  #----------------------------------------------------------------------------

  nl_out = open( nl_file_name, 'w' )

  # First order the namelists in the sequence they were writen. Next write 
  # the associated header. Delimit the namelist with "$"s and write the 
  # variables and their values separated by and equal, "=', sign. Multiple 
  # line outputs must be indented for clarity. No character should appear in 
  # column 1. Strings and Booleans need to be treated as exception cases.

  namelists = set_theoretic.complement( nl_dict.keys(), [ 'tailer', '_delimiter' ] )
  sequence = []
  for namelist in namelists: 
    sequence.append( nl_dict[namelist]['_nl_sequence'] )
  ordered_namelists = [ '' ]*len(sequence)
  for index in range(len(sequence)):
    sequence_index = sequence.index(index)
    ordered_namelists[index] = namelists[sequence_index]

  for namelist in ordered_namelists:
    header_text = string.rstrip( nl_dict[namelist]['_nl_header'] )
    if header_text != '': nl_out.write( header_text + '\n' )
    nl_out.write( ' ' + nl_dict['_delimiter'] + namelist + '\n' )
    variables = set_theoretic.complement(
      nl_dict[namelist].keys(), [ '_nl_header', '_nl_sequence' ] ) 
    for variable in variables:
      value = nl_dict[namelist][variable]

      # Convert everything to a list if it isn't already. However, there
      # are two exceptional cases:
      #   - a single instance of a string as in "'/link/efit '"; len = 13
      #   - a single boolean as in (1, 'Boolean'); len = 2

      try:
	vector_length = len(value)
      except TypeError:
	value = [ value ]
	vector_length = 1
      else:
	if vector_length == 2 and value[1] == 'Boolean':
	  value = [ value ]
	  vector_length = 1
	elif value[0] == "'":
	  value = [ value ]
	  vector_length = 1

      # Now get the proper string representation of each of the types

      try:
	chk = value[0][1]
      except ( TypeError, IndexError ):                 # It's a number
	str_value = map( repr, value )
      else:
	if types.StringType == type(value[0]):          # It's a string
	  str_value = value
	else:                                           # It's a boolean
	  str_value = map( fort_bool, value )

      # Finally write the variable and it's value out.
	    
      len_substrings = map( len, str_value )
      line_count = 0
      line = [ ' ' + variable + ' =' ]
      line_length = len(line[0]) + 1
      while str_value != []:
	if line_length + len_substrings[0] + 1 <= max_line_length:
	  len_one = len_substrings.pop(0)
	  line_length = line_length + len_one + 1
	  str_value_one = str_value.pop(0)
	  line.append( str_value_one )
	else:
	  nl_out.write( string.join(line) + '\n' )
	  line_count = line_count + 1
	  line = [ ' ' ]
	  line_length = 2

      if line != [ ' ' ]:
	nl_out.write( string.join(line) + '\n' )
    
    if ( nl_dict['_delimiter'] =='$' ):
      nl_out.write( ' $END\n' )
    else:
      nl_out.write( ' /\n' )
      
  # Finally write the tailer

  nl_out.write( nl_dict['tailer'] )

  nl_out.close()

#------------------------------------------------------------------------------
