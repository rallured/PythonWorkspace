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
"""
MODULE namelist : Read and write FORTRAN namelist files.
"""
import set_theoretic
import string,re
import string_to_number
import types, gzip

#-----------------------------------------------------------------------------
# Module: namelist
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
#    2-20-06  Modification by T. Osborne to keep variables in order from read to
#             write and deal with indexed array input. The later is treated as a
#             separate namelist variable e.g. c(3) is separate from c. 
#
#------------------------------------------------------------------------------

def read( namelist_name = '', text_block = '' , combine_arrays = 1 ):

  """
  Reads a namelist and converts into a dictionary which can then be modified
  by the user. Parameters:

      namelist_name = string containing file name to parse (default)
      text_block = string (block of text) to parse as a namelist
      combine_arrays = If 1 namelist elements referencing the same array
         are combined togeather into a single entry. e.g. if x(1) and x(4)
         appear in the namelist, and entry x is created with length 4. Values
         which are not set ( x(2) and x(3)) are set to 0. If combine_arrays=0
         then x(1) and x(4) will be separate in the namelist dictionary with
         'x(1)' and 'x(4)' used as the keys. Note that this is only for 1D arrays
         since it is impossible to tell the shape of an n>1D from the data in
         the namelist
         
  The options are mutually exclusive; specifying both can lead to unexpected
  results. Usage:

      import namelist
      k1 = namelist.read('k101555.02750')

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
	print 'IOError in namelist: No such file or directory: %s' % namelist_name
    else:
      try:
	namelist = open( namelist_name, 'r' )
	namelist_found = 'true'
      except IOError:
	print 'IOError in namelist: No such file or directory: %s' % namelist_name
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

    begin_rea = re.compile("^( )*&(\w)*\s",re.MULTILINE)
    begin_red = re.compile("^( )*\$(\w)*\s",re.MULTILINE)
    if begin_rea.search(namelist_string):
      begin_re = begin_rea
      end_re = re.compile("/",  re.MULTILINE)
      begin_delimiter = '&'
      end_delimiter = '/\012'
    elif begin_red.search(namelist_string): 
      begin_re = begin_red
      end_re   = re.compile("^( )*\$",  re.MULTILINE)
      begin_delimiter = '$'
      end_delimiter   = '$'
    else:
      raise IOError('NO NAMELISTS SECTIONS FOUND')

    namelist_dictionary['_delimiter'] = begin_delimiter

    # Used to strip blanks from indexed entries
    pindex_re = re.compile("\([^\)]+\)")

    # Used to strip comments imbedded within namelist blocks
    comment_re = re.compile("!(.)*\n",re.MULTILINE)

    namelist_string = namelist_string.replace( '\"', '\'' )
    # Used to search for quotes
    quote_re = re.compile("\'",re.MULTILINE)
    
    # Used to strip line feeds from sting variables.
    string_re = re.compile('\'(.|\n)*?\'')
    
    while namelist_string != '':
	    
      # Span of delimiter + namelistname
      delimiter_index_1, name_index_2 = begin_re.search( namelist_string ).span()

      # Header text is before delimiter
      header_text = string.strip( namelist_string[ : delimiter_index_1 ] )

      # Namelist name
      nl_name = namelist_string[delimiter_index_1:name_index_2].strip()[1:].upper()

      # namelist boday is after delimiter + namelist name
      namelist_string = namelist_string[ name_index_2 : ]

      # Next, peal off the body of the namelist. Remove line feeds and commas, 
      # and make sure equal signs are distinct from both variable names and 
      # values.

      # Find location on of end delimiter. Must take into account that the
      # end delimiter could occur in within a string inside the namnelist
      # Position of fist end delimiter
      k = end_re.search( namelist_string ).end()
      # Loop until there are an even number of quote marks before the end terminator
      while len( quote_re.findall( namelist_string[:k] ) )%2 :
        k =  end_re.search( namelist_string[k:] ).end() + k
      # Index of end delimiter
      delimiter_index_2 = k - 1

      # Namelist_body
      namelist_body   = namelist_string[ : delimiter_index_2 ]
      
      # Rest of namelist_string buffer
      namelist_string = namelist_string[ delimiter_index_2 : ]

      # Strip off from the end delimiter to the next line
      k = string.find( namelist_string, '\n' )
      namelist_string = namelist_string[k+1:]

      # Strip comments from namelist body
      namelist_body = comment_re.sub('\n',namelist_body)

      # The next bit is to deal with the possibility that some
      # of the strings have imbedded line feeds, this seems to be
      # possible for Layhey Fujitsu compiler namelist output.
      k = 0
      strings = []
      while 1:
        v = string_re.search( namelist_body[k:] )
        if v is None:break
        s = list(v.span())
        s[0] = s[0]+k
        s[1] = s[1]+k
        strings.append( namelist_body[s[0]:s[1]] )
        k = s[1]
      for s in strings:
        if '\n' in s or ' ' in s:
          v=s.replace('\n','').replace(' ','!')
          namelist_body = namelist_body.replace( s, v )
        
      # Strip linefeeds from namelist_body
      namelist_body = string.join( string.split( namelist_body, '\012' ), ' ' )
      
      # Strip blanks from indexed entries, e.g abc( 3, 2)=5
      # and replace , in multi indices with %
      m = 0
      lnl = len(namelist_body)
      while 1:
        kk = pindex_re.search( namelist_body[m:] )
        if kk is None or m >= lnl:break
        ki = kk.span()
        s = namelist_body[ ki[0]+m : ki[1]+m ].replace(' ','').replace(',','%')
        namelist_body = namelist_body[ :ki[0] + m ] + s + namelist_body[ ki[1] + m: ]
        m = m + ki[1]
      # Strip commas from namelist_body
      namelist_body = string.join( string.split( namelist_body, ',' ), ' ' )
      # Surround = with blanks in namelist_body
      namelist_body = string.join( string.split( namelist_body, '=' ), ' = ' )
      # Put comma back into multi indicies
      namelist_body = namelist_body.replace('%',',')
      
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
      namelist_dictionary[nl_name]['_nl_order']=[]

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
	  variable_value = [ ("%s" % string.join(variable_value,'')).replace('!',' ') ]

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
	namelist_dictionary[nl_name]['_nl_order'].append( variable_name )
      var = namelist_dictionary[nl_name]['_nl_order'][:]
      var.reverse()
      var = set_theoretic.unique(var)
      var.reverse()
      namelist_dictionary[nl_name]['_nl_order']=var
        
  # Combine data for a given array into a single element
  # This can only be done for 1D arrays
  if combine_arrays:
    for nl in namelist_dictionary.keys():
      if nl in ['tailer','_delimiter']:continue
      anames  = []  # List of array names
      karrays = {}  # Dict with key=array names and value=[ variables of that array ]

      # Get array names
      for v in namelist_dictionary[nl]['_nl_order']:
        if '(' in v and ')' in v and not ',' in v:
          aname = v.split('(')[0]
          if aname not in anames:
            anames.append(aname)
            karrays[aname] = []
      if not anames:continue

      # Get list of variables referencing a given array
      for v in namelist_dictionary[nl]['_nl_order']:
        if '(' in v and ')' in v and not ',' in v:
          aname = v.split('(')[0]
          karrays[aname].append(v)
          continue
        if v in anames:
          karrays[v].append(v)

      arrays  = {}
      # Combine values into a single array
      for a in karrays.keys():

        # Determine number of elements in the array
        n=1
        for v in karrays[a]:
          if '(' in v:
            i0 = int(v.split('(')[1].split(')')[0])
          else:
            i0 = 1
          y = namelist_dictionary[nl][v]
          if ( type(y) is types.ListType ):
            n = max( n, i0 + len(y) - 1 )
            atype = type( y[0] )
          else:
            n = max( n, i0 )
            atype = type( y )

        # Create a list to hold all the defined elements
        if atype is types.FloatType:
          arrays[a] = [ 0.0 ]*n
        elif atype is types.IntType:
          arrays[a] = [ 0 ]*n
        elif atype is types.StringType:
          arrays[a] = [ '' ]*n
        elif atype is types.TupleType:
          arrays[a] = [ (0, 'Boolean') ]*n

        # Assign the values to the list
        for v in karrays[a]:
          if '(' in v:
            i0 = int(v.split('(')[1].split(')')[0])
          else:
            i0 = 1
          y = namelist_dictionary[nl][v]
          if ( type(y) is types.ListType ):
            arrays[a][i0-1:i0-1+len(y)] = y[:]
          else:
            arrays[a][i0-1] = y

        # Delete explicit index references and replace with array type entry
        for v in karrays[a]:
          del namelist_dictionary[nl][v]
          k = namelist_dictionary[nl]['_nl_order'].index(v)
          del namelist_dictionary[nl]['_nl_order'][k]

      # Order combined arrays at end of namelist
      kn = arrays.keys()
      kn.sort()
      for k in kn:
        namelist_dictionary[nl][k] = arrays[k]
        namelist_dictionary[nl]['_nl_order'].append(k)

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

     import namelist
     namelist.write( k1, 'k1_out.nl', max_line_length = 128 )
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
    if '_nl_order' in nl_dict[namelist].keys():
      variables = nl_dict[namelist]['_nl_order']
    else:
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
