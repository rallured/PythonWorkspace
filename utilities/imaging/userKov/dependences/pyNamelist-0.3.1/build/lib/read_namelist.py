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
# Program: read_namelist
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
#   text_block = sting keyword containing a block of text known to hold a
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
#
#------------------------------------------------------------------------------

def read_namelist( namelist_name = '', text_block = '' ):

    import string
    import string_to_number

# Determine the souce of input: file or block of text. If both file and 
# text_block are specified default to the file and issue a warning. If a file,
# confirm that it is there.

    namelist_found = 'false'
    source = ''
    if namelist_name != '' and text_block == '':
	source = 'file'
	try:
	    namelist = open( namelist_name, 'r' )
	    namelist_found = 'true'
	except IOError:
	    print 'IOError: No such file or directory'
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

    print 'namelist_name =', namelist_name
    print 'source =', source

    namelist_dictionary = {}
    if namelist_found == 'true' or source == 'text':

# Stuff namelist_string from the appropriate source

	if namelist_found == 'true':
	    namelist_string = namelist.read()
	    namelist.close()
	else:
	    namelist_string = text_block

# The approach is to divide and conquer. Each namelist is assumed to be of the
# form ( header, namelist ) where header lines are delimited by ";" characters
# and the namelist is bounded by "$" characters. Finally, there may be a block
# of text following the last namelist in the file. This is saved in ( tailer ).
# The file is first parsed into these blocks and saved in a dictionary.

	namelist_count = 0
	while namelist_string != '':
	    namelist_count = namelist_count + 1
	    
# First peal off the header

	    dollar_index_1 = string.find( namelist_string, '$' )
	    header_text = namelist_string[ : dollar_index_1 ]
	    namelist_string =                                                 \
	       string.strip( namelist_string[ dollar_index_1 : ] )

# Extract the name of the namelist.

            nl_name = string.split( namelist_string )[0]
	    nl_name_len = len( nl_name )
	    nl_name = string.upper( nl_name[ 1 : ] )
	    namelist_string = namelist_string[ nl_name_len : ]

# Next, peal off the body of the namelist. Remove line feeds and make sure
# equal signs are distinct from both variable names and values.

            dollar_index_2 = string.find( namelist_string, '$' )
	    namelist_body = namelist_string[ : dollar_index_2 ]
	    namelist_body = string.join( string.split( namelist_body, 
	        '\012' ), ' ' )
	    namelist_body = string.join( string.split( namelist_body, '=' ),
		' = ' )
	    namelist_string = string.strip( 
		namelist_string[ dollar_index_2 : ] )

# Next, snip off the '$' or '$end' demarking the end of the namelist.

            nl_end = string.split( namelist_string )[0]
	    nl_end_len = len( nl_end )
	    namelist_string = namelist_string[ nl_end_len : ]

# Now stuff the pieces into the dictionary

	    namelist_dictionary[nl_name] = {}
	    namelist_dictionary[nl_name]['header'] = header_text
	    namelist_dictionary[nl_name]['variables'] = string.join( 
		string.split( namelist_body ) )

# Finally, check if there are any remaining '$' in the file. If not, what is 
# left is the tailer.

            if string.find( namelist_string, '$' ) == -1:
		namelist_dictionary['tailer'] = namelist_string
		namelist_string = ''

# Next extract the variable and their values and put them in the dictionary. 
# To do so continue the divide and conquer strategy. Here, "=" signs form the 
# delimiters between variable names and their values. Note that this method
# of parsing the files uniquely identifies variables, so that only one copy
# (the last occurance) of a variable is saved. We also clean up and convert
# the string values of the variables to python values

	for nl_name in namelist_dictionary.keys():
	    if nl_name == 'tailer': continue

# Convert the string into a list and initialize a dictionary for the namelist 
# and its variables.

	    variables_and_values = string.split( 
		namelist_dictionary[nl_name]['variables'] )
	    namelist_dictionary[nl_name]['variables'] = {}

# Make a list of the indices of the equal signs which form the delimiters 
# betweent the variable names and their values.

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

# 2) Replace commas with space and form a list of values

		variable_value = string.split( string.join( string.split( 
		    string.join( variable_value ), ',' ), ' ' ) )

# 3) Translate the values contained in the strings into Python numbers

		expanded_value = []
		for value in variable_value:
		    parsed_value = string_to_number.parse_number_string(value)
		    expanded_value = expanded_value + [ parsed_value['Value'] 
			] * parsed_value['Replications']

# 4) Save the list under the variable name. If the value is a scalar, save
#    the value itself rather than as a list containing one item

                if len( expanded_value ) == 1:
		    expanded_value = exapanded_value[0]
		namelist_dictionary[nl_name]['variables'][variable_name] =    \
		    expanded_value

    return namelist_dictionary
