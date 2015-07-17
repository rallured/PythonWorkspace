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
# Class: boolean
#
#   Python routine to create a boolean class object. Likely not much use other
#   than to demonstrate a non-trivial example of a class definition. Two 
#   methods are implemented, that of setting the value of the boolean object
#   (to its python value, either 0 or 1) and another to get the boolean value.
#   Various representions of a boolean value are supported. These include
#
#     1) Python: 0 = false, 1 = true
#     2) String: 'true' or 'false'
#     3) Dotted-Fortran: '.TRUE.' or '.FALSE.' 
#     4) Fortran: 'T' or 'F'
#     5) Switch: 'on' = true, 'off' = false
#     6) Answer: 'yes' = true, 'no' = false
#
#   Truth values are case insensitive.
#
#------------------------------------------------------------------------------
# Creating an instance:
#
#   import boolean
#   logical = boolean.boolean()
#   logical.set_value('.true.')
#
#------------------------------------------------------------------------------
# History (dd-mm-yy):
#
#   05-10-00  Began to write routine (mam)
#   06-10-00  Finished debugging. Overloaded the arithmatical logical operators
#             'and' = &, 'or' = |, 'xor' = ^, 'not' methods as well as the 
#             multiplication operator. (mam)
#
#------------------------------------------------------------------------------

class boolean:

    """
    Creates a boolean class object. Two methods are implemented

      1) set_value(truth_value)
      2) get_value(truth_value_type)

    Various representions of a boolean value are supported. These include

                              Truth Values
       Truth Value        --------------------
          Type              True      False
       ---------------------------------------
        'Python'              1         0
        'String'            'true'   'false'
        'Dotted_Fortran'   '.TRUE.' '.FALSE.' 
        'Fortran'            'T'       'F'
        'Switch'             'on'     'off'
        'Answer'             'yes'     'no'
    """

    def __mul__( self, other ):
	return self.get_value('python') * other

    def __rmul__( self, other ):
	return other * self.get_value('python')

    def __rnot__( self, other ):
	return not self.get_value('python')

    def __and__( self, other ):
	return self.get_value('python') & other

    def __rand__( self, other ):
	return other & self.get_value('python')

    def __or__( self, other ):
	return self.get_value('python') | other

    def __ror__( self, other ):
	return other | self.get_value('python')

    def __xor__( self, other ):
	return self.get_value('python') ^ other

    def __rxor__( self, other ):
	return other ^ self.get_value('python')

    def set_value( self, truth_value ):
	
	def test_for_valid_truth_value( self, truth_value ):

	    def set_truth_value( self, python_truth_value ):

		if python_truth_value == 1:
		    self.truth_value = 1
		elif python_truth_value == 0:
		    self.truth_value = 0
		else:
		    self.truth_value = 'Undefined'

		return self


	    def is_a_python_truth_value( truth_value ):

		if truth_value == 1: 
		    return_truth_value = 1
		elif truth_value == 0:
		    return_truth_value = 0
		else:
		    return_truth_value = 'Undefined'

		return return_truth_value


	    def is_a_string_truth_value( truth_value ):

		import string

		if string.lower( truth_value ) == 'true':
		    return_truth_value = 1
		elif string.lower( truth_value ) == 'false':
		    return_truth_value = 0
		else:
		    return_truth_value = 'Undefined'

		return return_truth_value


	    def is_a_switch_truth_value( truth_value ):

		import string

		if string.lower( truth_value ) == 'on':
		    return_truth_value = 1
		elif string.lower( truth_value ) == 'off':
		    return_truth_value = 0
		else:
		    return_truth_value = 'Undefined'

		return return_truth_value


	    def is_an_answer_truth_value( truth_value ):

		import string

		if string.lower( truth_value ) == 'yes':
		    return_truth_value = 1
		elif string.lower( truth_value ) == 'no':
		    return_truth_value = 0
		else:
		    return_truth_value = 'Undefined'

		return return_truth_value


	    def is_a_fortran_truth_value( truth_value ):

		import re
		import string

		fortran_re_tf = [ '[.][tT][rR][uU][eE][.]', '[tT]',
				  '[.][fF][aA][lL][sS][eE][.]', '[fF]' ]

		return_truth_value = 'Undefined'
		for fortran_re in fortran_re_tf:
		    match_out = re.compile(fortran_re).match(truth_value)
		    try:
			match_out.group()
		    except AttributeError:
			continue
		    else:
			if ( string.upper(truth_value) == 'T' or 
			     string.upper(truth_value) == '.TRUE.' ):
			    return_truth_value = 1
			    break
			elif ( string.upper(truth_value) == 'F' or 
			       string.upper(truth_value) == '.FALSE.' ):
			    return_truth_value = 0
			    break
		    print return_truth_value

		return return_truth_value


	    truth_value_tests = [ is_a_python_truth_value,
		is_a_string_truth_value, is_a_fortran_truth_value,
	        is_an_answer_truth_value, is_a_switch_truth_value ]

	    for truth_value_test in truth_value_tests:
		value_test_result = apply( truth_value_test, ( truth_value, ) )
		if value_test_result != 'Undefined':
		    self = set_truth_value( self, value_test_result )
		    break
	    else:
		self = set_truth_value( self, 'Undefined' )
		print 'ValueError: Illegal truth value specified.'

	    return self

	return test_for_valid_truth_value( self, truth_value ) 


    def get_value( self, truth_value_type ):

	import string

	truth_value_types = [ 
	    ( 'python',         1,        0         ),
	    ( 'string',         'true',   'false'   ),
	    ( 'dotted_fortran', '.TRUE.', '.FALSE.' ), 
	    ( 'fortran',        'T',      'F'       ),
	    ( 'switch',         'on',     'off'     ),
	    ( 'answer',         'yes',    'no'      ) ]

	truth_value_type_lc = string.lower(truth_value_type)
	for value_type in truth_value_types:
	    if truth_value_type_lc == value_type[0]:
		if self.truth_value == 1:
		    return value_type[1]
		else:
		    return value_type[2]
	else:
	    print 'ValueError: invalid truth value type'
	    return 'Undefined'
		
