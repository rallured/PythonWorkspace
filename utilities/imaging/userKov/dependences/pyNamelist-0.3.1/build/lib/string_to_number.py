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
import re,string
#------------------------------------------------------------------------------
# Routine: is_a_boolean
#
#   Determines if the passed string can be interpreted as a boolean number.
#   This is one of a number of routines used to parse FORTRAN data files, such 
#   as namelists and g-files, and translate them into PYTHON format.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being a boolean number
#
#------------------------------------------------------------------------------
# Outputs:
#
#   The routine returns a dictionary, the elements of which are:
#
#        'Truth_value':  either 'true' or 'false'
#        'Value':        input string if not Boolean, otherwise 0 or 1
#        'Number_type':  'Boolean'
#        'Replications': 1
#
#   Note: some of this information is redundant if the routine itself is 
#   called (clearly, if the truth_value is 'true', a Boolean value has been 
#   found). This, however, is useful in identifying a string of unknown type 
#   with the parse number_type routine below. Also, the key, 'Replications' 
#   is added to be consistent with the 'replicator' type (see below).
#   
#------------------------------------------------------------------------------
# History (yy-mm-dd):
# 
#   00-10-02  Wrote routine (mam)
#   00-10-03  Changed return value from tuple to dictionary for easier 
#             identification of the return values.
#   01-10-17  Modified representation of Booleans from 0 or 1 to (0, 'Boolean')
#             and (1, 'Boolean')
#   02-07-10  Modified to handle unquoted strings found in modern versions of
#             namelists
#
#------------------------------------------------------------------------------
# Notes: 
#
#   Boolean values are represented in Python by the integers 0 and 1 (false
#   and true respectively) and are not an implemented class of number as are
#   integers, long integers, and floating point numbers. This routine is 
#   useful for identifying and converting FORTRAN Boolean varibles into their 
#   Python equivalent.
#
#   However, to distinguish between the intergers, 0 and 1, and a 0 or 1 that 
#   represents a Boolean value, the output is saved as a two-tuple containing
#   the (Python) truth value, 0 or 1, and a tag 'Boolean', eg: (1, 'Boolean').
#   Tagging the value in some way is necessary in order to convert Python 
#   Boolean values back to FORTRAN truth values.
#
#------------------------------------------------------------------------------

def is_a_boolean( number_string ):

    boolean_re_tf = [ '[.][tT][rR][uU][eE][.]', '[tT]',
		      '[.][fF][aA][lL][sS][eE][.]', '[fF]' ]

    for boolean_re in boolean_re_tf:
	match_out = re.compile(boolean_re).match(number_string)
	try:
	    match_out.group()
	except AttributeError:
	    is_a_boolean = 'false'
	    x = number_string
	else:
	    true_or_false = string.upper(number_string)
	    if match_out.group() == number_string:
		is_a_boolean = 'true'
		if true_or_false == '.TRUE.' or true_or_false == 'T':
		    x = ( 1, 'Boolean' )
		    break
		else:
		    x = ( 0, 'Boolean' )
		    break

    return { 'Truth_value': is_a_boolean, 'Value': x, 'Number_type': 
	     'Boolean', 'Replications': 1 }

#------------------------------------------------------------------------------
# Routine: is_an_octal
#
#   Determines if the passed string can be interpreted as an octal number.
#   This is one of a number of routines used to parse FORTRAN data files, such 
#   as namelists and g-files, and translate them into PYTHON format.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being an octal number
#
#------------------------------------------------------------------------------
# Outputs:
#
#   The routine returns a dictionary, the elements of which are:
#
#        'Truth_value':  either 'true' or 'false'
#        'Value':        input string if not Octal, otherwise a Python octal
#                        number
#        'Number_type':  'Octal'
#        'Replications': 1
#
#   Note: some of this information is redundant if the routine itself is 
#   called (clearly, if the truth_value is 'true', an Octal value has been 
#   found. This, however, is useful in identifying a string of unknown type 
#   with the parse number_type routine below. Also, the key, 'Replications' 
#   is added to be consistent with the 'replicator' type (see below).
#   
#------------------------------------------------------------------------------
# History (yy-mm-dd):
# 
#   00-09-27  Started to write routine (mam)
#   00-09-28  Added check for a leading minus sign. Then rewrote everything
#             using the regular expression module, re. (mam)
#   00-10-03  Changed return value from tuple to dictionary for easier 
#             identification of the return values.
#
#------------------------------------------------------------------------------
# Notes: 
#
#   Octals are represented by STRINGs consisting of 1 or more digits in the 
#   range 0-7 preceded by 0 that may be further preceded by a minus sign. In
#   Python there are also LONG OCTALs which end in the letter 'L'
#
#------------------------------------------------------------------------------

def is_an_octal( number_string ):

    octal_re = '-?0[0-7]+[lL]?'

    match_out = re.compile(octal_re).match(number_string)

    try:
	match_out.group()
    except AttributeError:
	is_an_octal = 'false'
    else:
	if match_out.group() == number_string:
	    is_an_octal = 'true'
	else:
	    is_an_octal = 'false'

    return { 'Truth_value': is_an_octal, 'Value': number_string, 
	     'Number_type': 'Octal', 'Replications': 1 }

#------------------------------------------------------------------------------
# Routine: is_a_hexadecimal
#
#   Determines if the passed string can be interpreted as a hexidecimal number.
#   This is one of a number of routines used to parse FORTRAN data files, such 
#   as namelists and g-files, and translate them into PYTHON format.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being a hexadecimal number
#
#------------------------------------------------------------------------------
# Outputs:
#
#   The routine returns a dictionary, the elements of which are:
#
#        'Truth_value':  either 'true' or 'false'
#        'Value':        input string if not a hexadecimal, otherwise a PYTHON 
#                        hexadecimal value
#        'Number_type':  'Hex'
#        'Replications': 1
#
#   Note: some of this information is redundant if the routine itself is 
#   called (clearly, if the truth_value is 'true', a hexadecimal value has 
#   been found. This, however, is useful in identifying a string of unknown 
#   type with the parse number_type routine below. Also, the key, 
#   'Replications' is added to be consistent with the 'replicator' type (see 
#   below).
#
#------------------------------------------------------------------------------
# History (yy-mm-dd):
# 
#   00-09-27  Started to write routine (mam)
#   00-09-28  Added check for a leading minus sign. Then rewrote everything
#             using the regular expression module, re. (mam)
#   00-10-03  Changed return value from tuple to dictionary for easier 
#             identification of the return values.
#
#------------------------------------------------------------------------------
# Notes: 
#
#   Hexadecimal numbers are represented by STRINGs consisting of 1 or more 
#   digits in the range 0-7 preceded by 0 that may be further preceded by a 
#   minus sign. In Python there are also LONG OCTALs which end in the letter 
#   'L'.
#
#------------------------------------------------------------------------------

def is_a_hexadecimal( number_string ):

    hexadecimal_re = '-?0[xX][0123456789abcdefABCDEF]+[lL]?'

    match_out = re.compile(hexadecimal_re).match(number_string)

    try:
	match_out.group()
    except AttributeError:
	is_a_hexadecimal = 'false'
    else:
	if match_out.group() == number_string:
	    is_a_hexadecimal = 'true'
	else:
	    is_a_hexadecimal = 'false'

    return { 'Truth_value': is_a_hexadecimal, 'Value': number_string, 
	     'Number_type': 'Hex', 'Replications': 1 }

#------------------------------------------------------------------------------
# Routine: is_an_integer
#
#   Determines if the passed string can be interpreted as an integer. 
#   This is one of a number of routines used to parse FORTRAN data files, such 
#   as namelists and g-files, and translate them into PYTHON format.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being an integer
#
#------------------------------------------------------------------------------
# Outputs:
#
#   The routine returns a dictionary, the elements of which are:
#
#        'Truth_value':  either 'true' or 'false'
#        'Value':        input string if not an integer, otherwise a PYTHON
#                        integer
#        'Number_type':  'Integer'
#        'Replications': 1
#
#   Note: some of this information is redundant if the routine itself is 
#   called (clearly, if the truth_value is 'true', an integer value has been 
#   found. This, however, is useful in identifying a string of unknown type 
#   with the parse number_type routine below. Also, the key, 'Replications' 
#   is added to be consistent with the 'replicator' type (see below).
#
#------------------------------------------------------------------------------
# History (yy-mm-dd):
#
#   00-09-28  Wrote routine (mam)
#   00-10-03  Changed return value from tuple to dictionary for easier 
#             identification of the return values.
#
#------------------------------------------------------------------------------

def is_an_integer( number_string ):

    integer_re = '-?\d+'

    match_out = re.compile(integer_re).match(number_string)

    try:
	match_out.group()
    except AttributeError:
	is_an_integer = 'false'
	x = number_string
    else:
	if match_out.group() == number_string:
	    is_an_integer = 'true'
	    x = int(number_string)
	else:
	    is_an_integer = 'false'
	    x = number_string

    return { 'Truth_value': is_an_integer, 'Value': x, 'Number_type': 
	     'Integer', 'Replications': 1 }

#------------------------------------------------------------------------------
# Routine: is_a_long
#
#   Determines if the passed string can be interpreted as a long (python) 
#   integer. This is one of a number of routines used to parse FORTRAN data 
#   files, such as namelists and g-files, and translate them into PYTHON 
#   format.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being a long integer
#
#------------------------------------------------------------------------------
# Outputs:
#
#   The routine returns a dictionary, the elements of which are:
#
#        'Truth_value':  either 'true' or 'false'
#        'Value':        input string if not a long interger, otherwise a
#                        PYTHON long integer value
#        'Number_type':  'Long'
#        'Replications': 1
#
#   Note: some of this information is redundant if the routine itself is 
#   called (clearly, if the truth_value is 'true', a long integer value has 
#   been found. This, however, is useful in identifying a string of unknown 
#   type with the parse number_type routine below. Also, the key, 
#   'Replications' is added to be consistent with the 'replicator' type (see 
#   below).
#
#------------------------------------------------------------------------------
# History (yy-mm-dd):
#
#   00-09-28  Wrote routine (mam)
#   00-10-03  Changed return value from tuple to dictionary for easier 
#             identification of the return values.
#
#------------------------------------------------------------------------------
# Notes: 
#
#   In python, long integers are terminated with the letter, 'L'.
#
#------------------------------------------------------------------------------

def is_a_long( number_string ):

    long_re = '-?\d+[lL]?'

    match_out = re.compile(long_re).match(number_string)

    try:
	match_out.group()
    except AttributeError:
	is_a_long = 'false'
	x = number_string
    else:
	if match_out.group() == number_string:
	    is_a_long = 'true'
	    x = long(number_string)
	else:
	    is_a_long = 'false'
	    x = number_string

    return { 'Truth_value': is_a_long, 'Value': x, 'Number_type': 'Long',
	     'Replications': 1 }

#------------------------------------------------------------------------------
# Routine: is_a_real
#
#   Determines if the passed string can be interpreted as a real (i.e. a 
#   floating point number). This is one of a number of routines used to parse 
#   FORTRAN data files, such as namelists and g-files, and translate them into
#   PYTHON format.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being a floating point number
#
#------------------------------------------------------------------------------
# Outputs:
#
#   The routine returns a dictionary, the elements of which are:
#
#        'Truth_value':  either 'true' or 'false'
#        'Value':        input string if not a floating point number, 
#                        otherwise the Python value of the floating point 
#                         number
#        'Number_type':  'Float'
#        'Replications': 1
#
#   Note: some of this information is redundant if the routine itself is 
#   called (clearly, if the truth_value is 'true', a floating point value has 
#   been found. This, however, is useful in identifying a string of unknown 
#   type with the parse number_type routine below. Also, the key, 
#   'Replications' is added to be consistent with the 'replicator' type (see 
#   below).
#
#------------------------------------------------------------------------------
# History (yy-mm-dd):
#
#   00-09-28  Wrote routine (mam)
#   00-10-03  Changed return value from tuple to dictionary for easier 
#             identification of the return values.
#   00-10-03  Changed return value from tuple to dictionary for easier 
#             identification of the return values.
#
#------------------------------------------------------------------------------

def is_a_real( number_string ):

    real_re = '-?\d*[.]?\d*([eEdDgG][+-]?\d{1,3})?'

    match_out = re.compile(real_re).match(number_string)

    try:
	match_out.group()
    except AttributeError:
	is_a_real = 'false'
	x = number_string
    else:
	if match_out.group() == number_string:
	    is_a_real = 'true'
	    x = float(number_string)
	else:
	    is_a_real = 'false'
	    x = number_string

    return { 'Truth_value': is_a_real, 'Value': x, 'Number_type': 'Float',
	     'Replications': 1 }

#------------------------------------------------------------------------------
# Routine: is_a_string1
#
#   Determines if the passed string can be interpreted as a string literal. This
#   is one of a number of routines used to parse FORTRAN data files, such as 
#   namelists and g-files, and translate them into PYTHON format.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being a string literal
#
#------------------------------------------------------------------------------
# Outputs:
#
#   The routine returns a dictionary, the elements of which are:
#
#        'Truth_value':  either 'true' or 'false'
#        'Value':        input string if not a floating point number, 
#                        otherwise the Python value of the floating point 
#                         number
#        'Number_type':  'String'
#        'Replications': 1
#
#   Note: some of this information is redundant if the routine itself is 
#   called (clearly, if the truth_value is 'true', a floating point value has 
#   been found. This, however, is useful in identifying a string of unknown 
#   type with the parse number_type routine below. Also, the key, 
#   'Replications' is added to be consistent with the 'replicator' type (see 
#   below).
#
#------------------------------------------------------------------------------
# History (yy-mm-dd):
#
#   01-10-15  Added this routine to detect string literals which are among the
#             possible types of namelist inputs.
#
#------------------------------------------------------------------------------

def is_a_string1( number_string ):

    if number_string[0][0] == "'":
        x = string.join( number_string, '' )
	is_a_string = 'true'
    else:
	is_a_string = 'false'
	x = number_string

    return { 'Truth_value': is_a_string, 'Value': x, 'Number_type': 'String1',
	     'Replications': 1 }

#------------------------------------------------------------------------------
# Routine: is_a_string2
#
#   Older type FORTRAN namelists used quoted strings, making them easy to 
#   detect. This type of string is treated in is_a_string1. Modern namelist no
#   longer quote strings, hence another routine to treat these. This routine
#   determines if the passed string can be interpreted as an unquoted string 
#   literal. This is one of a number of routines used to parse FORTRAN data 
#   files, such as namelists and g-files, and translate them into PYTHON format.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being a string literal
#
#------------------------------------------------------------------------------
# Outputs:
#
#   The routine returns a dictionary, the elements of which are:
#
#        'Truth_value':  either 'true' or 'false'
#        'Value':        input string if not a floating point number, 
#                        otherwise the Python value of the floating point 
#                         number
#        'Number_type':  'String'
#        'Replications': 1
#
#   Note: some of this information is redundant if the routine itself is 
#   called (clearly, if the truth_value is 'true', a floating point value has 
#   been found. This, however, is useful in identifying a string of unknown 
#   type with the parse number_type routine below. Also, the key, 
#   'Replications' is added to be consistent with the 'replicator' type (see 
#   below).
#
#------------------------------------------------------------------------------
# History (yy-mm-dd):
#
#   02-07-10  Added this routine to detect quoted string literals which are 
#             among the possible types of namelist inputs.
#
#------------------------------------------------------------------------------

def is_a_string2( number_string ):

    valid_first_characters = \
      "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/"

    if number_string[0] in valid_first_characters:
        x = string.join( "'" + number_string + "'", '' )
	is_a_string = 'true'
    else:
	is_a_string = 'false'
	x = number_string

    return { 'Truth_value': is_a_string, 'Value': x, 'Number_type': 'String2',
	     'Replications': 1 }

#------------------------------------------------------------------------------
# Routine: is_a_replicator
#
#   Determines if the passed string can be interpreted as a "replicator", a
#   shorthand notation used in namelist to represent a series of identical 
#   values. This is one of a number of routines used to parse FORTRAN data 
#   files, such as namelists and g-files, and translate them into PYTHON 
#   format.
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being a replicator
#
#------------------------------------------------------------------------------
# Outputs:
#
#   The routine returns a dictionary, the elements of which are:
#
#        'Truth_value':  either 'true' or 'false'
#        'Value':        value of valid number type
#        'Number_type':  'Relicator'
#        'Replications': Integer value corresponding to the number of 
#                        replications of 'Value'
#
#------------------------------------------------------------------------------
# History (yy-mm-dd):
#
#   00-10-02  Wrote routine (mam)
#   00-10-03  Changed return value from tuple to dictionary for easier 
#             identification of the return values.
#
#------------------------------------------------------------------------------
# Notes:
#
#   FORTRAN replicators are of the following form: 
#
#                            (integer)*(number)
#
#   where (number) is one of the FORTRAN number types (Boolean, integer, float,
#   ...), where '*' can be thought of a as multiplier, and (integer) is an
#   the number of times (number) is replicated. Parsng a replicator then 
#   reduces to using '*' to split the string into two parts. The first is
#   checked to see if it corresponds to a standard number type, and the second
#   is checked to see if it is an integer'
#
#------------------------------------------------------------------------------

def is_a_replicator( number_string ):

    asteric_check = '.+[*].+'

    asteric_check_out = re.compile(asteric_check).match(number_string)

    try:
	asteric_check_out.group()
    except AttributeError:
	is_a_replicator = 'false'
	x = number_string
	n = -1
    else:
	split_string = string.split( number_string, '*' )
	integer_check = is_an_integer( split_string[0] )
	parsed_number = parse_number_string( split_string[1] )
	if ( integer_check['Truth_value'] == 'true' and 
            parsed_number['Truth_value'] != None ):
	    is_a_replicator = 'true'
	    x = parsed_number['Value']
            n = integer_check['Value']
	else:
	    is_a_replicator = 'false'
	    x = number_string
	    n = -1

    return { 'Truth_value': is_a_replicator, 'Value': x, 'Number_type': 
	     'Replicator', 'Replications': n }

#------------------------------------------------------------------------------
# Routine: parse_number_string
#
#   Determines if the passed string can be interpreted as one of a variety of
#   standard number types (octal, hexadecimal, integer, long, floating point).
#
#------------------------------------------------------------------------------
# Inputs:
#
#   number_string = string suspected of being some type of standard formatted
#     number
#
#------------------------------------------------------------------------------
# Outputs:
#
#   A two-tuple of the form ( number, number_type ) where number is the python
#   representation of the input, number_string, and number_type is type of
#   number (octal, hexadecimal, integer, long, floating point). In the event
#   that number_string does not correspond to any of the recognized forms, the
#   tuple ( number_string, None ) is returned.
#
#------------------------------------------------------------------------------
# History (yy-mm-dd):
#
#   00-09-29  Scrapped former routine and rewrote using the apply function 
#             (mam)
#   01-10-15  Included ability to detect string literals
#
#------------------------------------------------------------------------------
# Note:
#
#   is_a_string2 must follow is_a boolean as a character such as T or F will be
#   identified a character.
#
#------------------------------------------------------------------------------

def parse_number_string( number_string ):

    function_names = [ is_an_integer,  is_a_real, is_a_long, is_an_octal, 
      is_a_hexadecimal, is_a_boolean, is_a_replicator, is_a_string1,
      is_a_string2 ]

    for func in function_names:
	parsed_number = apply( func, ( number_string, ) )
	if parsed_number['Truth_value'] == 'true':
	    x = parsed_number
	    break
    else:
	print 'ERROR in string_to_number: the string, <',  number_string, 
	print '>, does not correspond to a standard type.'
	print 'Returning original string'
	print ''
	x = { 'Truth_value': None, 'Value': number_string, 'Number_type': 
	      None, 'Replications': 1 }
	    
    return x
