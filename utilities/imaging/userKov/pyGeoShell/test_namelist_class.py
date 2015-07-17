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
##################################################################################################
# CLASSES WRAPPER FOR NAMELIST
import namelist,copy,types
class Namelist:
    """
-------------------------------------------------------------------------------------
 class Namelist: Read, write and examine the values in a FORTRAN namelist file.
                 See __init__.__doc__ for details of sturcture instances.
                 
     self.namelist contains the namelist dictionary as in the namelist module
                   
 Methods (see individual .__doc__ 's for details):
     
     __init__(self,file=None) : file = namelist file to open; = None if file=None
     
     __getitem__(self,key) : get value of namelist variable key, assumed unique
     
     __setitem__(self,key,value):set namelist varialble key=value, assumed unique and exists
     
     __delitem__(self,key):delete namelist variable key, assumed unique
     
     __len__(self) : returns dictionary of number of variables in namelists
     
     __repr__(self) : prints variables and values to terminal
     
     __add__(self,other): combine, 'add', two instances. self and other cannot contain same namelists
     
     __sub__(self,other): 'subtract' two instances. Namelists in other are removed from self.

     write(self,file) : file = namelist file to write
     
     keys(self,key=None) : list namelist vairables for all namelists in the file
     
     order(self,nl): print order variable key would appear in file containing namelist nl
-------------------------------------------------------------------------------------
    """
    def __init__(
        # Create a namelist object by reading in a FORTRAN namellist file
        # self.namelist is a dictionary:
        # self.namelist[namelist_name] is a dictionary with keys = the variable names
        # in namelist namelist_name and values = values of these variables; and a few
        # special variables:self.namelist[namelist_name]['_nl_header']=text which preceeded
        # this namelist in the file, self.namelist[namelist_name]['_nl_sequence']=what order
        # this namelist_name namelist appeared in file, self.namelist[namelist_name]['_nl_order']
        # list of variables in self.namelist[namelist_name] in the order they were read in--this
        # is used to write the namelist variables out in the same order they were read in.
        # also self.namelist['_delimiter'] = namelist delimiter (& or $) and
        # self.namelist['tailer'] = text following all the namelists in file
        self,
        file = None # Namelist file name, if None, self.nameilst=None
        ):
        """
        __init__(
        # Create a namelist object by reading in a FORTRAN namellist file
        # self.namelist is a dictionary:
        # self.namelist[namelist_name] is a dictionary with keys = the variable names
        # in namelist namelist_name and values = values of these variables; and a few
        # special variables:self.namelist[namelist_name]['_nl_header']=text which preceeded
        # this namelist in the file, self.namelist[namelist_name]['_nl_sequence']=what order
        # this namelist_name namelist appeared in file, self.namelist[namelist_name]['_nl_order']
        # list of variables in self.namelist[namelist_name] in the order they were read in--this
        # is used to write the namelist variables out in the same order they were read in.
        # also self.namelist['_delimiter'] = namelist delimiter (& or $) and
        # self.namelist['tailer'] = text following all the namelists in file
        self,
        file = None # Namelist file name, if None, self.nameilst=None
        ) 
        """
        if file is None:
            self.namelist = None
        else:
            self.namelist = namelist.read(file)
        
    def write(
        # write a fortran namelist file from a namelist instance
        self,
        file = 'namelist', # Filename to write
        ):
        """
        write(
        # Write a fortran namelist file from a namelist instance
        self,
        file = 'namelist', # Filename to write
        ) 
        """
        namelist.write(self.namelist,file)

    def keys(
        # Print out the keys to all namelists with variables sorted alphabetically
        # or return keys dictionary
        self,
        nl = None, # If nl is none then print out all keys, else retrun keys dictionary
                   # for namelist nl
        ):
        """
        keys(
        # Print out the keys to all namelists with variables sorted alphabetically
        # or return keys dictionary
        self
        nl = None, # If nl is none then print out all keys, else retrun keys dictionary
                   # for namelist nl
        )
        """
        if nl is not None:
            s = self.namelist[nl].keys()
            s.sort()
            ss=s[:]
            for v in ss:
                if v in [ '_nl_header', '_nl_sequence', '_nl_order']:
                    s.remove(v)
            return s
        for nl in self.namelist.keys():
            if nl in ['_delimiter','tailer']:continue
            print '-'*80
            print nl
            print '-'*80
            names = self.namelist[nl].keys()
            names.sort()
            for v in names:
                if v in [ '_nl_header', '_nl_sequence', '_nl_order']:continue
                print v

    def __repr__(
        # Print out the keys and values for all namelists with variables sorted alphabetically
        self
        ):
        """
        __repr__(
        # print out the keys and values for all namelists with variables sorted alphabetically
        self
        )
        """
        rep = ''
        for nl in self.namelist.keys():
            if nl in ['_delimiter','tailer']:continue
            rep=rep+'-'*80+'\n'
            rep=rep+ str(nl)+'\n'
            rep=rep+ '-'*80+'\n'
            names = self.namelist[nl].keys()
            names.sort()
            for v in names:
                if v in [ '_nl_header', '_nl_sequence', '_nl_order']:continue
                rep=rep+ '%-10s = %s\n' % (v,self.namelist[nl][v])
        return rep

    def __getitem__(
        # Get value of namelist vaiable key. If a key has a single value,
        # e.g. v['xyz'], v.namelist[nl]['xyz'] is returned from the first nl
        # which has the key. If key has two values, e.g. v['nl1','xyz'], the
        # value of 'xyz' from namelist 'nl1' is returned.
        self,
        key,   # Variable name or namelist,variable
        ):
        """
        __getitem__(
        # Get value of namelist vaiable key. If a key has a single value,
        # e.g. v['xyz'], v.namelist[nl]['xyz'] is returned from the first nl
        # which has the key. If key has two values, e.g. v['nl1','xyz'], the
        # value of 'xyz' from namelist 'nl1' is returned.
        self,
        key,   # Variable name or namelist,variable
        )
        """
        if type(key) is types.StringType:
            # FOR SINGLE ARGUMENT RETURN VALUE IN FIRST NAMELIST FOUND
            for nl in self.namelist.keys():
                if nl in ['_delimiter','tailer']:continue
                if key in self.namelist[nl].keys():
                    return self.namelist[nl][key]
            raise KeyError (key+' is not in any namelist')
        else:
            if ( type(key[0]) is not types.StringType or
                 type(key[1]) is not types.StringType ):
                raise KeyError('Both values of key must be string type')
            return self.namelist[key[0]][key[1]]

    def __setitem__(
        # Set value of namelist vaiable key to value. If key has only one value,
        # e.g. v['xyz']=3 then the value is set in the first namelist with 'xyz';
        # if no namelist has 'xyz' an error is raised. If key has two values, e.g.
        # v['nl1','xyz']=3, 'xyz' in namelist 'nl1' is set to 3. If 'xyz' does not
        # exist in 'nl1' it will be added and written at the end of the namelist.
        self,
        key,    # Variable name or namelist,variable
        value,  # Set varialble key to this value
        ):
        """
        __setitem__(
        # Set value of namelist vaiable key to value. If key has only one value,
        # e.g. v['xyz']=3 then the value is set in the first namelist with 'xyz';
        # if no namelist has 'xyz' an error is raised. If key has two values, e.g.
        # v['nl1','xyz']=3, 'xyz' in namelist 'nl1' is set to 3. If 'xyz' does not
        # exist in 'nl1' it will be added and written at the end of the namelist.
        self,
        key,    # Variable name or namelist,variable
        value,  # Set varialble key to this value
        )
        """
        if type(key) is types.StringType:
            for nl in self.namelist.keys():
                if nl in ['_delimiter','tailer']:continue
                if key in self.namelist[nl].keys():
                    self.namelist[nl][key]=value
                    return
            raise KeyError (key+' is not in any namelist')
        else:
            nl = key[0]
            k  = key[1]
            if ( type(nl) is not types.StringType or
                 type(k) is not types.StringType ):
                raise KeyError('Both values of key must be string type')
            if ( '_nl_order' in self.namelist[nl].keys() and
                 k not in self.namelist[nl]['_nl_order']   ):
                self.namelist[nl]['_nl_order'].append(k)
            self.namelist[nl][k] = value
            

    def __delitem__(
        # Delete namelist vaiable key. If key has only one value,
        # e.g. v['xyz'] then the variable is deleted in the first namelist with 'xyz';
        # if no namelist has 'xyz' an error is raised. If key has two values, e.g.
        # v['nl1','xyz'], 'xyz' in namelist 'nl1' is deleted, and removed from _nl_order.
        self,
        key  # Variable name to delete
        ):
        """
        __delitem__(
        # Delete namelist vaiable key. If key has only one value,
        # e.g. v['xyz'] then the variable is deleted in the first namelist with 'xyz';
        # if no namelist has 'xyz' an error is raised. If key has two values, e.g.
        # v['nl1','xyz'], 'xyz' in namelist 'nl1' is deleted, and removed from _nl_order.
        self,
        key  # Variable name to delete
        )
        """
        if type(key) is types.StringType:
            for nl in self.namelist.keys():
                if nl in ['_delimiter','tailer']:continue
                if key in self.namelist[nl].keys():
                    del(self.namelist[nl][key])
                    if '_nl_order' in self.namelist[nl].keys():
                        del self.namelist[nl]['_nl_order'][ self.namelist[nl]['_nl_order'].index(key) ]
                    return
            raise KeyError (key+' is not in any namelist')
        else:
            nl = key[0]
            k  = key[1]
            if ( type(nl) is not types.StringType or
                 type(k) is not types.StringType ):
                raise KeyError('Both values of key must be string type')
            del self.namelist[nl][k]
            if '_nl_order' in self.namelist[nl].keys():
                del self.namelist[nl]['_nl_order'][ self.namelist[nl]['_nl_order'].index(k) ]
            
    def len(
        # Return dictionary with number of varialbles in each namelist,
        # if self.namelist is None return 0.
        self
        ):
        """
        len(
        # Return dictionary with number of varialbles in each namelist,
        # if self.namelist is None return 0.
        self
        )
        """
        if self.namelist is None:return 0
        slen = {}
        for nl in self.namelist.keys():
            if nl in ['_delimiter','tailer']:continue
            slen[nl] = 0
            for v in self.namelist[nl].keys():
                if v in [ '_nl_header', '_nl_sequence', '_nl_order']:continue
                slen[nl] = slen[nl] + 1
        return slen

    def order(
        # Return the order variable key in namelist nl, will be appear if
        # written to a file with .write
        self,
        nl,         # Namelist name
        key = None, # Varialble name, if None then _nl_order list is returned
        ):
        """
        order(
        # Return the order variable key in namelist nl, will be appear if
        # written to a file with .write
        self,
        nl,         # Namelist name
        key = None, # Varialble name, if None then _nl_order list is returned
        )
        """
        if key is None:
            return self.namelist[nl]['_nl_order']
        if '_nl_order' in self.namelist[nl].keys():
            return self.namelist[nl]['_nl_order'].index(key)
        raise KeyError('_nl_order not defined in namelist'+nl)
            
                

    def __add__(
        # 'Add' togeather two namelist instances. if the two Namelist instances
        # contain different namelists, adding combines these into a single Namelist instance.
        # if the same namelist is in both instances an error is raised.
        self,
        other,
        ):
        """
        __add__(
        # 'Add' togeather two namelist instances. if the two Namelist instances
        # contain different namelists, adding combines these into a single Namelist instance.
        # if the same namelist is in both instances an error is raised.
        self,
        other,
        )
        """
        if ( type( self  ) is not types.InstanceType or
             type( other ) is not types.InstanceType ):
            raise TypeError( 'unsupported operand type(s) for +: '+
                             str(type(self))+' '+str(type(other))   )
        if ( not isinstance(self,Namelist) or not isinstance(other,Namelist) ):
            raise TypeError( ' Both instances must be Namelist class of subclass of Namelist class ' +
                             self.__class.__name__ + ' ' + other.__class__.__name__ )
        if ( self.__class__.__name__ <> other.__class__.__name__):
            raise TypeError( ' Both instances must be of the same class ' +
                             self.__class__.__name__ + ' ' + other.__class__.__name__ )

        v = copy.deepcopy( self )
        v.namelist = copy.deepcopy( self.namelist )
        v.namelist[ 'tailer' ] = v.namelist[ 'tailer' ][:] + other.namelist[ 'tailer' ][:]
        vkeys =     v.namelist.keys()
        okeys = other.namelist.keys()
        
        nl_seq_max=-1
        for k in vkeys:
            if k == 'tailer' or k == '_delimiter':continue
            nl_seq_max=max( nl_seq_max, v.namelist[k]['_nl_sequence'] )
        nl_seq_max = nl_seq_max+1
            
        for k in okeys:
            if k == 'tailer' or k == '_delimiter':continue
            if k in vkeys:
                raise RuntimeError( 'both Namelist instances contain the same namelist : ' + k )
            v.namelist[k] = copy.deepcopy( other.namelist[k] )
            v.namelist[k]['_nl_sequence'] = v.namelist[k]['_nl_sequence'] + nl_seq_max
        return v

    def __sub__(
        # 'Subtract' two Namelist instances. if c=a-b, namelists in b will
        # be removed from a to produce c.
        self,
        other
        ):
        """
        __sub__(
        # 'Subtract' two Namelist instances. if c=a-b, namelists in b will
        # be removed from a to produce c.
        self,
        other
        )
        """
        if ( type( self  ) is not types.InstanceType or
             type( other ) is not types.InstanceType ):
            raise TypeError( 'unsupported operand type(s) for +: '+
                             str(type(self))+' '+str(type(other))   )
        if ( not isinstance(self,Namelist) or not isinstance(other,Namelist) ):
            raise TypeError( ' Both instances must be Namelist class of subclass of Namelist class ' +
                             self.__class.__name__ + ' ' + other.__class__.__name__ )
        if ( self.__class__.__name__ <> other.__class__.__name__):
            raise TypeError( ' Both instances must be of the same class ' +
                             self.__class.__name__ + ' ' + other.__class__.__name__ )
        v = Namelist( None )
        v.namelist = copy.deepcopy( self.namelist )
        vkeys =     v.namelist.keys()
        okeys = other.namelist.keys()
        for k in okeys:
            if k == 'tailer' or k == '_delimiter':continue
            if k in vkeys:del(v.namelist[k])
        nl_seq={}
        for k in v.namelist.keys():
            if k == 'tailer' or k == '_delimiter':continue
            nl_seq[ v.namelist[k]['_nl_sequence'] ] = k
        seq = nl_seq.keys()
        seq.sort()
        k=0
        for s in seq:
            v.namelist[ nl_seq[s] ]['_nl_sequence']=k
            k=k+1
        return v
                
if __name__=="__main__":
    a=Namelist("imp_OffAxis.txt")
    print a.kov
    