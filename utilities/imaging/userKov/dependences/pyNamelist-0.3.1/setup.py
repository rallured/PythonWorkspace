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

from distutils.core import setup, Extension
from distutils.command.install_headers import install_headers

setup (name = "pyNamelist",
       version = "0.3.1",
       description = """Allows python to read and write Fortran namelist files""",
       author = "Mike Makowski",
       author_email = "makowski@fusion.gat.com",
       url = "",
       license = "GPL",
       packages = [''],
       package_dir = {'': 'Lib'},
       extra_path = 'pyNamelist',
       )
