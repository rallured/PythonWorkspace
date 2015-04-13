import io
import tarfile
import base64
import pyjs9

# this is the default smartx server
srv = pyjs9.JS9('smartx.cfa.harvard.edu')

def solve_pzt(userid, cmd):
    # send command to server
    r = srv.send({'cmd': 'solve_pzt', 'args': [userid, cmd]}, msg='smartx')
    if 'stderr' in r:
        # if stderr was returned, it's an error
        raise ValueError(r['stderr'])
    else:
        # stdout is an encoded string containing a tar file
        encstr = r['stdout'].strip()
        # decode it into bytes
        bytes = base64.decodestring(encstr)
        # get a file object from the bytes
        fobj = io.BytesIO(bytes)
        # open the file object as a tar file
        tar = tarfile.open(fileobj=fobj)
        # extract all files onto disk
        tar.extractall()
        # close tar file
        tar.close()
        # everything is OK
        return 'OK'

def file_list(userid, tmpl):
    # send command to server
    r = srv.send({'cmd': 'file_list', 'args': [userid, tmpl]}, msg='smartx')
    if 'stderr' in r:
        # if stderr was returned, it's an error
        raise ValueError(r['stderr'])
    else:
        # stdout is the text we want
        return r['stdout'].strip()

if __name__ == '__main__':
    import datetime
    import time
    # srv = pyjs9.JS9()
    cmd = 'max_strain=10 save_premath=Math.dat @fit-5x10mm-gap02-flex.par'
    print "solve_pzt:\n    host: %s\n    cmd: %s" % (srv.host, cmd)
    t = time.strftime("%H:%M:%S:%MS", time.localtime())
    print "start time: %s" % t
    print solve_pzt("eric", cmd)
    t = time.strftime("%H:%M:%S:%MS", time.localtime())
    print "end time: %s" % t

    print "premath files:"
    print file_list("eric", "*.dat")

    cmd = 'max_strain=10 premath=Math.dat @fit-5x10mm-gap02-flex.par -o ./Y'
    t = time.strftime("%H:%M:%S:%MS", time.localtime())
    print "start time: %s" % t
    print "solve_pzt:\n    host: %s\n    cmd: %s" % (srv.host, cmd)
    print solve_pzt("eric", cmd)
    t = time.strftime("%H:%M:%S:%MS", time.localtime())
    print "end time: %s" % t
