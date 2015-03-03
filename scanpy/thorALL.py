""" code to run Thorlabs stages in conjunction with the VB allmounts project.
First, run the application ThorLabsSocket.
It won't put up a window until a TCP socket is established
Then init() to setup, r_move(pos), v_move(pos) or h_move(pos) to move, close() to close.
"""

import socket, time

def init():
    """Initializes socket connections. Make sure to run the ThorLabsAllMounts VB program before this;
    that is the program this method attempts to connect to.
    """
    global thorSocketVERT
    host = 'localhost'
    port = 8887
    bufsize = 256
    addr = (host, port)
    thorSocketVERT = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    thorSocketVERT.connect(addr)

def v_move(pos):
    #Move the vertical stage to a position given by pos.
    global thorSocketVERT
    thorSocketVERT.send("MOVEVERT "+str(pos)+"$")
    time.sleep(10) # wait for move, should really read back position

def r_move(pos):
    #Move the rotational stage to a position. Must be in the range [0, 360).
    global thorSocketVERT
    if (pos < 0 or pos >= 360):
        print "Not a valid move command! Must be 0 <= x < 360."
        print "Returning without moving."
        return
    thorSocketVERT.send("MOVEROT "+str(pos)+"$")
    time.sleep(10) # wait for move, should really read back position

def h_move(pos):
    #Move the horizontal stage to a position.
    global thorSocketVERT
    thorSocketVERT.send("MOVEHOR "+str(pos)+"$")
    time.sleep(10) # wait for move, should really read back position

  
def close():
    global thorSocketVERT
    thorSocketVERT.close()

