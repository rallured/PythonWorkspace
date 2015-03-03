#Serial program to link up to Daly's MCA and read out the bin number.
#An option exists to have the output in terms of an integer number of 
#millivolts rather than bin number.

######IMPORTANT: Do not input a pulse over 5V into the MCA!

#Authors: Ryan Daly, Sean Mattingly
#Date: Jan 2011

import serial, time, struct

__port = 'COM6'
__baudrate = 9600
__def_timeout = 1
__noise_factor = 100

#Housekeeping routines.
def CRB_open():
    return serial.Serial(__port, __baudrate, timeout = __def_timeout) 

def CRB_close(serial_obj):
    serial_obj.close()
    return

def CRB_start_collect(serial_obj):
    serial_obj.write('b')
    return

def CRB_end_collect(serial_obj):
    serial_obj.write('e')
    return

def CRB_noise_run(serial_obj, val):
    """
    
    Sends the commands for the CRB to perform a noise measurement. The way
    that Daly's MCU is currently programmed, it will obtain 100 * val measurements
    in order to perform a noise estimate.
   
    If the input value of val will return a lenght longer than 4096, it is
    broken up into several individual queries to the MCU, so that the pySerial
    buffer will not overflow. Still, be careful not to enter too large of a 
    readout value, as the resulting array may be too large for 32 bit systems 
    to handle.
 
    Inputs:
        serial_obj  -   Already open file object to the pySerial buffer.
        val         -   Integer value that gives the desired number of values.

    Oututs:
    noise        -    Array of bin values corresponding to the noise reading. 
    
    """
    serial_obj.flush()
    val = int(val)
    v_list = [None] * (val / 10)
    
    #Populate the max length elements
    for x in range(len(v_list) - 1):
        serial_obj.write('n')
        serial_obj.write(10)
        #Per Daly's instructions, wait for 50 ms for the mCU to read.
        time.sleep(.05)
        v_list[x] = serial_obj.read(1000)
    #Populate the remainder, if any exists.
    if val % 10 != 0:
        serial_obj.write('n')
        serial_obj.write(val % 10)
        time.sleep(.05)
        v_list.append(serial_obj.read(val % 10 * 100))
    output = []
    for x in v_list:
        for y in x:
            output.append(struct.unpack('1B', y)[0])
    return output    

def CRB_change_summed(serial_obj, val):
    """Change the discriminator on the summed values, ie, the discriminator for
    the output that is the sum over all the cathode pads.

    Inputs:
        serial_obj  - Already open file object on the pySerial buffer.
        val         - Integer value that gives the desired bin to cutoff at.

    In millivolts, val comes out to mv = 5000 / 256 * val
    """
    serial_obj.write('s')
    serial_obj.write(struct.pack('1B', val))
    return

def CRB_change_individual(serial_obj, val):
    """Change the discriminator on the individual discriminators, ie,
    each discriminator's value for each cathode pad.
    
    Inputs:
        serial_obj  - Already open file object on the pySerial buffer.
        val         - integer value that gives the desired bin to cutoff at.
                        Between 0 and 256.

    In millivolts, val comes out to mv = 5000 / 256 * val
    """
    serial_obj.write('i')
    serial_obj.write(struct.pack('1B', val))
    return

def CRB_read_data_packet(serial_obj):
    """Reads out a single 3 - tuple of data from the stream.
    Outputs:
    (pos, tot, ref_sum)
                    pos - position binary on / off
                    tot - summed height - Sum of A225s
                    ref_sum - reference summed height - A250

    """
    pos = serial_obj.read(1)
    tot = serial_obj.read(1)
    ref_sum = serial_obj.read(1)
    if len(pos) != 0:
        pos = struct.unpack('1B', pos)[0]
    if len(tot) != 0:
        tot = struct.unpack('1B', tot)[0]
    if len(ref_sum) != 0:
        ref_sum = struct.unpack('1B', ref_sum)[0]
    return (pos, tot, ref_sum)

def CRB_read(serial_obj):
    """Reads out all the data from the serial buffer. Returns a tuple of three arrays
    of the different information.

    Outputs:

    (pos, tot, ref_sum)
                    pos - Array of position binary values
                    tot - array of summed heights
                    ref_sum - array of summed reference heights"""
    p, t, r = [], [], []
    while serial_obj.inWaiting() > 0:
        pos = serial_obj.read(1)
        tot = serial_obj.read(1)
        ref_sum = serial_obj.read(1)
        if len(pos) != 0:
            p.append(struct.unpack('1B', pos)[0])
        if len(tot) != 0:
            t.append(struct.unpack('1B', tot)[0])
        if len(ref_sum) != 0:
            r.append(struct.unpack('1B', ref_sum)[0])
    return (p, t, r)
        

        
def auto_read(actime):
    """Activates the main data acquisition for a certain amount of time and reads out all the data.
    
    Inputs:
        actime -    Number of seconds to acquire data for.
    """
    
    c = CRB_open()
    CRB_start_collect(c)
    time.sleep(actime)
    CRB_end_collect(c)
    p, t, r = CRB_read(c)
    CRB_close(c)
    return (p, t, r)
        
        
        
#Old method used to read off the 2 channel MCA

def acquire(acqtime, output_mvolts=False):
    """Makes a connection to the MCA using the pyserial module, reads in data for the prescribed amount of time, and then outputs the data in either millivolt or bin number.
    
    Required Inputs:
        acqtime:
            Integer. Time in seconds to acquire data.
        
    Optional Inputs:
        output_mvolts:
            Boolean. Whether the ouput will be in millivolts                          
            (True) or in bin number (False). Default is bin 
            number (False).
    """
    firstSignal = []
    secondSignal = []
    firstV = []
    secondV = []
    # open serial communications with MCAPIC
    #Previous value was 'COM1'
    ser = serial.Serial('COM3', 9600, timeout=1)
    # set one second time out
    ser.timeout = 1
    #Flush to make sure we are on the first output after the interrupt, rather     
    #than the second pulse.
    ser.flushInput()
    ser.write('b') #'b' is bytecode to begin data acquisition
    time.sleep(acqtime)
    ser.write('e') #'e' is bytecode to end data acquisition
    #ensure that we are between events
    #if ser.read(1) == 0x01 | 0x02:
    #    ser.read(1)
    while ser.inWaiting() > 0:
        #Data acquisition loop.
        first = ser.read(1)
        second = ser.read(1)
        if len(first) == 1:
            firstSignal.append(struct.unpack('1B',first)[0])
        if len(second) == 1:
            secondSignal.append(struct.unpack('1B',second)[0])
    if output_mvolts == True:
        for entry in zip(firstSignal, secondSignal):
            #Convert to mv, mult by 1000, typecast to int and append to list
            firstV.append(int(entry[0] * .01953125 * 1000))
            secondV.append(int(entry[1] * .01953125 * 1000))
    elif output_mvolts == False:
        for entry in zip(firstSignal, secondSignal):
            firstV.append(entry[0])
            secondV.append(entry[1])
        
    #Terminate connection
    ser.close()
    return (firstV, secondV)
