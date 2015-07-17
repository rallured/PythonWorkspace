from ctypes import *

if __name__=="__main__":
    dll = windll.LoadLibrary("helloWorld.dll")
    a=c_float(3.1)
    b=5
    dll.test(a,c_int(b))
    print a, dll.test.restype," ",b, dll.test.argtypes
    

    '''
    fortran:
    subroutine test(y,x) 
    !MS$ ATTRIBUTES DLLEXPORT, ALIAS:'test', STDCALL::test 
    integer,intent(in):: x 
    real, intent(out):: y 
    y = x * 1.35 
    return 
    end subroutine test
    -------------------------------------------------------
    a=3.1
    result = dll.test(byref(c_float(a)),byref(c_int(2)))
    print result," ",a, dll.test.restype," ", dll.test.argtypes
    
    18351048   3.1 <class 'ctypes.c_long'>   None
    il primo numero cambia ogni volta se si rilancia piu' volte
    -----------------------------------------------------------
    a=c_float(3.1)
    result = dll.test(byref(a),byref(c_int(2)))
    print result," ",a, dll.test.restype," ", dll.test.argtypes
    
    18350968   c_float(3.0999999046325684) <class 'ctypes.c_long'>   None
    -------------------------------------------------------
    a=c_float(3.1)
    b=c_int(2)
    result = dll.test(byref(a),byref(b))
    print result," ",a," ",b," ", dll.test.restype," ", dll.test.argtypes
    18353848   c_float(3.0999999046325684)   c_long(2)   <class 'ctypes.c_long'>   None
    
    '''