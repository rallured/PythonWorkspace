

subroutine test(y,x) 
!MS$ ATTRIBUTES DLLEXPORT, ALIAS:'test', STDCALL::test 
integer,intent(in):: x 
real, intent(out):: y 
y = x * 1.35 
return 
end subroutine test