module gzipmod
use, intrinsic :: ISO_C_BINDING
implicit none 

interface 
   function gzOpen( path, mode ) bind(C, NAME='gzopen') 
     use, intrinsic :: ISO_C_BINDING
     character(c_char) :: path(*),mode(*) 
     integer(c_int) :: gzOpen 
   end function
end interface


interface 
   function gzGets( file,buf,len ) bind(C, NAME='gzgets') 
     use, intrinsic :: ISO_C_BINDING
     character(c_char) :: buf(*)  
     type(c_ptr) :: gzGets 
     integer(c_int),value :: file,len
   end function
end interface

interface 
   function gzClose( file ) bind(C, NAME='gzclose') 
     use, intrinsic :: ISO_C_BINDING
     integer(c_int) :: gzClose
     integer(c_int),value :: file
   end function
end interface

end module
