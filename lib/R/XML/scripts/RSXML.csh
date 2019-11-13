if(`test -n "-lxml2 -lz -lm"`) then

if(${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:-lxml2 -lz -lm
else
   setenv LD_LIBRARY_PATH -lxml2 -lz -lm
endif

endif