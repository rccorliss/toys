# toys
various small toy programs that wanted version control.
In the current configuration, compiling requires:
autoreconf -fvi
make
#look at the failed linker step that double-refs the _Dict.o.   Remove the second ref and run that command by hand
#move the _Dict...pcm file into .libs/
