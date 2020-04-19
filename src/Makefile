include ../config/makeconfig

MAKE = make

all: exe

exe: lib
	cd ./main; $(MAKE); cd ..
#exe: lib
#	cd ./main; $(MAKE) -f Makefile.proj; cd ..

#lib: module.lib control.lib io.lib mpi.lib bas.lib fedvr.lib dvr.lib sph.lib tdcc.lib ormas.lib wfn.lib hprod.lib prop.lib field.lib gbas.lib pes.lib lapack.lib util.lib futil.lib
lib: module.lib control.lib io.lib mpi.lib bas.lib fedvr.lib dvr.lib sph.lib tdcc.lib ormas.lib wfn.lib hprod.lib prop.lib field.lib gbas.lib pes.lib surff.lib lapack.lib util.lib futil.lib

module.lib:
	cd $(@:.lib=); $(MAKE); cd ..
control.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
io.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
mpi.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
bas.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
fedvr.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
dvr.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
sph.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
tdcc.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
ormas.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
wfn.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
hprod.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
prop.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
field.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
gbas.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
pes.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
surff.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
lapack.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
util.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
futil.lib: 
	cd $(@:.lib=); $(MAKE); cd ..

install: 
	mkdir -p $(TD1C_TOP)/bin
	ar rcs libtd1c.$(TARGET).a */*.o
	mv libtd1c.$(TARGET).a $(TD1C_TOP)/bin
	cp main/td1c.exe $(TD1C_TOP)/bin/td1c.$(TARGET).exe

clean:
	rm -f */*.o */*.mod
fullclean:
	rm -f */*.o */*.mod */*__genmod.* *~ */*~
