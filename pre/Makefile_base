#!/bin/sh
#----------------------------------------------------------
AUTHOR="Henrique de Andrade"
#------------------Variaveis Globais-----------------------
name=prepar
MPI=no
PARMETIS=no
METIS=
DEBUG=no
FC=gfortran
CO=gfortran
PATH_MOD=./src
PATH_LIB_METIS=./lib
PATH_INCLUDE=./include
# ----------------------------------------------------------
# Codigo com opcao de otimizacao codigo compilado intel
# -O2 e incompativel com os compiladores da gnu
# nao o contrario no e verdade
# logo usaremos para questao de desempenho libs compiladas
# com intel com o compilador intel
# e libs compilados com gnu com o campilar dor gnu
#make
mod_mef=src/Malloc.f

src_mef = src/Adjacency.f\
src/Filenames.f\
src/Main.f\
src/Mpi_comunicate.f\
src/Parallel_imp.f\
src/Partdmesh.f\
src/Pdmetis.f\
src/Propnode.f\
src/Rdata.f\
src/Read_mef.f\
src/Read_metis.f\
src/Time.f\
src/Vtk.f\
src/Write_mef.f\
src/Write_par.f\
src/Write_vtk.f
#src_mef=$(wilcard ./src/*.f)
#-------------------Flags necessarios--------------------------------
NFLAGS=-I$(PATH_INCLUDE) -L$(PATH_LIB_METIS)
LDFLAGS=
#--------------------compiladores------------------------------------
# intel ifort
ifeq ($(CO),ifort)
  OFLAGS  += -fc=ifort -module $(PATH_MOD) -fpp
  LDFLAGS +=
  ifeq ($(MPI),yes)
    LDFLAGS += -lmpich
    OFLAGS  += -DMPI
  endif
  ifeq ($(PARMETIS),yes)
    LDFLAGS += -lparmetis-x86
    OFLAGS  += -DPARMETIS
  endif
#---------------Biblioteca metis 32 ou 64-----------------------------
  ifeq ($(METIS),64)
    OFLAGS  += -DMETIS64
    LDFLAGS += -lmetis-x86_64
  else
    OFLAGS  +=
    LDFLAGS += -lmetis-x86
  endif
endif
# gnu gcc
ifeq ($(CO),gfortran)
  # OFLAGS  += -fc=gfortran -J$(PATH_MOD) -cpp
  # OFLAGS  += -J$(PATH_MOD) -fallow-argument-mismatch -cpp# Fedora gfortran 12.1.1
  OFLAGS  += -J$(PATH_MOD) -cpp
  LDFLAGS +=
  ifeq ($(MPI),yes)
    LDFLAGS += -lmpich
    OFLAGS  += -DMPI
  endif
  ifeq ($(PARMETIS),yes)
    LDFLAGS += -lparmetis-x86
    OFLAGS  += -DPARMETIS
  endif
#---------------Biblioteca metis 32 ou 64-----------------------------
  ifeq ($(METIS),64)
    OFLAGS  += -DMETIS64
    LDFLAGS += -lmetis-x86_64
  else
    OFLAGS  += -DMETIS5
    LDFLAGS += -lmetis
  endif
endif
#--------------------------------------------------------------------
#---------------------------Debug------------------------------------
ifeq ($(DEBUG),yes)
  OFLAGS += -g
else
  OFLAGS += -O2
endif
#--------------------------------------------------------------------
FFLAGS= $(NFLAGS) $(OFLAGS)

.SUFFIXES:
.SUFFIXES: .for .f .h .fi .o
mod_src  = $(mod_mef:%.f=%.o)
objs_src = $(mod_src) $(src_mef:%.f=%.o)

#par:
#	mpif77 $(FFLAGS) $(LDFFLAGS) -names as_is $(src_par)


#build:	$(objs_src)
build:	$(objs_src) $(mod_src)
	@mkdir -p bin
	$(FC) $(objs_src) $(FFLAGS) -o bin/$(name)  $(LDFLAGS)

tags:
	ctags -R src/*.f include/*.fi

.PHONY: cleantags
cleantags:
	@rm -v tags

.PHONY: clean
clean:
	@rm -v src/*.o
	@rm -v src/*.mod
	@rm -v bin/$(name)

.PHONY: cleanall
cleanall:
	@rm -v tags
	@rm -v src/*.o
	@rm -v bin/$(name)

.PHONY: cleanmod
cleanmod:
	@rm -v src/*.mod

.PHONY: help
help:
	@echo "Autor :$(AUTHOR)                              "
	@echo "Makefile para prepar para sitemas linux.      "
	@echo -e "\E[7;32mOpcoes:\E[1;0m                      "
	@echo "build         - compila o prepar              "
	@echo "tags          - gera os tags                  "
	@echo "cleantags     - limpa os tags                 "
	@echo "clean         - limpa os obj, bin e mod       "
	@echo "cleaall       - limpa tudo obj,bin,mod e tags "
	@echo "cleanmod      - limpa os modulos              "

# DO NOT DELETE
