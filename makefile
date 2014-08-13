CC=g++

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
ROOT_D := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

ifeq ($(CFG),debug)
	CXXFLAGS=-Wall -O2 -g -DDEBUG -DROOT_DIR=\"$(ROOT_D)\"
	LDFLAGS=-Wall -g -DDEBUG
	EXEC=$(BINPATH)/mindthegapDebug
else
	CXXFLAGS=-Wall -O2 -DROOT_DIR=\"$(ROOT_D)\"
	LDFLAGS=-Wall
	EXEC=$(BINPATH)/mindthegap
endif

INCPATH=$(ROOT_D)/inc
SRCPATH=$(ROOT_D)/src
OBJPATH=$(ROOT_D)/obj
BINPATH=$(ROOT_D)/bin

SRC=$(SRCPATH)/util.cpp \
	$(SRCPATH)/ctmf.cpp \
	$(SRCPATH)/medianBlur.cpp \
	$(SRCPATH)/posterize.cpp \
    $(SRCPATH)/lodepng.cpp \
    $(SRCPATH)/fitcurve.cpp \
    $(SRCPATH)/GGVecLib.cpp \
    $(SRCPATH)/svg.cpp \
    $(SRCPATH)/graph.cpp \
    $(SRCPATH)/codeImage.cpp \
    $(SRCPATH)/bitmap.cpp \
    $(SRCPATH)/main.cpp \
    $(SRCPATH)/info.cpp
OBJ=$(OBJPATH)/util.o \
	$(OBJPATH)/ctmf.o \
	$(OBJPATH)/medianBlur.o \
	$(OBJPATH)/posterize.o \
	$(OBJPATH)/lodepng.o \
	$(OBJPATH)/fitcurve.o \
	$(OBJPATH)/GGVecLib.o \
	$(OBJPATH)/svg.o \
	$(OBJPATH)/graph.o \
    $(OBJPATH)/codeImage.o \
    $(OBJPATH)/bitmap.o \
    $(OBJPATH)/main.o \
    $(OBJPATH)/info.o

INCLUDES=-I $(INCPATH)

default: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJPATH)/%.o: $(SRCPATH)/%.cpp $(INCPATH)/%.h
	$(CC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean cleanall

clean:
	rm -f $(OBJPATH)/*.o

cleanall: clean
	rm -f $(BINPATH)/*