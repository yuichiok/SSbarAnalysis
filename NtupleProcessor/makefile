# Make file for ZCAnalysis/NtupleProcessor package

CC=g++
#CC=g++-4.9
#CFLAGS=-std=c++14 -c -Wall
CFLAGS=-c -Wall -fdiagnostics-color
ROOTFLAGS=`root-config --cflags --glibs --libs`
OBJDIR=NtupleProcessor/obj
OBJS=$(OBJDIR)/main.o $(OBJDIR)/NtupleProcessor.o
HDRDIR=NtupleProcessor/include
SRCDIR=NtupleProcessor/src
LIBDIR=SSbarLibrary


NtupleProcessor : createObjDir NtupleProcessor.exe

createObjDir :
	mkdir -p NtupleProcessor/obj

NtupleProcessor.exe : $(OBJS)
	$(CC) $(OBJS) -o NtupleProcessor.exe $(ROOTFLAGS)

$(OBJDIR)/main.o : NtupleProcessor/main.cc $(HDRDIR)/NtupleProcessor.hh
	$(CC) $(CFLAGS) NtupleProcessor/main.cc -o $(OBJDIR)/main.o $(ROOTFLAGS)

$(OBJDIR)/NtupleProcessor.o : $(SRCDIR)/NtupleProcessor.cc $(HDRDIR)/NtupleProcessor.hh
	$(CC) $(CFLAGS) $(SRCDIR)/NtupleProcessor.cc -o $(OBJDIR)/NtupleProcessor.o $(ROOTFLAGS)

$(OBJDIR)/TreeIterator.o : $(SRCDIR)/TreeIterator.cc $(HDRDIR)/TreeIterator.hh $(HDRDIR)/EventAnalyzer.hh
	$(CC) $(CFLAGS) $(SRCDIR)/TreeIterator.cc -o $(OBJDIR)/TreeIterator.o $(ROOTFLAGS)

$(OBJDIR)/EventAnalyzer.o : $(SRCDIR)/EventAnalyzer.cc $(HDRDIR)/EventAnalyzer.hh
	$(CC) $(CFLAGS) $(SRCDIR)/EventAnalyzer.cc -o $(OBJDIR)/EventAnalyzer.o $(ROOTFLAGS)

clean :
	rm NtupleProcessor/obj/*.o NtupleProcessor.exe