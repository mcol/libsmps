DEPTH   := ../..

include $(DEPTH)/mkfhost/oopshome
include $(DEPTH)/mkincl

CC = g++

LIBSMPS	:= ../libsmps.a

LIBS	+= $(LIBSMPS)

TESTOBJS := \
	test-SmpsCore.o \
	test-SmpsTree.o \
	unit-tests.o \
	$(NULL)

HEADERS := \
	../smps.h \
	unit-tests.h \
	$(NULL)

DEPENDS := \
	$(HEADERS) \
	$(LIBSMPS) \
	$(NULL)

TESTEXE	:= xtests

# Allow the normal rules to run (so 'make' alone will compile the tests)
OBJECTS	:= $(TESTOBJS)
OUTPUT	:= $(TESTEXE)

include $(DEPTH)/mkfhost/rules.mk

INCLUDE	+= -I ../
