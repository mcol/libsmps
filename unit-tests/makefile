OOPS	:= ../..
DEPTH	:= ..
LIBDIR	:= $(DEPTH)/libs

include $(OOPS)/mkincl

LIBS	+= \
	-L $(LIBDIR) \
	-lsmps -lstoc \
	$(NULL)

TESTOBJS := \
	test-Smps.o \
	test-SmpsCore.o \
	test-SmpsTree.o \
	test-Tokenizer.o \
	test-Utils.o \
	unit-tests.o \
	$(NULL)

HEADERS := \
	$(DEPTH)/Smps.h \
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

include $(OOPS)/mkfhost/rules.mk

INCLUDE	+= -I $(DEPTH)
