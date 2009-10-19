OOPS	:= ../..
DEPTH	:= ..
LIBDIR	:= $(DEPTH)/libs

include $(OOPS)/mkincl

LIBS	+= \
	-L $(LIBDIR) \
	-lsmps++ -lsmps \
	$(NULL)

SOURCES := \
	test-Smps.cpp \
	test-SmpsCore.cpp \
	test-SmpsTree.cpp \
	test-Tokenizer.cpp \
	test-Utils.cpp \
	unit-tests.cpp \
	$(NULL)

HEADERS := \
	$(DEPTH)/Smps.h \
	unit-tests.h \
	$(NULL)

DEPENDS := \
	$(HEADERS) \
	$(LIBDIR)/libsmps++.a \
	$(NULL)

TESTEXE	:= xtests

# Allow the normal rules to run (so 'make' alone will compile the tests)
OUTPUT	:= $(TESTEXE)

include $(OOPS)/mkfhost/rules.mk

INCLUDE	+= -I $(DEPTH)
