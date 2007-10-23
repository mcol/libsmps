DEPTH   := ..

include $(DEPTH)/mkfhost/oopshome
include $(DEPTH)/mkincl

CC = g++

OBJECTS := \
	SmpsCore.o \
	SmpsTree.o \
	$(NULL)

TESTS	:= \
	SmpsCore.o \
	SmpsTree.o \
	unit-tests/test-SmpsCore.o \
	unit-tests/test-SmpsTree.o \
	unit-tests/unit-tests.o \
	$(NULL)

HEADERS := \
	smps.h \
	$(NULL)

DEPENDS := \
	$(HEADERS) \
	$(NULL)

LIBRARY := libsmps.a

include $(DEPTH)/mkfhost/rules.mk
