DEPTH   := ..

include $(DEPTH)/mkfhost/oopshome
include $(DEPTH)/mkincl

CC = g++

OBJECTS := \
	SmpsTree.o \
	$(NULL)

TESTS	:= \
	SmpsTree.o \
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
