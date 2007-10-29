DEPTH   := ..

include $(DEPTH)/mkincl

CC = g++

OBJECTS := \
	SmpsCore.o \
	SmpsTree.o \
	$(NULL)

HEADERS := \
	smps.h \
	$(NULL)

DEPENDS := \
	$(HEADERS) \
	$(NULL)

LIBRARY := libsmps.a

include $(DEPTH)/mkfhost/rules.mk
