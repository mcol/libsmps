DEPTH   := ..

include $(DEPTH)/mkincl

CC = g++

OBJECTS := \
	Node.o \
	Smps.o \
	SmpsCore.o \
	SmpsTree.o \
	Tokenizer.o \
	Utils.o \
	$(NULL)

HEADERS := \
	Node.h \
	Smps.h \
	Tokenizer.h \
	Utils.h \
	$(NULL)

DEPENDS := \
	$(HEADERS) \
	$(NULL)

LIBRARY := libsmps.a

include $(DEPTH)/mkfhost/rules.mk
