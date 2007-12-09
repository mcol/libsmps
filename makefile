DEPTH   := ..

include $(DEPTH)/mkincl

CC = g++

SOURCES	:= \
	Node.cpp \
	Smps.cpp \
	SmpsCore.cpp \
	SmpsTree.cpp \
	Tokenizer.cpp \
	Utils.cpp \
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
