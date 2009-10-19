OOPS	:= ..
DEPTH	:= .
LIBDIR	:= $(DEPTH)/libs

include $(OOPS)/mkincl

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
	$(LIBDIR)/libsmps.a \
	$(NULL)

APPS	:= \
	interface \
	cplex-driver \
	hopdm-driver \
	lpsolve-driver \
	oops-driver \
	$(NULL)

LIBRARY := $(LIBDIR)/libsmps++.a

include $(OOPS)/mkfhost/rules.mk
