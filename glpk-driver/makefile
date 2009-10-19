OOPS	:= ../..
DEPTH	:= ..
LIBDIR	:= $(DEPTH)/libs

GLPK	:= /home/mcolombo/local/

include $(OOPS)/mkincl

LIBS	+= \
	-L $(LIBDIR) \
	-L $(GLPK)/lib \
	-lsmps++ -lsmps \
	-lglpk \
	$(NULL)

DEFINE  = \
	$(INCLUDE) \
	$(NULL)

SOURCES	:= \
	SmpsGlpk.cpp \
	main.cpp \
	$(NULL)

HEADERS = \
	SmpsGlpk.h \
	$(NULL)

DEPENDS = \
	$(HEADERS) \
	$(LIBDIR)/libsmps++.a \
	$(NULL)

OUTPUT  = x2glpk

include $(OOPS)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE = \
	-I $(OOPS)/include \
	-I $(DEPTH)/interface \
	-I $(DEPTH) \
	-I $(GLPK)/include \
	$(NULL)
