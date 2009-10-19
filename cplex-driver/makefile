OOPS	:= ../..
DEPTH	:= ..
LIBDIR	:= $(DEPTH)/libs

include $(OOPS)/mkincl

CPLEX   = /usr/local/ilog/cplex

LIBS	+= \
	-L $(LIBDIR) \
	-L $(CPLEX)/lib/x86_debian4.0_4.1/static_pic \
	-linterface -lsmps++ -lsmps \
	-lilocplex -lcplex -lpthread \
	$(NULL)

DEFINE  = \
	$(INCLUDE) \
	$(NULL)

SOURCES	:= \
	SmpsCplex.cpp \
	main.cpp \
	$(NULL)

HEADERS = \
	SmpsCplex.h \
	$(NULL)

DEPENDS = \
	$(HEADERS) \
	$(LIBDIR)/libsmps++.a \
	$(LIBDIR)/libinterface.a \
	$(NULL)

OUTPUT  = x2cplex

include $(OOPS)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE = \
	-I $(CPLEX)/include \
	-I $(OOPS)/include \
	-I $(DEPTH)/interface \
	-I $(DEPTH) \
	$(NULL)
