OOPS	:= ../..
DEPTH	:= ..
LIBDIR	:= $(DEPTH)/libs

include $(OOPS)/mkincl

CPLEX   = /usr/local/ilog/cplex10

LIBS	+= \
	-L $(LIBDIR) \
	-L $(CPLEX)/lib/static_pic \
	-linterface -lsmps -lstoc \
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
	$(LIBDIR)/libsmps.a \
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
