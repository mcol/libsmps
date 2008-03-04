OOPS	:= ../..
DEPTH	:= ..
INTERFACE = $(DEPTH)/interface

include $(OOPS)/mkincl

CPLEX   = /usr/local/ilog/cplex10

LIBS	+= \
	-L $(CPLEX)/lib/static_pic \
	-L $(INTERFACE) \
	-L $(DEPTH) \
	-linterface -lsmps -lmps \
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
	$(INTERFACE)/libinterface.a \
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
