DEPTH   = ../..
INTERFACE = ../interface

include $(DEPTH)/mkincl

CPLEX   = /usr/local/ilog/cplex10

LIBDIRS += \
	-L $(CPLEX)/lib/static_pic \
	-L $(INTERFACE) \
	-L .. \
	$(NULL)

LIBS    += $(LIBDIRS) \
	-linterface -lsmps -lmps \
	-lilocplex -lcplex -lpthread \
	$(NULL)

DEFINE  = \
	$(INCLUDE) \
	$(NULL)

OBJECTS = \
	SmpsCplex.o \
	main.o \
	$(NULL)

HEADERS = \
	SmpsCplex.h \
	$(NULL)

DEPENDS = \
	$(HEADERS) \
	$(INTERFACE)/libinterface.a \
	$(NULL)
#	$(DEPTH)/libs/libinterface.a \

OUTPUT  = x2cplex

include $(DEPTH)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE = \
	-I $(CPLEX)/include \
	-I $(DEPTH)/include \
	-I ../interface \
	-I .. \
	$(NULL)
