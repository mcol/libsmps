OOPS	:= ../..
DEPTH	:= ..
LIBDIR	:= $(DEPTH)/libs

include $(OOPS)/mkincl

LIBS	+= \
	-L $(LIBDIR) \
	-L $(LPSOLVE) \
	-linterface -lsmps -lstoc \
	-llpsolve55 -ldl \
	$(NULL)

DEFINE  = \
	$(INCLUDE) \
	$(NULL)

SOURCES	:= \
	SmpsLpsolve.cpp \
	main.cpp \
	$(NULL)

HEADERS = \
	SmpsLpsolve.h \
	$(NULL)

DEPENDS = \
	$(HEADERS) \
	$(LIBDIR)/libsmps.a \
	$(LIBDIR)/libinterface.a \
	$(NULL)

OUTPUT  = x2lpsolve

include $(OOPS)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE = \
	-I $(OOPS)/include \
	-I $(DEPTH)/interface \
	-I $(DEPTH) \
	-I $(LPSOLVE) \
	$(NULL)
