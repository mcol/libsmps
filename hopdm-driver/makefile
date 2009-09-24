OOPS	:= ../..
DEPTH	:= ..
LIBDIR	:= $(DEPTH)/libs

include $(OOPS)/mkincl

LIBS    += \
	-L $(LIBDIR) \
	-linterface -lsmps -lstoc \
	-lhopdminterface -lhopdm \
	$(NULL)

DEFINE  = \
	$(INCLUDE) \
	$(NULL)

SOURCES	:= \
	main.cpp \
	SmpsHopdm.cpp \
	$(NULL)

HEADERS = \
	SmpsHopdm.h \
	$(NULL)

DEPENDS = \
	$(HEADERS) \
	$(LIBDIR)/libsmps.a \
	$(LIBDIR)/libinterface.a \
	$(NULL)

OUTPUT  = x2hopdm

include $(OOPS)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE = \
	-I $(OOPS)/include \
	-I $(DEPTH)/interface \
	-I $(DEPTH) \
	$(NULL)
