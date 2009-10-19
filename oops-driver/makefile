OOPS	:= ../..
DEPTH	:= ..
LIBDIR	:= $(DEPTH)/libs

include $(OOPS)/mkincl

LIBS	+= \
	$(LIBOOPS) \
	$(LAPACKLIB) \
	$(METISLIB) \
	-L $(LIBDIR) \
	-linterface -lsmps++ -lsmps \
	$(NULL)

SOURCES	:= \
	generateSmps.cpp \
	main.cpp \
	SmpsOops.cpp \
	$(NULL)

HEADERS := \
	SmpsOops.h \
	$(DEPTH)/Smps.h \
	$(NULL)

DEPENDS := \
	$(HEADERS) \
	$(LIBOOPS) \
	$(LIBDIR)/libsmps++.a \
	$(LIBDIR)/libinterface.a \
	$(NULL)

OUTPUT	:= x2oops

include $(OOPS)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE += \
        -I $(DEPTH) \
	-I $(DEPTH)/interface \
        $(NULL)
