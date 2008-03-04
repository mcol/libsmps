OOPS	:= ../..
DEPTH	:= ..
LIBDIR	:= $(DEPTH)/libs

include $(OOPS)/mkincl

DEFINE  = \
        -DUNDERSCORE_F \
	$(INCLUDE) \
	$(NULL)

SOURCES	:= \
	interface.cpp \
	options.cpp \
	$(NULL)

HEADERS = \
	interface.h \
	options.h \
	$(NULL)

DEPENDS = \
	$(HEADERS) \
	$(LIBDIR)/libsmps.a \
	$(NULL)

LIBRARY = $(LIBDIR)/libinterface.a

include $(OOPS)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE += \
	-I $(OOPS)/include \
	-I $(DEPTH) \
        $(NULL)
