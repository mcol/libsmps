OOPS	:= ../..
DEPTH	:= ..

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
	$(NULL)

LIBRARY = libinterface.a

include $(OOPS)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE += \
	-I $(OOPS)/include \
	-I $(DEPTH) \
        $(NULL)
