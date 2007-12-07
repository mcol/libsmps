DEPTH   := ../..

include $(DEPTH)/mkincl

CC = g++

OBJECTS := \
	generateSmps.o \
	main.o \
	SmpsOops.o \
	$(NULL)

HEADERS := \
	SmpsOops.h \
	../Smps.h \
	$(NULL)

DEPENDS := \
	$(HEADERS) \
	$(LIBOOPS) \
	../libsmps.a \
	../interface/libinterface.a \
	$(NULL)

LIBS	+= \
	$(LIBOOPS) \
	$(LAPACKLIB) \
	$(METISLIB) \
	-L .. -lsmps -lmps \
	-L ../interface -linterface \
	$(NULL)

OUTPUT	:= x2oops

include $(DEPTH)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE += \
        -I .. \
	-I ../interface \
        $(NULL)
