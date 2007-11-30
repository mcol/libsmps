DEPTH   := ../..

include $(DEPTH)/mkincl

CC = g++

OBJECTS := \
	generateSmps.o \
	main.o \
	options.o \
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
	$(NULL)

LIBS	+= \
	$(LIBOOPS) \
	$(LAPACKLIB) \
	$(METISLIB) \
	-L .. -lsmps -lmps \
	$(NULL)

OUTPUT	:= x2oops

include $(DEPTH)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE += \
        -I .. \
        $(NULL)
