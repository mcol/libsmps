OOPS	:= ../..
DEPTH	:= ..

include $(OOPS)/mkincl

LIBS	+= \
	$(LIBOOPS) \
	$(LAPACKLIB) \
	$(METISLIB) \
	-L $(DEPTH) \
	-L ../interface -linterface \
	-lsmps -lmps \
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
	../libsmps.a \
	../interface/libinterface.a \
	$(NULL)

OUTPUT	:= x2oops

include $(OOPS)/mkfhost/rules.mk

# has to be below rules.mk because it adds!
INCLUDE += \
        -I $(DEPTH) \
	-I $(DEPTH)/interface \
        $(NULL)
