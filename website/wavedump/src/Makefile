release_hdr := $(shell sh -c './mkreleasehdr.sh')

OPTIMIZATION?=-O2
CFLAGS=$(OPTIMIZATION) -Wall -g -std=gnu99
CFLAGS+= -fPIC -DLINUX -std=gnu99

LDLIBS=-lCAENComm -lCAENDigitizer -lm -lpthread -lhdf5

PREFIX?=/usr/local
INSTALL_BIN=$(PREFIX)/bin
INSTALL=install

all: wavedump

.PHONY: all

Makefile.dep:
	$(CC) $(CFLAGS) -MM *.c > Makefile.dep

-include Makefile.dep

wavedump: wavedump.o fft.o flash.o  keyb.o  spi.o WDconfig.o  X742CorrectionRoutines.o release.o

install: all
	@mkdir -p $(INSTALL_BIN)
	$(INSTALL) wavedump $(INSTALL_BIN)

clean:
	rm -f wavedump *.o
