DIR_LIB=lib
DIR_INC=lib
FTYPEBITS=32
GTYPEBITS=8
LIB_NAME=findr
LIB_NAMEFULL="Fast Inference of Networks from Directed Regulations"
LIB_FNAME=lib$(LIB_NAME).so
AUTHOR="Lingfei Wang"
AUTHOR_EMAIL="Lingfei.Wang.github@outlook.com"
URL_LIB="https://github.com/lingfeiwang/findr"
URL_BIN="https://github.com/lingfeiwang/findr-bin"
URL_PYTHON="https://github.com/lingfeiwang/findr-python"
URL_R="https://github.com/lingfeiwang/findr-R"
URL_DOC="https://github.com/lingfeiwang/findr/blob/master/doc.pdf"
VERSION1=0
VERSION2=1
VERSION3=0
LICENSE=AGPL-3
LICENSE_FULL="GNU Affero General Public License, Version 3"
LICENSE_URL="https://www.gnu.org/licenses/agpl-3.0"
#UNAME=$$(uname)
ifdef INCLUDE_MAKEFILE_BEFORE
#Input package info here
include $(INCLUDE_MAKEFILE_BEFORE)
endif

include Makefile.flags
ifndef PREFIX
PREFIX=/usr/local
endif
ifndef DIR_BUILD
DIR_BUILD=.
endif
ifndef DIR_SRC
DIR_SRC=.
endif
DIR_INSTALL_PREFIX=$(PREFIX)
DIR_BUILD_LIB_PREFIX=$(PREFIX)
DIR_INSTALL_BIN=$(DIR_INSTALL_PREFIX)/bin
DIR_BUILD_LIB=$(addsuffix /lib,$(DIR_BUILD_LIB_PREFIX)) ../$(DIR_LIB)/$(DIR_BUILD)
DIR_BUILD_INC=$(addsuffix /include/$(LIB_NAME),$(DIR_BUILD_LIB_PREFIX)) ../$(DIR_INC)/$(DIR_SRC)
#ifdef DIR_SRC_GSL
#LDFLAGS:=$(LDFLAGS) -Wl,--whole-archive $(DIR_SRC_GSL)/.libs/libgsl.a $(DIR_SRC_GSL)/cblas/.libs/libgslcblas.a -Wl,--no-whole-archive
#CFLAGS:=$(CFLAGS) -I$(DIR_SRC_GSL)
#else
#LDFLAGS:=$(LDFLAGS) -lgsl -lgslcblas
#endif
CC=gcc
LD=gcc
#INSTALL=install
COMMA=,
OPTFLAGS=-O3 -DNDEBUG=1 -DGSL_RANGE_CHECK_OFF=1 -DHAVE_INLINE=1
#CFLAGS:=$(CFLAGS) -fopenmp -ggdb -fPIC -Wall -Wextra -Wconversion -Wsign-conversion -Wundef -Wendif-labels -std=gnu99 -pedantic-errors $(addprefix -I,$(DIR_BUILD_INC)) $(OPTFLAGS)
#LDFLAGS:=$(LDFLAGS) -fopenmp -lc -lm $(addprefix -L,$(DIR_BUILD_LIB)) -l$(LIB_NAME) -Wl,-rpath="$$ORIGIN" $(addsuffix ",$(addprefix -Wl$(COMMA)-rpath=",$(DIR_BUILD_LIB)))

BIN_C=$(wildcard $(DIR_SRC)/*.c)
BIN_PRODUCT=$(addsuffix .o,$(basename $(BIN_C)))
BIN_DPRODUCT=$(DIR_BUILD)/$(LIB_NAME)
BIN_UNINSTALL=$(DIR_INSTALL_BIN)/$(LIB_NAME)

.PHONY: all clean distclean install uninstall

all: $(BIN_DPRODUCT)

$(DIR_BUILD):
	mkdir -p $@

$(BIN_PRODUCT):

$(BIN_DPRODUCT): $(BIN_PRODUCT) $(DIR_BUILD)
	$(LD) -o $@ $(BIN_PRODUCT) $(LDFLAGS)

clean:
	$(RM) $(BIN_PRODUCT)

distclean: clean
	$(RM) $(BIN_DPRODUCT) Makefile.flags

install: SHELL:=/bin/bash
install: all
	@umask 0022 && mkdir -p $(DIR_INSTALL_BIN) && \
	cp $(BIN_DPRODUCT) $(DIR_INSTALL_BIN)/ && \
	chmod 0755 $(DIR_INSTALL_BIN)/$(notdir $(BIN_DPRODUCT))

uninstall:
	$(RM) -R $(BIN_UNINSTALL)

TMP_FILE=.tmp
Makefile.flags:
	# Testing gcc
	$(CC) --version &> /dev/null || ( echo "GCC not found. Please download the latest GCC or specify its location in CC variable in Makefile."; exit 1; )
	gver=$$($(CC) --version | grep -o gcc) ; \
	if ! [ -n "$$gver" ]; then echo "Invalid GCC version. Please download the latest GCC."; exit 1; fi
	# Testing test method
	cflags="$(CFLAGS) $(CFLAGS_EXTRA) -DLIBINFONAME=$(LIB_NAME) -DLIBINFOVERSION=$(VERSION1).$(VERSION2).$(VERSION3) -fopenmp -ggdb -fPIC -Wall -Wextra -Wconversion -Wsign-conversion -Wundef -Wendif-labels -std=gnu99 -pedantic-errors $(addprefix -I,$(DIR_BUILD_INC)) $(OPTFLAGS)" ; \
	ldflags="$(LDFLAGS) -fopenmp -lc -lm"; \
	$(LD) $$ldflags -lc --shared -o $(TMP_FILE) &> /dev/null || \
	( echo "Linking with default flags failed."; exit 1; ) ; \
	echo "Testing -Wl,--no-as-needed" ; \
	$(LD) -Wl,--no-as-needed $$ldflags --shared -o $(TMP_FILE) &> /dev/null && \
	ldflags="-Wl,--no-as-needed $$ldflags"; \
	echo "Testing -Wl,-rpath" ; \
	$(LD) -Wl,-rpath="$$$$ORIGIN" $$ldflags --shared -o $(TMP_FILE) &> /dev/null && \
	ldflags="-Wl,-rpath=\""'$$$$'"ORIGIN\" $(addsuffix \",$(addprefix -Wl$(COMMA)-rpath=\",$(DIR_BUILD_LIB))) $$ldflags $(addprefix -L,$(DIR_BUILD_LIB))"; \
	echo "CFLAGS=$$cflags" > $@ && \
	echo "LDFLAGS=$$ldflags -L. -l$(LIB_NAME) -lgsl -lgslcblas" >> $@
	$(RM) $(TMP_FILE)

ifdef INCLUDE_MAKEFILE_AFTER
include $(INCLUDE_MAKEFILE_AFTER)
endif












