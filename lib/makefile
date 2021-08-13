SMLDIR := sml
DMRGDIR := dmrg
INCDIR := ../include

CC = gcc
CFLAGS = -Wall -Wextra -Ofast -mtune=native -march=native -Xpreprocessor -fopenmp -I $(INCDIR)
LFLAGS = -llapack -lomp

.PHONY: all sml dmrg

all: sml dmrg

sml : libsml.a
dmrg: libdmrg.a

###SML
libsml.a: $(addprefix $(SMLDIR)/,$(notdir $(patsubst %.c,%.o,$(wildcard ../$(SMLDIR)/*.c))))
	ar r $@ $^

$(SMLDIR)/%.o: ../$(SMLDIR)/%.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(SMLDIR)/$(@F) -c ../$(SMLDIR)/$(<F)
	
###DMRG
libdmrg.a: $(addprefix $(DMRGDIR)/,$(notdir $(patsubst %.c,%.o,$(wildcard ../$(DMRGDIR)/*.c))))
	ar r $@ $^

$(DMRGDIR)/%.o: ../$(DMRGDIR)/%.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(DMRGDIR)/$(@F) -c ../$(DMRGDIR)/$(<F)

clean:
	rm -rf $(SMLDIR) $(DMRGDIR)
