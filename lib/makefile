SMLDIR   := sml
DMRGDIR  := dmrg
ONDIR    := onsite_op
EXACTDIR := exact
INCDIR   := ../include

# Mac OS clang
CC = gcc
CFLAGS = -Wall -Wextra -Ofast -Xpreprocessor -fopenmp -I $(INCDIR)

# Linux
#CC = gcc
#CFLAGS = -Wall -Wextra -Ofast -fopenmp -I $(INCDIR)

.PHONY: all sml dmrg

all: sml dmrg op ed

sml : libsml.a
dmrg: libdmrg.a
op  : libonsite_op.a
ed  : libexact.a

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
	
###ONSITE_OP
libonsite_op.a: $(addprefix $(ONDIR)/,$(notdir $(patsubst %.c,%.o,$(wildcard ../$(ONDIR)/*.c))))
	ar r $@ $^

$(ONDIR)/%.o: ../$(ONDIR)/%.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(ONDIR)/$(@F) -c ../$(ONDIR)/$(<F)
	
###EXACT
libexact.a: $(addprefix $(EXACTDIR)/,$(notdir $(patsubst %.c,%.o,$(wildcard ../$(EXACTDIR)/*.c))))
	ar r $@ $^

$(EXACTDIR)/%.o: ../$(EXACTDIR)/%.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(EXACTDIR)/$(@F) -c ../$(EXACTDIR)/$(<F)

clean:
	rm -rf $(SMLDIR) $(DMRGDIR) $(ONDIR) $(EXACTDIR)
