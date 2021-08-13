SMLDIR := sml
INCDIR := ../include

CC = gcc
CFLAGS = -Wall -Wextra -Ofast -mtune=native -march=native -Xpreprocessor -fopenmp -I $(INCDIR)
LFLAGS = -llapack -lomp

.PHONY: all sml

all: sml

sml  : libsml.a

###SML
libsml.a: $(addprefix $(SMLDIR)/,$(notdir $(patsubst %.c,%.o,$(wildcard ../$(SMLDIR)/*.c))))
	ar r $@ $^

$(SMLDIR)/%.o: ../$(SMLDIR)/%.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $(SMLDIR)/$(@F) -c ../$(SMLDIR)/$(<F)

clean:
	rm -rf $(SMLDIR)
