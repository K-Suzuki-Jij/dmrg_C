PROGNAME := 1D_KLM_TVF_EXACT.out
SRCDIR   := src
OUTDIR   := build
INCDIR   := ../../../include
LIBDIR   := ../../../lib
TARGET   := $(OUTDIR)/$(PROGNAME)
SRCS     := $(wildcard $(SRCDIR)/*.c)
OBJS     := $(addprefix $(OUTDIR)/,$(patsubst %.c,%.o,$(SRCS)))

# Mac OS clang
CC = gcc
CFLAGS = -Ofast -Xpreprocessor -fopenmp -I $(INCDIR) -I include
LFLAGS = -ldmrg -lonsite_op -lsml -llapack -lomp -L $(LIBDIR)

# Linux
#CC = gcc
#CFLAGS = -Ofast -fopenmp -I $(INCDIR) -I include
#LFLAGS = -ldmrg -lonsite_op -lsml -llapack -lm -L $(LIBDIR)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(OUTDIR)/%.o:%.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf $(OUTDIR)
