PROGNAME := 1D_KLM_VF_EXACT.out
SRCDIR   := src
OUTDIR   := build
INCDIR   := ../../../include
LIBDIR   := ../../../lib
TARGET   := $(OUTDIR)/$(PROGNAME)
SRCS     := $(wildcard $(SRCDIR)/*.c)
OBJS     := $(addprefix $(OUTDIR)/,$(patsubst %.c,%.o,$(SRCS)))

#CC = icc
#CFLAGS = -Wall -g -qopenmp -I $(INCDIR)
#LFLAGS = -mkl -lsmlicc -lonsiteicc -ldmrgicc -lexacticc

#CC = icc
#CFLAGS = -Wall -fast -parallel -qopenmp -I $(INCDIR)
#LFLAGS = -mkl -lonsiteicc -ldmrgicc -lexacticc -lsmlicc

CC = gcc
CFLAGS = -Ofast -Xpreprocessor -fopenmp -I $(INCDIR) -I include
LFLAGS = -llapack -lomp -lsml -lexact -lonsite_op -L $(LIBDIR)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(OUTDIR)/%.o:%.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf $(OUTDIR)
