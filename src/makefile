CC      =  gcc
CFLAGS  =  -g -Wall  -D DEBUG #-O2 
LDFLAGS = -lz -lm

OBJS = gfa.o opt.o paf.o sdict.o eg.o 
PROG = calcuts purge_dups pbcstat split_fa get_seqs  ngscstat

.SUFFIXS:.c .o

all:$(PROG)

calcuts: calcuts.o 
	$(CC)  $(CFLAGS) $^ -o $@ $(LDFLAGS) 
purge_dups:opt.o sdict.o asset.o cov.o paf.o purge_dups.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
pbcstat: asset.o pbcstat.o paf.o sdict.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
ngscstat: asset.o ngscstat.o bamlite.o sdict.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
split_fa: split_fa.o  
	$(CC) $(CFLAGS) $^ -o $@  $(LDFLAGS) 
get_seqs: get_seqs.o sdict.o build_graph.o cdict.o bed.o graph.o
	$(CC) $(CFLAGS) $^ -o $@  $(LDFLAGS) 

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

	
