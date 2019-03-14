CC=g++
CCFLAGS=-Wall `root-config --cflags --libs` -fopenmp -g
EXECS=Module_SNBurst.exe Plotting_SNBurst.exe MakeSNTheoryDistributions.exe CNN_burst_trigger.exe
# Plotting_SNBurst_Efficiency.exe Module_BetterSNBurst.exe
NERROR=100
OPTERROR= -fmax-errors=$(NERROR)

%.exe: %.C
	$(CC) -o $@ $^ $(CCFLAGS) $(OPTERROR)

all: $(EXECS)

clean:
	rm -rf *.exe

cleanpdf:
	rm -rf *.pdf

