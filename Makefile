#Arquivo de entrada
IN = cholesky.c

#Executável
OUT = cholesky

#Flag para warnings
CFLAGS = -Wall -pg

#Flags para otimização de espaço
OFLAGS = -Os -fconserve-stack -ffunction-sections -fdata-sections
OFLAGS += -fno-force-addr -fno-inline-functions -fbranch-probabilities -finline-limit=1
OFLAGS += -fno-schedule-insns -fno-optimize-sibling-calls -fno-if-conversion -fprofile-arcs

#Compilador
CC = gcc

#Compila sem nenhum tipo de otimização
all:
	$(CC) $(CFLAGS) $(IN) -o $(OUT)

#Compila com as flags de otimização de espaço
opt:
	$(CC) $(CFLAGS) $(OFLAGS) $(IN) -o $(OUT)

#Gera o código assembly
assembly:
	$(CC) $(CFLAGS) $(OFLAGS) -S $(IN)

gprof:
	./$(OUT)
	$(shell gprof $(OUT) gmon.out > analysis.txt)
	@echo ""
	@echo "Profiling concluído. Verifique \"analysis.txt\""
	@echo ""
	@echo "Para saber mais sobre o GNU GCC Profiling Tool, visite:"
	@echo "http://www.thegeekstuff.com/2012/08/gprof-tutorial"


clean:
	rm -rf *.s *.out $(OUT)
